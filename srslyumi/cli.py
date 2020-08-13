"""
Process SRSLY UMI from BCL to FASTQ
"""

import argparse
import sys
import os
import os.path
import stat
import subprocess
import xml.etree.ElementTree as ET

DEFAULT_BCL_SUBDIR = "Data/Intensities/BaseCalls"
DEFAULT_SAMPLE_SHEET = "SampleSheet.csv"
DEFAULT_RUN_INFO = "RunInfo.xml"
DEFALUT_BCL_CMD_FN = "runbcl2fastq.sh"
BCL2FASTQ_TEMPLATE = """#!/bin/sh
bcl2fastq --input-dir {BCL_DIR} \\
    --runfolder-dir {OUT_DIR} \\
    --output-dir {OUT_DIR}
"""


def read_run_info(fn, umi_bp, index_bp):
    tree = ET.parse(fn)
    root = tree.getroot()
    reads = tree.find("./Run/Reads")
    if reads is None:
        raise ValueError("No <RunInfo><Run><Reads> element")
    read_data = {}
    for read in reads:
        read_data[read.attrib["Number"]] = read.attrib
    if len(read_data) != 4 or not all([str(x) in read_data for x in range(1, 5)]):
        raise ValueError("SRSLY setup requires 4 index reads")
    i7len = str(umi_bp + index_bp)
    if read_data["2"]["NumCycles"] != i7len:
        raise ValueError("i7 Index Read should have {} cycles".format(i7len))
    if read_data["2"]["IsIndexedRead"] != "Y":
        raise ValueError("RunInfo.xml Read 2 should be an index read")
    if read_data["3"]["IsIndexedRead"] != "Y":
        raise ValueError("RunInfo.xml Read 3 should be an index read")

    reads.clear()
    ET.SubElement(reads, "Read", attrib=read_data["1"])

    r2attrib = dict(read_data["2"])
    r2attrib["NumCycles"] = str(index_bp)
    ET.SubElement(reads, "Read", attrib=r2attrib)

    r3attrib = {"NumCycles": str(umi_bp), "Number": "3", "IsIndexedRead": "N"}
    ET.SubElement(reads, "Read", attrib=r3attrib)

    r4attrib = dict(read_data["3"])
    r4attrib["Number"] = "4"
    ET.SubElement(reads, "Read", attrib=r4attrib)

    r5attrib = dict(read_data["4"])
    r5attrib["Number"] = "5"
    ET.SubElement(reads, "Read", attrib=r5attrib)
    return tree, read_data


def read_csv(fn, field_delim=",", line_end="\r\n"):
    lines = []
    with open(fn) as f:
        for line in f:
            lines.append(line.rstrip(line_end).split(field_delim))
    return lines


def find_block(csv, name):
    """For an Illumina SampleSheet.csv, return a tuple of the index of the
    line containing the header specified by name, and the index of the
    line just past the end of the data block. `range(*r)` will index
    all lines for the block, starting at the header line.
    """
    start = None
    end = None

    def blockend(f):
        maxfieldlen = max([0] + [len(x) for x in f])
        if len(f) > 0 and len(f[0]) > 0 and f[0][0] == "[" and f[0][-1] == "]":
            return True
        return 0 == maxfieldlen

    for i, fields in enumerate(csv):
        if len(fields) > 0 and fields[0] == name:
            start = i
        elif start is not None and blockend(fields):
            return start, i
    if start is not None:
        end = len(csv)
    return start, end


def rewrite_sample_sheet(sample_sheet_fn, umi_bp, index_bp):
    c = read_csv(sample_sheet_fn)
    width = max([len(x) for x in c])

    rstart, rend = find_block(c, "[Reads]")
    if rend - rstart != 3:
        raise ValueError("SampleSheet.csv does not have exactly two reads")
    iline = rend - 1
    newline = [str(umi_bp)] + [""] * (width - 1)
    c = c[0:iline] + [newline] + c[iline:]

    dstart, dend = find_block(c, "[Data]")
    if dstart is None or dend - dstart < 3:
        raise ValueError("No samples in [Data] section")
    dataheader = c[dstart + 1]
    if "index" not in dataheader or "index2" not in dataheader:
        raise ValueError("[Data] header does not contain index and index2")
    icol = dataheader.index("index")
    i2col = dataheader.index("index2")
    for drow in range(dstart + 2, dend):
        i = c[drow][icol]
        i2 = c[drow][i2col]
        if len(i[index_bp:].rstrip("N")) > 0 or len(i2[index_bp:].rstrip("N")) > 0:
            raise ValueError("Adapters are longer than {} bp".format(index_bp), i, i2)
        c[drow][icol] = i[:index_bp]
        c[drow][i2col] = i2[:index_bp]
        if len(c[drow][icol]) < index_bp or len(c[drow][i2col]) < index_bp:
            raise ValueError("Adapters are shorter than {}bp".format(index_bp), i, i2)

    sstart, send = find_block(c, "[Settings]")
    c = (
        c[0:send]
        + [
            ["TrimUMI", "1"] + [""] * (width - 2),
            ["Read2UMIStartFromCycle", "1"] + [""] * (width - 2),
            ["Read2UMILength", str(umi_bp)] + [""] * (width - 2),
        ]
        + c[send:]
    )

    return "\n".join([",".join(x) for x in c]) + "\n"


def write_run_dir(indir, sample_sheet_fn, run_info_fn, umi_bp, index_bp, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ri, ri_read_data = read_run_info(run_info_fn, umi_bp, index_bp)
    ri_out_fn = os.path.join(outdir, DEFAULT_RUN_INFO)
    if os.path.exists(ri_out_fn):
        raise ValueError("Output RunInfo.xml already exists, will not overwrite")
    ri.write(ri_out_fn)

    si = rewrite_sample_sheet(sample_sheet_fn, umi_bp, index_bp)
    si_out_fn = os.path.join(outdir, DEFAULT_SAMPLE_SHEET)
    if os.path.exists(si_out_fn):
        raise ValueError("Output SampleSheet.xml already exists, will not overwrite")
    with open(si_out_fn, mode="w") as sout:
        sout.write(si)

    bclcmd_fn = os.path.join(outdir, DEFALUT_BCL_CMD_FN)
    if os.path.exists(bclcmd_fn):
        raise ValueError("Output bcl2fastq.sh already exists, will not overwrite")
    with open(bclcmd_fn, mode="w") as cout:
        cout.write(
            BCL2FASTQ_TEMPLATE.format(
                BCL_DIR=os.path.abspath(indir), OUT_DIR=os.path.abspath(outdir)
            )
        )
        for sample, num in get_run_samples(si_out_fn):
            oldname = sample + "_S" + str(num) + "_R3_001.fastq.gz"
            newname = sample + "_S" + str(num) + "_R2_001.fastq.gz"
            oldpath = os.path.join(outdir, oldname)
            newpath = os.path.join(outdir, newname)
            cout.write("mv {} {}\n".format(oldpath, newpath))
    os.chmod(bclcmd_fn, os.stat(bclcmd_fn).st_mode | stat.S_IEXEC)
    return bclcmd_fn


def get_run_samples(ss_fn):
    c = read_csv(ss_fn)
    dstart, dend = find_block(c, "[Data]")
    name_idx = c[dstart + 1].index("Sample_Name")
    return [(c[i][name_idx], i - dstart - 1) for i in range(dstart + 2, dend)]


def srslyumi(indir, sample_sheet_fn, run_info_fn, umi_bp, index_bp, outdir):
    bclcmd_fn = write_run_dir(
        indir, sample_sheet_fn, run_info_fn, umi_bp, index_bp, outdir
    )
    subprocess.call(bclcmd_fn)


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument(
        "rundir", help=("Root directory of the run (typically containing RunInfo.xml)"),
    )
    ap.add_argument(
        "--input-dir",
        help=(
            "Directory containing BCLs "
            "(default: <rundir>/Data/Intensities/BaseCalls)"
        ),
    )
    ap.add_argument(
        "--run-info", help="path to RunInfo.xml (default: <rundir>/RunInfo.xml",
    )
    ap.add_argument(
        "--sample-sheet",
        help="path to SampleSheet.csv (default: <rundir>/SampleSheet.csv",
    )
    ap.add_argument(
        "--umi-bp",
        type=int,
        default=10,
        help="number of cycles in the index UMI (default:10)",
    )
    ap.add_argument(
        "--index-bp",
        type=int,
        default=7,
        help="number of cycles in the index barcode (default:7)",
    )
    ap.add_argument("outdir", help="output directory for FASTQ")

    a = ap.parse_args()

    if a.input_dir is None:
        a.input_dir = os.path.join(a.rundir, DEFAULT_BCL_SUBDIR)
    if a.sample_sheet is None:
        a.sample_sheet = os.path.join(a.rundir, DEFAULT_SAMPLE_SHEET)
    if a.run_info is None:
        a.run_info = os.path.join(a.rundir, DEFAULT_RUN_INFO)

    srslyumi(a.input_dir, a.sample_sheet, a.run_info, a.umi_bp, a.index_bp, a.outdir)
