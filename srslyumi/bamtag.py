"""
Move bcl2fastq read name UMIs to a BAM tag
"""
import argparse
import re
import sys

import pysam


UMI_REGULAR_EXPRESSION = "^[acgtnAGCTN_+-]+$"


def looks_like_umi(umi):
    """
    >>> looks_like_umi("ACGT")
    True
    >>> looks_like_umi("AGCT+TCGA")
    True
    >>> looks_like_umi("AGCT_TCGA")
    True
    >>> looks_like_umi("AGCT:TCGA")
    False
    >>> looks_like_umi("1029")
    False
    >>> looks_like_umi("M02607")
    False
    >>> looks_like_umi("1GATATCCGTCGGGCC")
    False
    """
    return re.match(UMI_REGULAR_EXPRESSION, umi) is not None


def picard_friendly(umi):
    """
    Converts a UMI string into a format that Picard Tools can understand

    >>> picard_friendly("ACGT+ACGT")
    'ACGT-ACGT'
    >>> picard_friendly("ACGT_ACGT")
    'ACGT-ACGT'
    >>> picard_friendly("TTTAA")
    'TTTAA'
    """
    return re.sub("[^actgnACTGN]", "-", umi)


def take_fragment(umi, fragment_index):
    """
    >>> take_fragment("AGCT-TCGA", 0)
    'AGCT'
    >>> take_fragment("AGCT-TCGA", 1)
    'TCGA'
    >>> take_fragment("AGCT_TCGA", 0)
    'AGCT_TCGA'
    >>> take_fragment("AGCT", 1)
    Traceback (most recent call last):
      ...
    srslyumi.bamtag.FragmentIndexOutOfBounds: ('AGCT', 1)
    """
    frags = umi.split("-")
    if fragment_index < 0 or fragment_index >= len(frags):
        raise FragmentIndexOutOfBounds(umi, fragment_index)
    return frags[fragment_index]


def bamtag(inbam, out, sam_tag, keep_symbols, fragment_index, quiet):
    num_missing_umis = 0
    reads = 0
    for read in inbam.fetch(until_eof=True):
        reads += 1
        name, _, umi = read.query_name.rpartition(":")
        if looks_like_umi(umi):
            if not keep_symbols:
                umi = picard_friendly(umi)
            if fragment_index is not None:
                umi = take_fragment(umi, fragment_index)
            read.query_name = name
            read.set_tag(sam_tag, umi, "Z")
        else:
            num_missing_umis += 1
        out.write(read)
    if not quiet and num_missing_umis > 0:
        msg = "WARNING: {} of {} reads did not have a UMI_like string\n"
        sys.stderr.write(msg.format(num_missing_umis, reads))


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument(
        "--sam-tag", default="RX", help="SAM tag to use for the UMI infomation"
    )
    ap.add_argument(
        "--keep-symbols",
        action="store_true",
        help="do not replace special characters in UMI with '_' (needed for Picard)",
    )
    ap.add_argument(
        "--take-fragment",
        type=int,
        help="split UMI on '-' character, use only this fragment (0-based)",
    )
    ap.add_argument(
        "-b", "--binary", action="store_true", help="write BAM instead of sam"
    )
    ap.add_argument("-q", "--quiet", action="store_true", help="don't report warnings")
    ap.add_argument("-o", default=sys.stdout, help="output SAM/BAM (default: STDOUT)")
    ap.add_argument(
        "inputbam", default=sys.stdin, help="input SAM/BAM (default: STDIN)"
    )
    a = ap.parse_args()

    if (isinstance(a.o, str) and a.o.endswith("bam")) or a.binary:
        mode = "wb"
    else:
        mode = "w"

    inbam = pysam.AlignmentFile(a.inputbam)
    out = pysam.AlignmentFile(a.o, mode=mode, template=inbam)

    bamtag(inbam, out, a.sam_tag, a.keep_symbols, a.take_fragment, a.quiet)

    out.close()
    inbam.close()
