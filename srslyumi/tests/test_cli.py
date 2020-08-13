from contextlib import contextmanager
import filecmp

# Python 2.7 io.SringIO can only take unicode writes, but many packges
# will write with str instead, which is outside our control. So we try
# the 2.7 import if we can, otherwise we go to the Python 3 package.
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import os.path
import tempfile
import unittest

from mock import patch

import xml.etree.ElementTree as ET
import sys

from xmldiff.main import diff_files as xml_diff_files

from srslyumi.cli import main
from srslyumi.cli import find_block
from srslyumi.cli import read_csv
from srslyumi.cli import read_run_info
from srslyumi.cli import rewrite_sample_sheet
from srslyumi.cli import get_run_samples


THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def f(fn):
    return os.path.join(THIS_DIR, fn)


@contextmanager
def capture_cli(argv):
    """This context manager captures the standard output and error of a
    function, as well as setting the command line arguments of
    sys.argv. This is meant for testing functions that are used as a
    command line interface.
    """
    new_out, new_err = StringIO(""), StringIO("")
    with patch.object(sys, "argv", argv):
        with patch.object(sys, "stdout", new_out):
            with patch.object(sys, "stderr", new_err):
                yield sys.stdout, sys.stderr


class Main(unittest.TestCase):
    def test_usage(self):
        argv = ["foobar-run"]
        with capture_cli(argv) as (stdout, stderr):
            with self.assertRaises(SystemExit) as e:
                main()
        assert e.exception.code == 2
        self.assertEqual("", stdout.getvalue())
        self.assertNotEqual("", stderr.getvalue())

    def test_integration(self):
        test_dir = tempfile.mkdtemp()
        argv = ["srslyumi", "--index-bp", "8", "--umi-bp", "5", f("TestInput"), test_dir]
        origpath = os.environ["PATH"]
        os.environ["PATH"] = f("") + os.pathsep + origpath
        with capture_cli(argv) as (stdout, stderr):
            # follow code path of srslyumi that creates the output directory
            os.rmdir(test_dir)
            main()
        os.environ["PATH"] = origpath
        dircmp = filecmp.dircmp(f("TestOutput"), test_dir)
        # the paths in this command will differ since the directory name is random
        self.assertTrue("runbcl2fastq.sh" in dircmp.diff_files)
        # Depending on the Python version, the XML can be serialized in different ways
        self.assertTrue("RunInfo.xml" in dircmp.diff_files + dircmp.same_files)
        self.assertTrue("SampleSheet.csv" in dircmp.same_files)
        self.assertEqual(0, len(dircmp.funny_files))

        with self.assertRaises(ValueError) as e:
            with capture_cli(argv) as (stdout, stderr):
                main()
        os.remove(os.path.join(test_dir, "RunInfo.xml"))
        with self.assertRaises(ValueError) as e:
            with capture_cli(argv) as (stdout, stderr):
                main()
        os.remove(os.path.join(test_dir, "RunInfo.xml"))
        os.remove(os.path.join(test_dir, "SampleSheet.csv"))
        with self.assertRaises(ValueError) as e:
            with capture_cli(argv) as (stdout, stderr):
                main()


class RunInfo(unittest.TestCase):
    def test_missing_file(self):
        with self.assertRaises(EnvironmentError) as e:
            ri = read_run_info("missing", 5, 8)

    def test_01_good(self):
        ri, read_data = read_run_info(f("RunInfo-input01-good.xml"), 5, 8)
        self.assertNotEqual(None, ri)
        self.assertEqual(4, len(read_data))
        rifile = StringIO(ET.tostring(ri.getroot()).decode("utf-8"))

        diff = xml_diff_files(f("RunInfo-output01-good.xml"), rifile)
        self.assertEqual([], diff)

    def test_02_no_reads(self):
        with self.assertRaises(ValueError) as e:
            read_run_info(f("RunInfo-input02-bad-no-reads.xml"), 5, 8)

    def test_03_single_index(self):
        with self.assertRaises(ValueError) as e:
            read_run_info(f("RunInfo-input03-bad-single-index.xml"), 5, 8)

    def test_04_no_umi_cycles(self):
        with self.assertRaises(ValueError) as e:
            read_run_info(f("RunInfo-input04-bad-no-umi-cycles.xml"), 5, 8)

    def test_05_r2notindex(self):
        with self.assertRaises(ValueError) as e:
            read_run_info(f("RunInfo-input05-bad-r2notindex.xml"), 5, 8)

    ## Removing this test, as it's overly restrictive
    # def test_06_shorti5(self):
    #     with self.assertRaises(ValueError) as e:
    #         read_run_info(f("RunInfo-input06-bad-shorti5.xml"), 5, 8)

    def test_07_r3notindex(self):
        with self.assertRaises(ValueError) as e:
            read_run_info(f("RunInfo-input07-bad-r3notindex.xml"), 5, 8)


class SampleSheet(unittest.TestCase):
    def test_missing_file(self):
        with self.assertRaises(EnvironmentError):
            read_csv("missing")

    def test_01_good_sample_sheet(self):
        c = read_csv(f("SampleSheet-input01-good.csv"))
        self.assertEqual(33, len(c))
        maxfields = max([len(x) for x in c])
        self.assertEqual(10, maxfields)

    def test_missing_data(self):
        bad_sample_sheets = [
            "SampleSheet-input03-bad-no-samples.csv",
            "SampleSheet-input04-bad-nodataheader.csv",
            "SampleSheet-input05-bad-nodatasection.csv",
            "SampleSheet-input06-noindex2.csv",
            "SampleSheet-input07-longadapters.csv",
            "SampleSheet-input08-shortadapters.csv",
            "SampleSheet-input09-bad-fourreads.csv",
        ]
        for ss_fn in bad_sample_sheets:
            with self.assertRaises(ValueError) as e:
                rewrite_sample_sheet(f(ss_fn), 5, 8)

    def test_data_detection_01(self):
        c = read_csv(f("SampleSheet-input01-good.csv"))
        r = find_block(c, "[Data]")
        self.assertEqual((None, None), find_block(c, "[missing]"))
        self.assertEqual((0, 10), find_block(c, "[Header]"))
        self.assertEqual((11, 14), find_block(c, "[Reads]"))
        self.assertEqual((15, 17), find_block(c, "[Settings]"))
        self.assertEqual((18, 33), find_block(c, "[Data]"))

    def test_data_detection_02(self):
        c = read_csv(f("SampleSheet-input02-good-no-blanks.csv"))
        r = find_block(c, "[Data]")
        self.assertEqual((None, None), find_block(c, "[missing]"))
        self.assertEqual((0, 10), find_block(c, "[Header]"))
        self.assertEqual((10, 13), find_block(c, "[Reads]"))
        self.assertEqual((13, 15), find_block(c, "[Settings]"))
        self.assertEqual((15, 30), find_block(c, "[Data]"))

    def test_rewrite(self):
        fn = f("SampleSheet-input01-good.csv")
        out = rewrite_sample_sheet(fn, 5, 8)
        with open(f("SampleSheet-output01-good.csv")) as efile:
            expected = efile.read()
        self.assertEqual(expected.split("\n"), out.split("\n"))

    def test_get_run_samples(self):
        samples = [
            ("SR2167_hu", 1),
            ("SR2168_hu", 2),
            ("SR2169_hu", 3),
            ("SR2170_hu", 4),
            ("SR2171_hu", 5),
            ("SR2172_hu", 6),
            ("SR2173_hu", 7),
            ("SR2174_hu", 8),
            ("SR2175_hu", 9),
            ("SR2176_hu", 10),
            ("SR2177_hu", 11),
            ("SR2178_hu", 12),
            ("SR2179_hu", 13),
        ]
        o = get_run_samples(f("SampleSheet-input01-good.csv"))
        self.assertEqual(samples, o)
        o = get_run_samples(f("SampleSheet-output01-good.csv"))
        self.assertEqual(samples, o)
