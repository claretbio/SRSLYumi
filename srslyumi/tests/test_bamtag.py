import filecmp
import unittest
import tempfile

from srslyumi.tests.test_cli import capture_cli
from srslyumi.tests.test_cli import f

from srslyumi.bamtag import main


class TestBamTag(unittest.TestCase):
    def test_help(self):
        argv = ["bamtag"]
        with capture_cli(argv) as (stdout, stderr):
            with self.assertRaises(SystemExit) as e:
                main()

    def test_expected(self):
        with tempfile.NamedTemporaryFile(mode="w+") as out:
            argv = ["bamtag", "-o", out.name, f("bamtag-in-01.sam")]
            with capture_cli(argv) as (stdout, stderr):
                main()
            self.assertEqual("", stderr.getvalue())
            self.assertEqual("", stdout.getvalue())
            out.seek(0)
            with open(f("bamtag-out-01.sam")) as expected:
                self.assertListEqual(out.readlines(), expected.readlines())

    def test_noumi(self):
        with tempfile.NamedTemporaryFile(mode="w+") as out:
            argv = ["bamtag", "-o", out.name, f("bamtag-in-02.sam")]
            with capture_cli(argv) as (stdout, stderr):
                main()
            self.assertEqual("", stdout.getvalue())
            self.assertEqual(
                "WARNING: 3 of 3 reads did not have a UMI_like string\n",
                stderr.getvalue(),
            )
            out.seek(0)
            with open(f("bamtag-out-02.sam")) as expected:
                self.assertListEqual(out.readlines(), expected.readlines())

    def test_take_fragment(self):
        with tempfile.NamedTemporaryFile(mode="w+") as out:
            argv = [
                "bamtag",
                "-o",
                out.name,
                "--take-fragment",
                "1",
                f("bamtag-in-03.sam"),
            ]
            with capture_cli(argv) as (stdout, stderr):
                main()
            self.assertEqual("", stdout.getvalue())
            self.assertEqual("", stderr.getvalue())
            out.seek(0)
            with open(f("bamtag-out-03.sam")) as expected:
                self.assertListEqual(out.readlines(), expected.readlines())

    def test_bam(self):
        with tempfile.NamedTemporaryFile(mode="w+") as out:
            argv = ["bamtag", "-b", "-o", out.name, f("bamtag-in-01.sam")]
            with capture_cli(argv) as (stdout, stderr):
                main()
            self.assertEqual("", stdout.getvalue())
            self.assertEqual("", stderr.getvalue())
            self.assertTrue(
                filecmp.cmp(out.name, f("bamtag-out-01.bam")), "BAM output differs"
            )
