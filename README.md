![build](https://github.com/claretbio/SRSLYumi/workflows/build/badge.svg)

# SRSLY UMI processing

SRSLY UMIs are attached to the i7 index, and require a bit of handling
to make it through `bcl2fastq`. This package helps guide that process.

# SRSLY UMI dual-index sequencing runs

Illumina sequencing performs read cycles for the i7 and i5 indices in
between the fragment reads. SRSLY UMIs are attached to the i7
reads, like this:

```ditaa
Order of read cycles
+------------------------+--------------+------------+----------------+-----------------------+
| Forward fragment read  | i7 index     |   UMI      | i5 index       | Reverse fragment read |
+------------------------+--------------+------------+----------------+-----------------------+

Standard bcl2fastq processing
+------------------------+---------------------------+----------------+-----------------------+
| RunInfo.xml Read 1     | RunInfo.xml R2            | RunInfo.xml R3 | RunInfo.xml R4        |
+------------------------+---------------------------+----------------+-----------------------+
| FASTQ output R1        | FASTQ header, without UMI                  | FASTQ output R2       |
+------------------------+--------------------------------------------+-----------------------+

Reformatted bcl2fastq processing
+------------------------+--------------+------------+----------------+-----------------------+
| new RunInfo.xml Read 1 | R2 IsIndex   | R3         | R4 IsIndex     | R5                    |
+------------------------+--------------+------------+----------------+-----------------------+
| fixed FASTQ output R1  | FASTQ header, with UMI                     | fixed FASTQ output R3 |
+------------------------+--------------------------------------------+-----------------------+
```

However, bcl2fastq can't insert UMIs into the fragment name in the
FASTQ header unless it is part of etiher the output R1 or output
R2. (Note that index reads as defined in the RunInfo.xml do not count
as output reads.

So to solve this, we define a new `RunInfo.xml` that defines five
reads instead of four:

With standard `bcl2fastq` processing with the `TrimUMI` option, this
results in the UMI in the fragment name in the FASTQ files. However,
it has two side effects: the UMI includes both the 5bp of the UMI as
well as followed by the first 5bp of the actual read2. This should be
compatible with most UMI analysis. Additionally, the proper R2 file is
labeled as R3. Post-processing can easily rename the R3 to R2.

# Using this package

After installation of this python package, the `srslyumi` command will
take two arguments: 1) an existing run directory, and 2) an output
directory for the FASTQ and `bcl2fastq` reports. Inside this output
directory, a new `RunInfo.xml` and `SampleSheet.csv` will be created,
along with a `bcl2fastq2.sh` command that can be used to rerun the
process. Note that at the end of this command, the `_R3_` files are
renamed to `_R2_`.

# Developing this package further

When your working directory is the root of this repository, the same
directory that contains `setup.py`, you can run `pip install -e .` to
install the package in a form that lets you edit your code and run it
as a python package at the same time.

## Testing and test coverage

During delevopment `tox`, will setup testing virtual environments for
Python 2.7 and Python 3.6 and run all tests. Before code can accepted
to the main repository it must pass test on Python 2.7, 3.5, 3.6, 3.7,
and 3.8, which will run on GitHub automatically.

For quick tests in your current Python environment, run `pytest`,
though you may need to install the test dependencies as listed under
the tox section of `pyproject.toml`.

To run quick tests in your current environment, run `pytest`

To assess code coverage of the tests, run `pytest --cov --report=html`.
