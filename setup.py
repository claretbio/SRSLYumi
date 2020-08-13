import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="srslyumi",
    version="0.3",
    author="Charles Vaske",
    author_email="charles.vaske@claretbio.com",
    description="process SRSLY UMIs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/claretbio/SRSLYumi",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "srslyumi=srslyumi.cli:main",
            "srslyumi-bamtag=srslyumi.bamtag:main",
        ]
    },
    install_requires=['pysam>=0.15.3'],
    python_requires=">=2.7",
)
