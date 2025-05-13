from setuptools import setup, find_packages

setup(
    name="atlasdavinci",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.0.0",
        "pysam>=0.16.0"
    ],
    entry_points={
        'console_scripts': [
            'atlasdavinci=atlasdavinci.run:main',
        ],
    },
    author="Atlas Bioinformatics",
    author_email="atlas.bioinfo@example.com",
    description="A tool for processing Nanopore sequencing BAM files with Davinci algorithm",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/atlasbioinfo/atlasdavinci",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.7",
) 