# atlas Davinci

This script is a simplified version of Davinci algorithm from the paper published in Nature (DOI: 10.1038/s41586-022-05135-9). It supports direct mutation calling using Davinci algorithm from Nanopore sequencing BAM files after alignment.

## Features

- Direct processing of Nanopore BAM files
- Implementation of core Davinci algorithm
- Simplified workflow compared to the original paper
- Variant calling optimized for Nanopore data
- Automatic BAM file validation and preprocessing
- Batch processing of read structures
- Comprehensive logging and progress tracking

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/atlasdavinci.git
cd atlasdavinci
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Install the package:
```bash
pip install .
```

## Dependencies

- Python 3.6 or higher
- numpy >= 1.19.0
- pandas >= 1.0.0
- pysam >= 0.16.0

## Usage

Basic usage:
```bash
python -m atlasdavinci.run input.bam -o output.bam
```

### Command Line Arguments

- `bamfile`: Input BAM file to process (required)
- `-o, --output`: Output BAM file name (default: input_bam_folded.bam)
- `-t, --temp-dir`: Temporary directory for intermediate files (default: auto-generated)
- `-v, --verbose`: Increase output verbosity

### Example

```bash
# Basic usage with default settings
python -m atlasdavinci.run sample.bam

# Specify output file and enable verbose logging
python -m atlasdavinci.run sample.bam -o processed_sample.bam -v

# Use custom temporary directory
python -m atlasdavinci.run sample.bam -t /path/to/temp/dir
```

## Output

The script generates a new BAM file containing the processed reads with additional structural information. The output file includes:
- Original alignment information
- Processed structural data
- Quality metrics and statistics

## Logging

The script provides detailed logging information including:
- BAM file statistics (total reads, mapped reads, mapping rate)
- Processing progress
- Error messages and warnings
- Final output file location

## Citation

If you use this software in your research, please cite the original Davinci paper:
```
[Paper citation details to be added]
```

## License

This project is licensed under the terms of the included LICENSE file.

