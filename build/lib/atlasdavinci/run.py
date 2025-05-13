import sys,pysam
from atlasdavinci.utils import (
    check_bam_file, 
    check_bam_tags,
    check_and_create_bam_index,
    count_bam_reads,
    fix_bam_tags_DF_DC,
    batch_fold_structures,
    merge_fold_files
)
import random,string,logging,os
import numpy as np
from datetime import datetime
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description='Process BAM files with Davinci algorithm for Nanopore sequencing data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('bamfile', 
                      help='Input BAM file to process')
    parser.add_argument('-o', '--output',
                      help='Output BAM file name (default: input_bam_folded.bam)',
                      default=None)
    parser.add_argument('-t', '--temp-dir',
                      help='Temporary directory for intermediate files',
                      default=None)
    parser.add_argument('-v', '--verbose',
                      help='Increase output verbosity',
                      action='store_true')
    return parser.parse_args()

def main():
    logo = r'''
          _   _             ____  _       _        __
     /\  | | | |           |  _ \(_)     (_)      / _|
    /  \ | |_| | __ _ ___  | |_) |_  ___  _ _ __ | |_ ___
   / /\ \| __| |/ _` / __| |  _ <| |/ _ \| | '_ \|  _/ _ \
  / ____ \ |_| | (_| \__ \ | |_) | | (_) | | | | | || (_) |
 /_/    \_\__|_|\__,_|___/ |____/|_|\___/|_|_| |_|_| \___/

        `-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
        `=`,'=/     `=`,'=/     `=`,'=/     `=`,'=/
            y==/        y==/        y==/        y==/
        ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,-<=`.
        ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
    '''
    description_text = '''{}
     Process BAM files with Davinci algorithm for Nanopore sequencing data.
    '''.format(logo)
    parser = argparse.ArgumentParser(
        description=description_text,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('bamfile', 
                      help='Input BAM file to process')
    parser.add_argument('-o', '--output',
                      help='Output BAM file name (default: input_bam_folded.bam)',
                      default=None)
    parser.add_argument('-t', '--temp-dir',
                      help='Temporary directory for intermediate files',
                      default=None)
    parser.add_argument('-v', '--verbose',
                      help='Increase output verbosity',
                      action='store_true')
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Set up temporary directory
    if args.temp_dir is None:
        date_str = datetime.now().strftime("%Y%m%d")
        random_chars = ''.join(random.choices(string.ascii_letters + string.digits, k=4))
        temp_dir = f"{date_str}{random_chars}_TEMP"
    else:
        temp_dir = args.temp_dir
    
    bamfile = args.bamfile
    
    if not check_bam_file(bamfile):
        logging.error(f"Error: {bamfile} is not a valid BAM file or does not have an index")
        sys.exit(1)
    
    total_reads, mapped_reads = count_bam_reads(bamfile)
    logging.info(f"\nBAM file statistics for {bamfile}:")
    logging.info(f"Total reads: {total_reads:,}")
    logging.info(f"Mapped reads: {mapped_reads:,}")
    logging.info(f"Mapping rate: {(mapped_reads/total_reads*100):.2f}%")

    if not check_bam_tags(bamfile):
        logging.info("\nTags are missing, running fix_bam_tags_DF_DC...")
        random_str = ''.join(random.choices(string.ascii_letters, k=6))
        output_bam = bamfile.replace('.bam', f'{random_str}.bam')
        fix_bam_tags_DF_DC(bamfile, output_bam)
        os.replace(output_bam, bamfile)
        check_and_create_bam_index(bamfile)
    
    # Run batch folding
    temp_dir = batch_fold_structures(bamfile)
    if temp_dir is None:
        logging.error("Folding process failed")
        sys.exit(1)
    logging.info(f"Folding results are in: {temp_dir}")
    
    # Merge results into new BAM file
    output_bam = merge_fold_files(temp_dir, bamfile)
    if output_bam is None:
        logging.error("Failed to merge folding results")
        sys.exit(1)
    
    # Rename output file if specified
    if args.output:
        os.rename(output_bam, args.output)
        output_bam = args.output
        
    logging.info(f"Successfully created BAM file with folding information: {output_bam}")

if __name__ == "__main__":
    main() 