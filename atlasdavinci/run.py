import sys,pysam
from .utils.checkbam import check_bam_file, check_bam_tags,check_and_create_bam_index
from .utils.statsbam import count_bam_reads
from .utils.bamprocess import fix_bam_tags_DF_DC, batch_fold_structures, merge_fold_files
import random,string,logging,os
import numpy as np
from datetime import datetime
import subprocess

def main():
    logging.basicConfig(level=logging.INFO)
    date_str = datetime.now().strftime("%Y%m%d")
    random_chars = ''.join(random.choices(string.ascii_letters + string.digits, k=4))
    temp_dir = f"{date_str}{random_chars}_TEMP"
    bamfile = sys.argv[1]
    
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
    logging.info(f"Successfully created BAM file with folding information: {output_bam}")

if __name__ == "__main__":
    main() 