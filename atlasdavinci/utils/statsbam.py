def count_bam_reads(bam_file):
    """
    Count total reads and mapped reads in a BAM file
    
    Args:
        bam_file (str): Path to BAM file
        
    Returns:
        tuple: (total_reads, mapped_reads)
    """
    import pysam
    
    total_reads = 0
    mapped_reads = 0
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                total_reads += 1
                if not read.is_unmapped:
                    mapped_reads += 1
                    
        return total_reads, mapped_reads
        
    except Exception as e:
        print(f"Error processing BAM file: {e}")
        return 0, 0
