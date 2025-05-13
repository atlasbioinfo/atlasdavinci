import os
import pysam

def check_file_exists(file_path):
    """
    Check if a file exists at the given path
    
    Args:
        file_path (str): Path to the file to check
        
    Returns:
        bool: True if file exists, False otherwise
    """
    return os.path.exists(file_path)

def check_bam_file(file_path):
    """
    Check if a file is a valid BAM file and has an index
    
    Args:
        file_path (str): Path to the BAM file to check
        
    Returns:
        bool: True if valid BAM file with index, False otherwise
    """
    if not check_file_exists(file_path):
        return False
        
    try:
        # Try to open as BAM file
        bam = pysam.AlignmentFile(file_path, "rb")
        # Check if index is older than BAM file
        bai_path = file_path + '.bai'
        alt_bai_path = file_path.replace('.bam', '.bai')
        
        if os.path.exists(bai_path):
            if os.path.getmtime(bai_path) < os.path.getmtime(file_path):
                os.remove(bai_path)
        elif os.path.exists(alt_bai_path):
            if os.path.getmtime(alt_bai_path) < os.path.getmtime(file_path):
                os.remove(alt_bai_path)
        check_and_create_bam_index(file_path)
        bam.close()
        return True
        
    except (ValueError, OSError):
        return False

def check_bam_tags(bamfile):
    """
    Check if the BAM file has DF and DC tags in mapped reads
    
    Args:
        bamfile (str): Path to the BAM file to check
        
    Returns:
        bool: True if all mapped reads have both DF and DC tags, False otherwise
    """
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            has_rf = read.has_tag('DF')
            has_ac = read.has_tag('DC')
            if not (has_rf and has_ac):
                return False
            break
    return True

def check_and_create_bam_index(bamfile):
    """
    Check if BAM index exists, create one if not
    
    Args:
        bamfile (str): Path to the BAM file to check/index
        
    Returns:
        bool: True if index exists or was created successfully, False if indexing failed
    """
    # Check if index exists
    if os.path.exists(bamfile + '.bai') or os.path.exists(bamfile.replace('.bam', '.bai')):
        return True
        
    try:
        # Create index if it doesn't exist
        pysam.index(bamfile)
        return True
    except Exception:
        return False
