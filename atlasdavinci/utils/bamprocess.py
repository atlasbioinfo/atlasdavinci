import numpy as np
import pysam
import re
import os
import subprocess
import logging
from datetime import datetime
import random
import string

def remove_soft_clipping(read):
    """
    Remove soft clipping from read sequence and aligned pairs
    
    Args:
        read: pysam.AlignedSegment object
        
    Returns:
        tuple: (trimmed_query_seq, trimmed_pairs)
    """
    tquery_seq = read.query_sequence
    pairs = np.array(read.get_aligned_pairs())
    
    if read.cigartuples[0][0] == 4:  # Soft clip at start
        pairs = pairs[read.cigartuples[0][1]:]
        tquery_seq = tquery_seq[read.cigartuples[0][1]:]
    if read.cigartuples[-1][0] == 4:  # Soft clip at end
        pairs = pairs[:-read.cigartuples[-1][1]]
        tquery_seq = tquery_seq[:-read.cigartuples[-1][1]]
        
    return tquery_seq, pairs

def fix_bam_tags_DF_DC(bamfile, output_bam):
    with pysam.AlignmentFile(bamfile, "rb") as in_bam:
        header = in_bam.header.copy()
        with pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
            for read in in_bam:
                if read.is_unmapped:
                    continue
                
                # Get reference sequence
                ref_seq = read.get_reference_sequence()
                read.set_tag('DF', ref_seq)
                
                # Process read sequence and create constraint string
                tquery_seq, pairs = remove_soft_clipping(read)
                
                
                constraint = ['0'] * len(ref_seq)
                
                # Process MD tag for mutations and deletions
                mdtag = read.get_tag("MD")
                tread = tquery_seq
                tref = ref_seq
                tpos = [p - read.reference_start for p in read.get_reference_positions()]
                
                for match in re.finditer(r'(\d+)|(\^[A-Z]+)|([A-Z])', mdtag):
                    if match.group(1):  # Match
                        tpos = tpos[int(match.group(1)):]
                        tread = tread[int(match.group(1)):]
                        tref = tref[int(match.group(1)):]
                    elif match.group(2):  # Deletion
                        deletion = match.group(2)[1:]
                        pos = tpos[0] - 1
                        constraint[pos] = '1'  # Mark deletion position
                        tref = tref[len(deletion):]
                    elif match.group(3):  # Mutation
                        pos = tpos[0]
                        constraint[pos] = '1'  # Mark mutation position
                        tpos = tpos[1:]
                        tread = tread[1:]
                        tref = tref[1:]
                
                
                if any(op[0] == 1 for op in read.cigartuples):  # If there are insertions
                    mask = (pairs[:, 0] != None) & (pairs[:, 1] == None)
                    mask_indices = np.where(mask)[0]
                    
                    if len(mask_indices) > 0:
                        splits = np.where(np.diff(mask_indices) > 1)[0] + 1
                        insertion_groups = np.split(mask_indices, splits)
                        
                        for group in insertion_groups:
                            if len(group) > 0:
                                # Mark position before insertion
                                pos = pairs[min(group)-1, 1] - read.reference_start
                                if pos >= 0 and pos < len(constraint):
                                    constraint[pos] = '1'
                
                constraint="".join(constraint)
                read.set_tag('DC', constraint)
                out_bam.write(read)
    
    
    # Index the new BAM file
    pysam.index(output_bam)

def batch_fold_structures(bamfile):
    """
    Process folding for all reads in a BAM file
    
    Args:
        bamfile (str): Path to the BAM file
        
    Returns:
        str: Path to the temporary directory containing folding results
    """
    # Create temporary directory with date and random string
    date_str = datetime.now().strftime("%Y%m%d")
    random_chars = ''.join(random.choices(string.ascii_letters + string.digits, k=4))
    temp_dir = f"{date_str}{random_chars}_TEMP"
    
    # Create subdirectories
    constr_dir = os.path.join(temp_dir, "constr")
    folded_dir = os.path.join(temp_dir, "folded")
    # posteriors_dir = os.path.join(temp_dir, "posteriors_output")
    
    os.makedirs(constr_dir, exist_ok=True)
    os.makedirs(folded_dir, exist_ok=True)
    # os.makedirs(posteriors_dir, exist_ok=True)
    
    # List to store read IDs
    id_list = []
    
    # Extract reads and create constraint files
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            rf_tag = read.get_tag('DF')
            ac_tag = read.get_tag('DC')
            if ac_tag.find("0") == -1:
                continue
            read_id = read.query_name
            id_list.append(read_id)
            
            constr_file = os.path.join(constr_dir, f"{read_id}.bpseq")
            with open(constr_file, 'w') as f:
                for i, (base, bit) in enumerate(zip(rf_tag, ac_tag), 1):
                    constraint = '0' if bit == '1' else '-1'
                    f.write(f"{i}\t{base}\t{constraint}\n")
    
    # Write read IDs to a file for xargs
    id_file = os.path.join(temp_dir, "read_ids.txt")
    with open(id_file, 'w') as f:
        for read_id in id_list:
            f.write(f"{read_id}\n")
    
    # Run folding in parallel using xargs
    num_cores = os.cpu_count() or 4  # Use all available cores or default to 4
    fold_script = os.path.join(os.path.dirname(__file__), "fold_single_read.py")
    
    cmd = f"cat {id_file} | xargs -P {num_cores} -I {{}} python {fold_script} {{}} {temp_dir}"
    logging.info(f"Running folding in parallel using {num_cores} cores")
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        logging.info("Folding completed successfully")
        return temp_dir
    except subprocess.CalledProcessError as e:
        logging.error(f"Error during folding: {e}")
        return None
    
def merge_fold_files(temp_dir, bamfile):
    """
    Merge folding results and add them to BAM file
    
    Args:
        temp_dir (str): Directory containing folding results
        bamfile (str): Path to input BAM file
        
    Returns:
        str: Path to the new BAM file with folding information
    """
    constr_dir = os.path.join(temp_dir, "constr")
    folded_dir = os.path.join(temp_dir, "folded")
    # posteriors_dir = os.path.join(temp_dir, "posteriors_output")
    
    # Get list of files in each directory
    constr_files = set(f.split('.')[0] for f in os.listdir(constr_dir) if f.endswith('.bpseq'))
    fold_files = set(f.split('.')[0] for f in os.listdir(folded_dir) if f.endswith('.fold'))
    # post_files = set(f.split('.')[0] for f in os.listdir(posteriors_dir) if f.endswith('.txt'))
    
    # Check if all directories have the same number of files
    if not (constr_files == fold_files):
        logging.error("Mismatch in number of files between directories")
        return None
    
    tnames = bamfile.split(".")
    tnames[0] = tnames[0] + '_folded'
    output_bam = ".".join(tnames)
    
    with pysam.AlignmentFile(bamfile, "rb") as in_bam:
        header = in_bam.header.copy()
        
        with pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
            for read in in_bam:
                if read.is_unmapped:
                    out_bam.write(read)
                    continue
                    
                read_id = read.query_name
                if read_id not in constr_files:
                    out_bam.write(read)
                    continue
                
                fold_file = os.path.join(folded_dir, f"{read_id}.fold")
                with open(fold_file) as f:
                    f.readline()
                    f.readline()
                    f.readline()
                    structure = f.readline().strip()
                read.set_tag('DS', structure)
                out_bam.write(read)
    

    pysam.index(output_bam)
    logging.info(f"Created new BAM file with folding information: {output_bam}")
    return output_bam
