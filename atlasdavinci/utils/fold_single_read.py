#!/usr/bin/env python3
import os
import sys
import subprocess
import logging

def fold_single_read(read_id, temp_dir):
    """
    Process folding for a single read
    
    Args:
        read_id (str): Read ID
        temp_dir (str): Base directory for temporary files
    """
    constr_dir = os.path.join(temp_dir, "constr")
    folded_dir = os.path.join(temp_dir, "folded")
    # posteriors_dir = os.path.join(temp_dir, "posteriors_output")
    
    constr_file = os.path.join(constr_dir, f"{read_id}.bpseq")
    fold_file = os.path.join(folded_dir, f"{read_id}.fold")
    # posteriors_file = os.path.join(posteriors_dir, f"{read_id}.txt")
    
    # Run contrafold
    try:
        pipe = ["contrafold", "predict", "--constraints", constr_file, 
                "--parens", fold_file]
        
        result = subprocess.run(pipe, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"Error running contrafold for {read_id}: {result.stderr}")
            return False
            
        # # Convert to dotbracket format
        # db_file = os.path.join(folded_dir, f"{read_id}.db")
        # pipe = ["fold2dotbracketFasta.py", "--input_file", fold_file,
        #         "--tag", read_id,
        #         "--output_file", db_file]
        
        # result = subprocess.run(pipe, capture_output=True, text=True)
        # if result.returncode != 0:
        #     logging.error(f"Error converting to dotbracket for {read_id}: {result.stderr}")
        #     return False
            
        # # Convert to element string
        # element_file = os.path.join(folded_dir, f"{read_id}.txt")
        # pipe = ["rnaConvert.py", db_file, '-T', 'element_string', 
        #         '--force', '--to-file', '--filename', element_file]
        
        # result = subprocess.run(pipe, capture_output=True, text=True)
        # if result.returncode != 0:
        #     logging.error(f"Error converting to element string for {read_id}: {result.stderr}")
        #     return False
            
        return True
        
    except Exception as e:
        logging.error(f"Error processing {read_id}: {str(e)}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: fold_single_read.py <read_id> <temp_dir>")
        sys.exit(1)
        
    read_id = sys.argv[1]
    temp_dir = sys.argv[2]
    
    logging.basicConfig(level=logging.INFO)
    success = fold_single_read(read_id, temp_dir)
    sys.exit(0 if success else 1) 