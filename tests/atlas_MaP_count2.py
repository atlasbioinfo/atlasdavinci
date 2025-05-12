import argparse,os
import logging,re,tqdm,logging
import pysam, gzip
import numpy as np
import subprocess
transNuc={"A":"T","T":"A","C":"G","G":"C","N":"N","-":"-"}

def run_RNA_MaP(input_file, output, del_thred, insertion_thred, all, gz):
    logging.info("run RNA MaP")
    bamfile = pysam.AlignmentFile(input_file, "rb")
    
    
    ref_stats = {}
    cov_file = input_file + '.cov'
    logging.info(f"Reading coverage file: {cov_file}")
    
    with open(cov_file) as f:
        for line in f:
            ref, pos, cov = line.strip().split('\t')
            pos = int(pos)
            if ref not in ref_stats:
                ref_stats[ref] = {
                    'start': pos,
                    'end': pos
                }
            else:
                ref_stats[ref]['end'] = pos

    with gzip.open(output + ".gz", "wt") if gz else open(output, "wt") as out:
        out.write(",".join(["ref", "pos", "coverageF", "mutF", "delF", "insF", 
                           "coverageR", "mutR", "delR", "insR"])+"\n")
        
        for tgeno in bamfile.header["SQ"]:
            ref_name = tgeno["SN"]
            if ref_name not in ref_stats:
                continue
            start_pos = ref_stats[ref_name]["start"]
            end_pos = ref_stats[ref_name]["end"]
            array_size = end_pos - start_pos
            mutCountF = np.zeros(array_size, dtype=int)
            mutCountR = np.zeros(array_size, dtype=int)
            coverageF = np.zeros(array_size, dtype=int)
            coverageR = np.zeros(array_size, dtype=int)
            mutInfoF = [[] for _ in range(array_size)]
            mutInfoR = [[] for _ in range(array_size)]
            delInfoF = [[] for _ in range(array_size)]
            delInfoR = [[] for _ in range(array_size)]
            insInfoF = [[] for _ in range(array_size)]
            insInfoR = [[] for _ in range(array_size)]
            
            for read in bamfile.fetch(tgeno["SN"]):
                if read.is_unmapped:
                    continue
                # if read.query_name != "V350236254L3C002R04500982377:0:0:0:0":
                #     continue
                # if read.flag != 147:
                #     continue
                # print(read)
                ForR = -1 if read.is_read1 == read.is_reverse else 1    
                tpos = [p - start_pos for p in read.get_reference_positions()]
                tquery_seq_trim=read.query_sequence
                pairs = np.array(read.get_aligned_pairs())
                # print(tquery_seq_trim)
                # print(read.get_reference_sequence())
                if any(op[0] == 1 for op in read.cigartuples):
                    if read.cigartuples[0][0] == 4:
                        pairs = pairs[read.cigartuples[0][1]:]
                    if read.cigartuples[-1][0] == 4:
                        pairs = pairs[:-read.cigartuples[-1][1]]
                    mask = (pairs[:, 0] != None) & (pairs[:, 1] == None)
                    mask_indices = np.where(mask)[0]
                    
                    splits = np.where(np.diff(mask_indices) > 1)[0] + 1
                    insertion_groups = np.split(mask_indices, splits)
                                           
                    tquery_array = np.array(list(tquery_seq_trim))
                    for group in insertion_groups:
                        mask_indices_pos = pairs[group, 0]
                        valid_indices = mask_indices_pos[mask_indices_pos != None].astype(int)
                        
                        extracted_nucleotides = "".join(tquery_array[valid_indices])
                        if (max(group)-min(group)) > insertion_thred:
                            continue
                        #取insertion上游第一个碱基的位置
                        tindex = pairs[min(group)-1,1] - start_pos
                        insInfoF[tindex].append(extracted_nucleotides)
                    
                    mask_indices_pos = pairs[mask_indices, 0]
                    valid_indices = mask_indices_pos[mask_indices_pos != None].astype(int)
                    mask = np.ones(len(tquery_array), dtype=bool)
                    mask[valid_indices] = False
                    tquery_array = tquery_array[mask]
                    tquery_seq_trim = ''.join(tquery_array)
                    if read.cigartuples[0][0] == 4:
                        tquery_seq_trim = tquery_seq_trim[read.cigartuples[0][1]:]
                    if read.cigartuples[-1][0] == 4:
                        tquery_seq_trim = tquery_seq_trim[:-read.cigartuples[-1][1]]
                else:
                    if read.cigartuples[0][0] == 4:
                        pairs = pairs[read.cigartuples[0][1]:]
                        tquery_seq_trim = tquery_seq_trim[read.cigartuples[0][1]:]
                    if read.cigartuples[-1][0] == 4:
                        pairs = pairs[:-read.cigartuples[-1][1]]
                        tquery_seq_trim = tquery_seq_trim[:-read.cigartuples[-1][1]]
                # print(tquery_seq_trim)
                tref = read.get_reference_sequence()
                mdtag = read.get_tag("MD")
                tread = tquery_seq_trim
                ttcov = []
                del_pos = []
                del_pattern = []
                mismatch_pos = []
                mismatch_pattern = []
                
                for match in re.finditer(r'(\d+)|(\^[A-Z]+)|([A-Z])', mdtag):
                    if match.group(1):
                        ttcov.extend(tpos[:int(match.group(1))])
                        tpos = tpos[int(match.group(1)):]
                        tread = tread[int(match.group(1)):]
                        tref = tref[int(match.group(1)):]
                    elif match.group(2): 
                        deletion = match.group(2)[1:]
                        if len(deletion) <= del_thred:
                            ttcov.append(tpos[0]-1)
                            del_pos.append(tpos[0]-1)
                            del_pattern.append(deletion)
                            if ForR == 1:
                                delInfoF[tpos[0]-1].append(deletion)
                            else:
                                delInfoR[tpos[0]-1].append(deletion)
                        tref = tref[len(deletion):]
                    elif match.group(3):
                        mut_pattern = tref[0].upper()+"->"+tread[0].upper()
                        mismatch_pos.append(tpos[0])
                        mismatch_pattern.append(mut_pattern)
                        ttcov.append(tpos[0])
                        if ForR == 1:
                            mutInfoF[tpos[0]].append(mut_pattern)
                            mutCountF[tpos[0]] += 1
                        else:
                            mutInfoR[tpos[0]].append(mut_pattern)
                            mutCountR[tpos[0]] += 1
                        tpos = tpos[1:]
                        tread = tread[1:]
                        tref = tref[1:]

                if ForR == 1:
                    coverageF[ttcov] += 1
                else:
                    coverageR[ttcov] += 1

                # print("Insertion nucleotides:", insert_nucleotides)
                # print("Insertion positions:", insert_pos)
                # print("del_pos:", del_pos)
                # print("mismatch_pos:", mismatch_pos)
                # print("mismatch_pattern:", mismatch_pattern)
                # print(tref)
                # print(tpos)
                # input()
            
            for i in range(len(mutCountF)):
                if not args.all and str(coverageF[i]) == "0" and str(coverageR[i]) == "0":
                    continue
                out.write(",".join([
                    tgeno["SN"],
                    str(i + start_pos + 1),
                    str(coverageF[i]),
                    "0" if not mutInfoF[i] else str(mutCountF[i])+":"+";".join(mutInfoF[i]),
                    "0" if not delInfoF[i] else str(len(delInfoF[i]))+":"+";".join(delInfoF[i]),
                    "0" if not insInfoF[i] else str(len(insInfoF[i]))+":"+";".join(insInfoF[i]),
                    str(coverageR[i]),
                    "0" if not mutInfoR[i] else str(mutCountR[i])+":"+";".join(mutInfoR[i]),
                    "0" if not delInfoR[i] else str(len(delInfoR[i]))+":"+";".join(delInfoR[i]),
                    "0" if not insInfoR[i] else str(len(insInfoR[i]))+":"+";".join(insInfoR[i]),
                ])+"\n")
         
    bamfile.close()

def run_DNA_MaP(input_file, output):
    pass

def check_bam(input_file):
    import os,sys
    import subprocess

    # Check if file exists
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} does not exist")
        sys.exit(1)
    else:
        logging.info(f"Found input file {input_file}")

    # Check if file is BAM format 
    if not input_file.endswith('.bam'):
        logging.error(f"Input file {input_file} is not a BAM file")
        sys.exit(1)
    else:
        logging.info("Input file is in BAM format")
        
    # Check if BAM index exists
    if not os.path.exists(input_file + '.bai'):
        logging.error(f"BAM index file {input_file}.bai does not exist. Please index the BAM file first")
        sys.exit(1)
    else:
        logging.info("Found BAM index file")

    
    cov_file = input_file + '.cov'
    if not os.path.exists(cov_file):
        logging.info(f"Coverage file {cov_file} not found, generating it now...")
        try:
            with pysam.AlignmentFile(input_file, "rb") as bam, open(cov_file, 'w') as out:
                for ref in bam.references:
                    for column in bam.pileup(ref, min_base_quality=0):
                        out.write(f"{ref}\t{column.pos+1}\t{column.n}\n")
            logging.info(f"Successfully generated coverage file: {cov_file}")
        except Exception as e:
            logging.error(f"Failed to generate coverage file: {e}")
            sys.exit(1)
    else:
        logging.info("Found coverage file")

    return True


def testMaP(input_file, dna):
    import pysam
    limit=10000
    bamfile = pysam.AlignmentFile(input_file, "rb")
    logging.info("First 5 references in BAM header:")
    for i, sq in enumerate(bamfile.header["SQ"][:5]):
        logging.info(f"  {sq['SN']}: {sq['LN']:,} bp")
    total_reads = 0
    mapped_reads = 0
    single_end = 0
    paired_end = 0
    r1_forward = 0
    r1_reverse = 0 
    r2_forward = 0
    r2_reverse = 0
    
    forward_count=0
    reverse_count=0

    for read in bamfile.fetch():
        total_reads += 1
        if total_reads > limit:
            break
        ForR = -1 if read.is_read1 == read.is_reverse else 1
        if ForR == 1:
            forward_count += 1
        else:
            reverse_count += 1
        
        if not read.is_unmapped:
            mapped_reads += 1
            
            if read.is_paired:
                paired_end += 1
                if read.is_read1:
                    if read.is_reverse:
                        r1_reverse += 1
                    else:
                        r1_forward += 1
                else:
                    if read.is_reverse:
                        r2_reverse += 1
                    else:
                        r2_forward += 1
            else:
                single_end += 1
    
    logging.info("=" * 60)
    logging.info("Read Statistics:")
    logging.info("-" * 60)
    logging.info(f"Total reads sampled:     {total_reads:>10,d}")
    logging.info(f"Mapped reads:            {mapped_reads:>10,d}  ({mapped_reads/total_reads*100:>6.1f}%)")
    logging.info(f"Single-end reads:        {single_end:>10,d}  ({single_end/total_reads*100:>6.1f}%)")
    logging.info(f"Paired-end reads:        {paired_end:>10,d}  ({paired_end/total_reads*100:>6.1f}%)")
    logging.info(f"Forward reads:           {forward_count:>10,d}  ({forward_count/total_reads*100:>6.1f}%)")
    logging.info(f"Reverse reads:           {reverse_count:>10,d}  ({reverse_count/total_reads*100:>6.1f}%)")
    logging.info("-" * 60)
    if paired_end > 0:
        logging.info("Paired-end Read Mapping Statistics:")
        logging.info(f"Read 1 Forward:          {r1_forward:>10,d}  ({r1_forward/paired_end*100:>6.1f}%)")
        logging.info(f"Read 1 Reverse:          {r1_reverse:>10,d}  ({r1_reverse/paired_end*100:>6.1f}%)")
        logging.info(f"Read 2 Forward:          {r2_forward:>10,d}  ({r2_forward/paired_end*100:>6.1f}%)")
        logging.info(f"Read 2 Reverse:          {r2_reverse:>10,d}  ({r2_reverse/paired_end*100:>6.1f}%)")
    logging.info("=" * 60)
    bamfile.close()

def run(args):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    check_bam(args.input_bam)
    if args.output is None:
        args.output = os.path.basename(args.input_bam) + ".AtlasMaP"
    if args.test:
        testMaP(args.input_bam, args.dna)
    else:
        if args.dna:
            run_DNA_MaP(args.input_bam, args.output)
        else:
            run_RNA_MaP(args.input_bam, args.output, args.del_thred, args.insertion_thred, args.all, args.gz)

if __name__ == "__main__":

    logo=r'''      
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
This module is used to process DMS-MaP, SHAPE-MaP output to complete the 
first step after mapping, that is, to obtain mutation, insertion, deletion 
information. Process BAM files and output into AtlasMaP format.'''.format(logo)

    parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_bam', type=str, help='Path to the input indexed BAM file.')
    parser.add_argument('--dna', action='store_true', help='DNA model or RNA model, default is RNA model', default=False)
    parser.add_argument('--test', action='store_true', help='test mode, default is False', default=False)
    parser.add_argument('-o', '--output', type=str, default=None, help='output file name, default is input_file.AtlasMaP')
    parser.add_argument('-dt', '--del_thred', type=int, default=5, help='Deletions longer than this threshold will be ignored (default: 5)')
    parser.add_argument('-it', '--insertion_thred', type=int, default=5, help='Insertions longer than this threshold will be ignored (default: 5)')
    parser.add_argument('--all', action='store_true', help='output position with no coverage, default is False', default=False)
    parser.add_argument('-gz', action='store_true', help='output gzipped file, default is False', default=False)

    args = parser.parse_args()
    run(args)
