import numpy as np
import pysam
import re
import os
import subprocess
import logging
import tempfile
import shutil
from datetime import datetime
from collections import defaultdict
import random
import string

def read_fasta(fasta_file):
    """Read a FASTA file and return a dict of {name: sequence}."""
    seqs = {}
    name = ""
    seq_parts = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq_parts).upper().replace('U', 'T')
                name = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name:
        seqs[name] = "".join(seq_parts).upper().replace('U', 'T')
    return seqs


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


def fix_bam_tags_DF_DC(bamfile, output_bam, filter_ga=True, reference=None):
    """
    Add DF (reference sequence) and DC (constraint string) tags to BAM reads.

    Uses aligned pairs + reference FASTA when provided (required for Nanopore BAMs
    without MD tags). Falls back to MD tag parsing when reference is not provided.

    Args:
        bamfile (str): Input BAM file
        output_bam (str): Output BAM file
        filter_ga (bool): Filter G->A mutations (DMS artifact). Default True.
        reference (str): Path to reference FASTA file. Required for BAMs without MD tags.
    """
    # Load reference sequences if provided
    ref_seqs = {}
    if reference:
        ref_seqs = read_fasta(reference)

    with pysam.AlignmentFile(bamfile, "rb") as in_bam:
        header = in_bam.header.copy()
        with pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
            for read in in_bam:
                if read.is_unmapped:
                    continue
                # Skip supplementary and secondary alignments
                if (read.flag & 2048) or (read.flag & 256):
                    continue

                ref_name = read.reference_name
                seq_len = read.reference_length  # aligned reference length for this read

                # Determine if we should use aligned-pairs approach
                has_md = read.has_tag("MD")

                if reference and ref_name in ref_seqs:
                    # --- Aligned-pairs approach (works for all BAMs) ---
                    full_ref_seq = ref_seqs[ref_name]
                    ref_len = len(full_ref_seq)
                    read.set_tag('DF', full_ref_seq)

                    constraint = ['0'] * ref_len

                    for query_pos, ref_pos in read.get_aligned_pairs():
                        if ref_pos is None or query_pos is None:
                            continue
                        if ref_pos < 0 or ref_pos >= ref_len:
                            continue
                        read_base = read.query_sequence[query_pos].upper()
                        ref_base = full_ref_seq[ref_pos].upper()

                        if ref_base == read_base:
                            # Match — constraint stays '0' (unconstrained/wild-type)
                            constraint[ref_pos] = '0'
                        else:
                            # G->A artifact filter
                            if filter_ga and ref_base == 'G' and read_base == 'A':
                                pass  # skip, leave as '0'
                            else:
                                constraint[ref_pos] = '1'  # mutation

                    constraint = "".join(constraint)
                    read.set_tag('DC', constraint)
                    out_bam.write(read)

                elif has_md:
                    # --- MD tag approach (legacy, for BAMs with MD tags) ---
                    ref_seq = read.get_reference_sequence()
                    read.set_tag('DF', ref_seq)

                    tquery_seq, pairs = remove_soft_clipping(read)
                    constraint = ['0'] * len(ref_seq)

                    mdtag = read.get_tag("MD")
                    tread = tquery_seq
                    tref = ref_seq
                    tpos = [p - read.reference_start for p in read.get_reference_positions()]

                    for match in re.finditer(r'(\d+)|(\^[A-Z]+)|([A-Z])', mdtag):
                        if match.group(1):  # Match
                            match_len = int(match.group(1))
                            tpos = tpos[match_len:]
                            tread = tread[match_len:]
                            tref = tref[match_len:]
                        elif match.group(2):  # Deletion
                            deletion = match.group(2)[1:]
                            pos = tpos[0] - 1
                            if 0 <= pos < len(constraint):
                                constraint[pos] = '1'
                            tpos = [p + len(deletion) for p in tpos]
                            tref = tref[len(deletion):]
                        elif match.group(3):  # Mutation
                            ref_base = match.group(3).upper()
                            pos = tpos[0]
                            if filter_ga and ref_base == 'G' and len(tread) > 0 and tread[0].upper() == 'A':
                                tpos = tpos[1:]
                                tread = tread[1:]
                                tref = tref[1:]
                                continue
                            if 0 <= pos < len(constraint):
                                constraint[pos] = '1'
                            tpos = tpos[1:]
                            tread = tread[1:]
                            tref = tref[1:]

                    if any(op[0] == 1 for op in read.cigartuples):
                        mask = (pairs[:, 0] != None) & (pairs[:, 1] == None)
                        mask_indices = np.where(mask)[0]
                        if len(mask_indices) > 0:
                            splits = np.where(np.diff(mask_indices) > 1)[0] + 1
                            insertion_groups = np.split(mask_indices, splits)
                            for group in insertion_groups:
                                if len(group) > 0:
                                    pos = pairs[min(group)-1, 1] - read.reference_start
                                    if pos >= 0 and pos < len(constraint):
                                        constraint[pos] = '1'

                    constraint = "".join(constraint)
                    read.set_tag('DC', constraint)
                    out_bam.write(read)
                else:
                    logging.warning(f"Read {read.query_name}: no MD tag and no reference provided, skipping")
                    continue

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
            # Skip supplementary and secondary alignments
            if (read.flag & 2048) or (read.flag & 256):
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
    
    cmd = f"cat {id_file} | xargs -P {num_cores} -I {{}} python3 {fold_script} {{}} {temp_dir}"
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


def forgi_vectorize_and_cluster(bamfile, output_dir, gamma=6, num_clusters=3):
    """
    Extract unique mutation profiles from BAM, fold each with CONTRAfold,
    vectorize with forgi, then run PCA + K-means clustering.

    Args:
        bamfile (str): Path to the folded BAM file (must have DF and DC tags)
        output_dir (str): Directory to write output files
        gamma (float): CONTRAfold gamma parameter
        num_clusters (int): Number of K-means clusters

    Returns:
        str: Path to the clusters CSV file, or None on failure
    """
    import pandas as pd
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    from sklearn.metrics import pairwise_distances_argmin_min
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)

    # --- Step 1: Extract unique DC bit profiles per reference ---
    # Group reads by reference name
    ref_profiles = defaultdict(lambda: defaultdict(int))
    ref_seqs = {}

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            if (read.flag & 2048) or (read.flag & 256):
                continue
            try:
                dc_tag = read.get_tag('DC')
                df_tag = read.get_tag('DF')
            except KeyError:
                continue
            ref_name = read.reference_name
            if ref_name not in ref_seqs:
                ref_seqs[ref_name] = df_tag.upper().replace('T', 'U')
            ref_profiles[ref_name][dc_tag] += 1

    if not ref_profiles:
        logging.warning("No reads with DC tags found")
        return None

    all_cluster_files = []

    for ref_name, profiles in ref_profiles.items():
        logging.info(f"Processing reference: {ref_name} ({len(profiles)} unique profiles)")
        ref_seq = ref_seqs[ref_name]
        seq_len = len(ref_seq)

        # Convert DC tags to bit-style profiles: '1' in DC means mutation
        unique_bits = {}
        for dc_tag, count in profiles.items():
            bits = dc_tag  # DC tag is already '0'/'1' string
            if len(bits) != seq_len:
                continue
            unique_bits[bits] = count

        if not unique_bits:
            logging.warning(f"  No valid profiles for {ref_name}, skipping")
            continue

        # --- Step 2: CONTRAfold fold + forgi vectorize each unique profile ---
        tmpdir = tempfile.mkdtemp(prefix='forgi_')
        bpseq_path = os.path.join(tmpdir, 'c.bpseq')
        fold_path = os.path.join(tmpdir, 'f.fold')
        post_path = os.path.join(tmpdir, 'p.txt')
        db_path = os.path.join(tmpdir, 'f.db')
        txt_path = os.path.join(tmpdir, 'f.txt')

        fold2db = os.path.join(os.path.dirname(__file__), '..', '..', 'tests',
                               'fold2dotbracketFasta.py')
        if not os.path.exists(fold2db):
            # Try to find it in PATH
            fold2db = shutil.which('fold2dotbracketFasta.py') or 'fold2dotbracketFasta.py'

        sizer_file = os.path.join(output_dir, f"{ref_name}_sizer.tab")
        forgi_file = os.path.join(output_dir, f"{ref_name}_forgi.tab")

        seen = defaultdict(int)
        processed = 0
        failed = 0

        sorted_profiles = sorted(unique_bits.items(), key=lambda x: x[1], reverse=True)

        with open(sizer_file, 'w') as outS, open(forgi_file, 'w') as serF:
            for t, (profile, count) in enumerate(sorted_profiles, 1):
                # Build constraint: '1' in DC (mutation) -> constraint '0' (unpaired)
                #                    '0' in DC (wild-type) -> constraint '-1' (unconstrained)
                state = []
                for pos, (base, bit) in enumerate(zip(ref_seq, profile), 1):
                    if bit == '1':
                        state.append((str(pos), base, '0'))
                    else:
                        state.append((str(pos), base, '-1'))

                ticks_string = ''.join([x[2] for x in state])
                seen[ticks_string] += 1
                if seen[ticks_string] > 1:
                    continue

                bit_prefix = f'bit_{t}'

                with open(bpseq_path, 'w') as f:
                    for i, base, tick in state:
                        f.write(f'{i}\t{base}\t{tick}\n')

                try:
                    result = subprocess.run(
                        ['contrafold', 'predict', '--constraints', bpseq_path,
                         '--parens', fold_path, '--posteriors', '0.0', post_path,
                         '--gamma', str(gamma)],
                        capture_output=True, text=True, timeout=120
                    )
                    if result.returncode != 0:
                        failed += 1
                        continue

                    subprocess.run(
                        ['python3', fold2db, '--input_file', fold_path,
                         '--tag', bit_prefix, '--output_file', db_path],
                        capture_output=True, text=True, timeout=120
                    )

                    subprocess.run(
                        ['rnaConvert.py', db_path, '-T', 'element_string',
                         '--force', '--to-file', '--filename', txt_path],
                        capture_output=True, text=True, timeout=120
                    )

                    es_file = txt_path + '001.element_string'
                    with open(es_file) as f:
                        for line_num, line in enumerate(f, 1):
                            if line_num == 3:
                                serF.write(f'{bit_prefix}\t{line.strip()}\n')
                except (FileNotFoundError, subprocess.TimeoutExpired) as e:
                    failed += 1
                    continue

                outS.write(f'{bit_prefix}\t{count}\t{profile}\n')
                processed += 1

        shutil.rmtree(tmpdir, ignore_errors=True)

        logging.info(f"  Folded {processed} unique structures ({failed} failed)")

        if processed < 3:
            logging.warning(f"  Too few structures ({processed}) for PCA/clustering, skipping {ref_name}")
            continue

        # --- Step 3: PCA on forgi vectors ---
        data = []
        headers = []
        with open(forgi_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                head, digits = parts[0], parts[1]
                vec = [float(d) for d in list(digits)]
                data.append(vec)
                headers.append(head)

        if len(data) < 3:
            logging.warning(f"  Too few forgi vectors ({len(data)}) for PCA, skipping {ref_name}")
            continue

        x = np.asarray(data, dtype='float64')
        n_components = min(2, x.shape[0], x.shape[1])
        pca = PCA(n_components=n_components)
        principal_components = pca.fit_transform(x)

        logging.info(f"  PCA explained variance ratio: {pca.explained_variance_ratio_.tolist()}")

        pca_file = os.path.join(output_dir, f"{ref_name}_pca.csv")
        pca_df = pd.DataFrame(
            data=principal_components,
            columns=[f'PC{i+1}' for i in range(n_components)]
        )
        pca_df['structure'] = headers
        pca_df.to_csv(pca_file, index=False)

        # PCA plot
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(principal_components[:, 0], principal_components[:, 1],
                   c='#ff7f00', s=10.0)
        ax.set_xlabel('Principal Component 1', fontsize=15)
        ax.set_ylabel('Principal Component 2', fontsize=15)
        ax.set_title(f'PCA on forgi vectors - {ref_name}', fontsize=15)
        plt.savefig(os.path.join(output_dir, f"{ref_name}_pca.png"), dpi=600)
        plt.savefig(os.path.join(output_dir, f"{ref_name}_pca.pdf"))
        plt.close(fig)

        # --- Step 4: K-means clustering ---
        actual_clusters = min(num_clusters, len(data))
        df = pd.read_csv(pca_file, index_col='structure')

        kmeans = KMeans(n_clusters=actual_clusters, random_state=0, n_init=10).fit(df)
        closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, df)

        cluster_colors = ['#4daf4a', '#ff7f00', '#f781bf', '#377eb8', '#e41a1c',
                          '#984ea3', '#a65628', '#999999', '#ffff33', '#dede00']

        for i in range(actual_clusters):
            count = kmeans.labels_.tolist().count(i)
            logging.info(f"  Cluster {i+1}: {count} structures")

        for k, i in enumerate(closest):
            logging.info(f"  Representative cluster {k+1}: {df.index[i]}")

        # Cluster plot
        fig, ax = plt.subplots(figsize=(6, 6))
        colors = [cluster_colors[k % len(cluster_colors)] for k in kmeans.labels_]
        ax.scatter(df['PC1'], df['PC2'], c=colors, s=50, alpha=0.50)
        for k, i in enumerate(closest):
            ax.scatter(df.iloc[i]['PC1'], df.iloc[i]['PC2'],
                       c=cluster_colors[k % len(cluster_colors)],
                       s=50, edgecolors='#e41a1c')
            ax.text(df.iloc[i]['PC1'] + 1.5, df.iloc[i]['PC2'] + 1.5,
                    df.index[i], ha='left', color='black', fontsize=8)
        ax.set_xlabel('PC 1', fontsize=15)
        ax.set_ylabel('PC 2', fontsize=15)
        ax.set_title(f'K-means clustering - {ref_name}', fontsize=15)
        plt.savefig(os.path.join(output_dir, f"{ref_name}_clusters.png"), dpi=600)
        plt.savefig(os.path.join(output_dir, f"{ref_name}_clusters.pdf"))
        plt.close(fig)

        # Cluster CSV output
        color_names = ['green', 'orange', 'pink', 'blue', 'red',
                       'purple', 'brown', 'grey', 'yellow', 'olive']
        df['cluster_id'] = kmeans.labels_.tolist()
        df['cluster_name'] = [str(k + 1) for k in kmeans.labels_.tolist()]
        df['cluster_colour'] = [color_names[k % len(color_names)] for k in kmeans.labels_.tolist()]

        clusters_csv = os.path.join(output_dir, f"{ref_name}_clusters.csv")
        df.to_csv(clusters_csv, index=True)
        all_cluster_files.append(clusters_csv)

        logging.info(f"  Clustering results saved to {clusters_csv}")

    if all_cluster_files:
        return all_cluster_files[0]
    return None
