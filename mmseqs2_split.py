import pandas as pd
import matplotlib.pyplot  as plt
import os
import subprocess
import random
import tempfile

def create_fasta_from_csv(csv_path, seq_col, fasta_path):
    """Create FASTA file from CSV using index as headers"""
    df = pd.read_csv(csv_path)
    with open(fasta_path, 'w') as f:
        for idx, row in df.iterrows():
            f.write(f">{idx}\n{row[seq_col]}\n")

def cluster_and_split(fasta_path, seq_id, split_ratios, coverage=0.8, cov_mode=1, threads=4, cluster_mode=0, seed=None,):
    """Cluster sequences and return split assignments with statistics"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # mmseqs2 clustering with configurable parameters
        cmds = [
            ["mmseqs", "createdb", fasta_path, f"{tmpdir}/DB"],
            ["mmseqs", "cluster", f"{tmpdir}/DB", f"{tmpdir}/clu", f"{tmpdir}/tmp",
             "--min-seq-id", str(seq_id), 
             "-c", str(coverage), 
             "--cov-mode", str(cov_mode), 
             "--threads", str(threads),
             "--cluster-mode", str(cluster_mode)],
            ["mmseqs", "createtsv", f"{tmpdir}/DB", f"{tmpdir}/DB", f"{tmpdir}/clu", f"{tmpdir}/clusters.tsv"]
        ]
        
        # Run mmseqs2 commands
        for cmd in cmds:
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"MMseqs2 failed: {e.stderr}")
                raise
        
        # Read cluster assignments and calculate statistics
        clusters = {}
        with open(f"{tmpdir}/clusters.tsv") as f:
            for line in f:
                rep, member = line.strip().split('\t')
                clusters.setdefault(rep, []).append(member)

        # Calculate cluster statistics
        num_clusters = len(clusters)
        total_members = sum(len(members) for members in clusters.values())
        average_size = total_members / num_clusters if num_clusters > 0 else 0

        # Split clusters with reproducible shuffling
        cluster_list = list(clusters.keys())
        rng = random.Random(seed)
        rng.shuffle(cluster_list)
        
        # Calculate split boundaries
        n = len(cluster_list)
        split_order = ['train', 'valid', 'test']
        ratios = [split_ratios[split] for split in split_order]
        
        train_end = int(ratios[0] * n)
        valid_end = train_end + int(ratios[1] * n)
        
        splits = {
            'train': cluster_list[:train_end],
            'valid': cluster_list[train_end:valid_end],
            'test': cluster_list[valid_end:]
        }

        # Create mapping from sequence ID to split
        cluster_id_map = {}
        split_map = {}
        for split_name, reps in splits.items():
            for rep in reps:
                for member in clusters[rep]:
                    cluster_id_map[member] = rep
                    split_map[member] = split_name

        stats = {
            'num_clusters': num_clusters,
            'average_cluster_size': average_size,
            'coverage': coverage,
            'cov_mode': cov_mode,
            'cluster_mode': cluster_mode,
            'split_train': split_ratios['train'],
            'split_valid': split_ratios['valid'],
            'split_test': split_ratios['test']
        }

        return split_map, cluster_id_map, stats

def add_splits_to_csv(input_csv_path, output_csv_path=None, seq_col="SEQ", seq_ids=[0.9], coverage=0.8, 
                      cov_mode=0, cluster_mode=0, threads=4, seed=None, split_ratios=None):
    """Add split columns with fully configurable parameters"""

    df = pd.read_csv(input_csv_path)

    # Set default output path
    if output_csv_path is None:
        output_csv_path = input_csv_path

    # Set default split ratios
    if split_ratios is None:
        split_ratios = {'train': 0.7, 'valid': 0.15, 'test': 0.15}
    
    # Validate split ratios
    required_keys = {'train', 'valid', 'test'}
    if not required_keys.issubset(split_ratios.keys()):
        missing = required_keys - split_ratios.keys()
        raise ValueError(f"Missing split ratio keys: {missing}")
    if abs(sum(split_ratios.values()) - 1.0) > 1e-6:
        raise ValueError(f"Split ratios must sum to 1.0, got {sum(split_ratios.values())}")

    # Initialize statistics log if file doesn't exist
    stats_path = "clustering_stats.csv"
    if not os.path.exists(stats_path):
        with open(stats_path, 'w') as f:
            f.write("seq_id,coverage,cov_mode,cluster_mode,train%,valid%,test%,num_clusters,avg_cluster_size\n")

    # Cluster sequences
    for seq_id in seq_ids:
        print(f"Processing seq_id: {seq_id}")
        # Store mmseqs DB and clu files in a temp file
        with tempfile.NamedTemporaryFile(mode='w') as fasta_file:
            # Create fasta to be used by mmseqs
            create_fasta_from_csv(input_csv_path, seq_col, fasta_file.name)

            # Get split and cluster mappings for data points in csv, and clustering stats
            split_map, cluster_map, stats = cluster_and_split(
                                                                fasta_file.name, 
                                                                seq_id,
                                                                split_ratios,
                                                                seed=seed,
                                                                coverage=coverage,
                                                                cov_mode=cov_mode,
                                                                cluster_mode=cluster_mode,
                                                                threads=threads
                                                            )
            
            # Create columns for cluster ID and split set
            col_name = (
                f"_seqid{int(seq_id*100)}_"
                f"cov{int(coverage*100)}_"
                f"cov_mode{cov_mode}_"
                f"cluster_mode{cluster_mode}_"
                f"tr{int(split_ratios['train']*100)}_va{int(split_ratios['valid']*100)}_ts{int(split_ratios['test']*100)}"
            )
            df["cluster"+col_name] = df.index.astype(str).map(cluster_map)
            df["split"+col_name] = df.index.astype(str).map(split_map)

            # Write clustering stats
            with open(stats_path, 'a') as f:
                f.write(
                    f"{seq_id},{coverage},{cov_mode},{cluster_mode}"
                    f"{split_ratios['train']},{split_ratios['valid']},{split_ratios['test']},"
                    f"{stats['num_clusters']},{stats['average_cluster_size']}\n"
                )

    # Save csv with cluster mappings and splits
    df.to_csv(output_csv_path, index=False)
    return df

def analyze_existing_splits(df, cluster_column):
    """Calculate cluster statistics from existing split column"""
    # Reconstruct cluster assignments
    clusters = df.groupby(cluster_column).apply(lambda x: list(x.index.astype(str)))
    
    stats = {
        'num_clusters': len(clusters),
        'avg_cluster_size': clusters.apply(len).mean(),
        'singleton_clusters': sum(1 for c in clusters if len(c) == 1),
        'max_cluster_size': clusters.apply(len).max()
    }
    
    # Plot cluster size distribution
    plt.figure()
    pd.Series(clusters.apply(len)).hist(bins=30)
    plt.title(f"Cluster Sizes - {cluster_column}")
    plt.show()
    
    return stats




# Example usage with custom splits
"""
add_splits_to_csv(
    "input.csv",
    "output_with_splits.csv",
    seq_col="SEQ",
    seq_ids=[0.9],
    split_ratios={'train': 0.8, 'valid': 0.1, 'test': 0.1},
    coverage=0.8,
    cov_mode=0,
    cluster_mode=0
    seed=42
)

add_splits_to_csv(
    "input.csv",
    seq_col="SEQ",
    seq_ids=[0.7],
    split_ratios={'train': 0.8, 'valid': 0.1, 'test': 0.1},
    coverage=0.7,
    cov_mode=2,
    cluster_mode=2
    seed=42
)

df = pd.read_csv("output_with_splits.csv")

stats = analyze_existing_splits(df, 'cluster_seqid90_cov80_mode0_tr80_va10_ts10')
print(stats)
"""




# Info about mmseqs2 params
"""
SEQUENCE IDENTITY THRESHOLD
---------------------------
The minimum fraction of aligned residues required for two sequences to cluster together
"""

"""
COVERAGE THRESHOLD
------------------
The minimum span of sequence required for two sequences to cluster together
"""

""" 
COVERAGE TYPES
--------------
- A query is the cluster representative
- A target is a target member for the cluster
--------------
Mode	Name	                    Condition Checked Against Threshold	        
0	    Bidirectional	        	(Aligned/QueryLen AND Aligned/TargetLen) >= -c	    
1	    Target Coverage	        	Aligned/TargetLen >= -c	                            
2	    Query Coverage	        	Aligned/QueryLen >= -c	                            
3	    Target-in-Query Coverage    TargetLen >= QueryLen * -c
4       ?                           ?
5       ?                           ?
"""

"""
CLUSTER MODES
-------------
Mode    Algorithm                       Description
0       Greedy Set Cover [Default]      Form a cluster (seq with most alignment is representative for cluster). 
                                        Subsequently, form a new cluster where a new representative is chosen 
                                        to meet threshold requirements in as few clusters as possible.
1       Connected component             Uses transitivity of the relations to form larger clusters with more 
                                        remote homologies. All sequences that are reachable in a breadth first
                                        search starting at the representative with the most connections are 
                                        added to a cluster.
2       Greedy Cluster (CD-HIT analog)  Sort sequences by length, and add all sequences that meet threshold in 
                                        a cluster. Create a new cluster with longest representative with the 
                                        remaining pool of sequences until all sequences are clustered.
"""

""" 
COMMON PARAM COMBOS
-------------------
Application	            Parameters
Domain-level clustering	--cov-mode 2 -c 0.7
Full-length orthologs	--cov-mode 1 -c 0.9
Short motif detection	--cov-mode 3 -c 0.3 
"""

"""
PARAM COMBOS SPECIFICALLY MOTIVATED FOR KINETIC PARAMETER PREDS
---------------------------------------------------------------
Application                                     Parameters
Conservative Approach (Full-Length Homologs)    -c 0.8 --cov-mode 2
Balanced Approach (Functional Motifs)           -c 0.6 --cov-mode 5
Lenient Approach (Broad Similarity)             -c 0.4 --cov-mode 3
"""

"""
VERIFYING DATA SPLIT
--------------------
Analyze cluster sizes:
import pandas as pd
clusters = pd.read_csv("clusters.tsv", sep="\t", header=None)
print(clusters[0].value_counts().describe())  # Cluster size distribution
- Aim for median cluster size 3-10 (adjust coverage if clusters are too large/small)

Compare Ki/Kd values within vs between clusters:
within_cluster_var = df.groupby("cluster_id")["Ki"].var().mean()
between_cluster_var = df["Ki"].var()
print(f"Variance reduction: {1 - within_cluster_var/between_cluster_var:.1%}")
- Good clusters reduce kinetic variance within clusters by >30%
"""

