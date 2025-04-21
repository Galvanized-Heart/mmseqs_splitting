# MMSeqs2 Sequence Clustering and Data Splitting Tool

This tool uses MMSeqs2 to perform sequence clustering and create train/validation/test splits for machine learning applications while avoiding data leakage between similar sequences.

## Installation

### Quick Installation (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install mmseqs2
```

### Manual Installation

For other platforms or if you need a specific version, follow the instructions at the [MMSeqs2 GitHub repository](https://github.com/soedinglab/MMseqs2).

#### Basic steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/soedinglab/MMseqs2.git
   cd MMseqs2
   ```

2. Compile and install:
   ```bash
   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
   make
   make install
   ```

## Dependencies

- Python 3.6+
- pandas
- matplotlib
- MMSeqs2

Install Python dependencies:
```bash
pip install pandas matplotlib
```

## Usage

### Basic Usage

```python
from mmseqs2_split import add_splits_to_csv, analyze_existing_splits

# Create splits based on sequence clustering
add_splits_to_csv(
    "input.csv",
    "output_with_splits.csv",
    seq_col="SEQ",
    seq_ids=[0.7],
    split_ratios={'train': 0.8, 'valid': 0.1, 'test': 0.1},
    coverage=0.7,
    cov_mode=2,
    cluster_mode=0,
    seed=42
)

# Analyze the resulting clusters
df = pd.read_csv("output_with_splits.csv")
stats = analyze_existing_splits(df, 'cluster_seqid70_cov70_cov_mode2_cluster_mode0_tr80_va10_ts10')
print(stats)
```

### Example Input and Output

**Input CSV:**
```
,CMPD_SMILES,SEQ,Ki
0,Nc1ncnc2c1ncn2[C@@H]1O[C@H](C[C@@H](N)CC[C@H](N)C(=O)O)[C@@H](O)[C@H]1O,SLSPAVQTFWKWLQEEGVITAKTPVKASVVTEGLGLVALKDISRNDVILQVPKRLWINPDAVAASEIGRVCSELKPWLSVILFLIRERSREDSVWKHYFGILPQETDSTIYWSEEELQELQGSQLLKTTVSVKEYVKNECLKLEQEIILPNKRLFPDPVTLDDFFWAFGILRSRAFSRLRNENLVVVPMADLINHSAGVTTEDHAYEVKGAAGLFSWDYLFSLKSPLSVKAGEQVYIQYDLNKSNAELALDYGFIEPNENRHAYTLTLEISESDPFFDDKLDVAESNGFAQTAYFDIFYNRTLPPGLLPYLRLVALGGTDAFLLESLFRDTIWGHLELSVSRDNEELLCKAVREACKSALAGYHTTIEQDRELKEGNLDSRLAIAVGIREGEKMVLQQIDGIFEQKELELDQLEYYQERRLKDLGLCGENGDILENLYFQ,4.3279021421
```

**Output CSV:**
```
,CMPD_SMILES,SEQ,Ki,cluster_seqid70_cov70_cov_mode2_cluster_mode0_tr80_va10_ts10,split_seqid70_cov70_cov_mode2_cluster_mode0_tr80_va10_ts10
0,Nc1ncnc2c1ncn2[C@@H]1O[C@H](C[C@@H](N)CC[C@H](N)C(=O)O)[C@@H](O)[C@H]1O,SLSPAVQTFWKWLQEEGVITAKTPVKASVVTEGLGLVALKDISRNDVILQVPKRLWINPDAVAASEIGRVCSELKPWLSVILFLIRERSREDSVWKHYFGILPQETDSTIYWSEEELQELQGSQLLKTTVSVKEYVKNECLKLEQEIILPNKRLFPDPVTLDDFFWAFGILRSRAFSRLRNENLVVVPMADLINHSAGVTTEDHAYEVKGAAGLFSWDYLFSLKSPLSVKAGEQVYIQYDLNKSNAELALDYGFIEPNENRHAYTLTLEISESDPFFDDKLDVAESNGFAQTAYFDIFYNRTLPPGLLPYLRLVALGGTDAFLLESLFRDTIWGHLELSVSRDNEELLCKAVREACKSALAGYHTTIEQDRELKEGNLDSRLAIAVGIREGEKMVLQQIDGIFEQKELELDQLEYYQERRLKDLGLCGENGDILENLYFQ,4.3279021421,16610,train
```

## Key Parameters

### Sequence Identity (`--min-seq-id`)
The minimum fraction of aligned residues required for two sequences to be considered similar. Values range from 0.0 to 1.0:
- 0.9-1.0: Very strict, nearly identical sequences
- 0.7-0.8: Moderate similarity, likely similar function
- 0.5-0.6: Remote homologs
- <0.5: Very distant relationships

### Coverage (`-c`)
The minimum fraction of sequence that must be aligned. Values range from 0.0 to 1.0:
- 0.8-1.0: Full-length homologs
- 0.5-0.7: Domain or motif similarity
- <0.5: Short motif detection

### Coverage Mode (`--cov-mode`)
Controls how coverage is calculated:
- 0: Bidirectional (both query and target coverage must meet threshold)
- 1: Target Coverage (useful for finding full-length homologs)
- 2: Query Coverage (useful for domain-level clustering)
- 3: Target-in-Query Coverage (useful for short motif detection)

### Cluster Mode (`--cluster-mode`)
Determines the clustering algorithm:
- 0: Greedy Set Cover (MMSeqs2 default) - Forms clusters to minimize the total number needed
- 1: Connected component - Uses transitivity to form larger clusters with remote homologs
- 2: Greedy Cluster (CD-HIT analog) - Length-sorted greedy approach

### Split Ratios
Controls how clusters are assigned to train/validation/test sets:
- Default: `{'train': 0.7, 'valid': 0.15, 'test': 0.15}`

## Recommended Parameter Combinations

### Application-Specific Settings

| Application | Parameters |
|-------------|------------|
| Domain-level clustering | `--cov-mode 2 -c 0.7` |
| Full-length orthologs | `--cov-mode 1 -c 0.9` |
| Short motif detection | `--cov-mode 3 -c 0.3` |

### For Kinetic Parameter Predictions

| Approach | Parameters | Use case |
|----------|------------|----------|
| Conservative | `-c 0.8 --cov-mode 2` | Full-Length Homologs |
| Balanced | `-c 0.6 --cov-mode 5` | Functional Motifs |
| Lenient | `-c 0.4 --cov-mode 3` | Broad Similarity |

## Analyzing Cluster Quality

Use the `analyze_existing_splits` function to evaluate clustering quality:

```python
stats = analyze_existing_splits(df, 'cluster_seqid70_cov70_cov_mode2_cluster_mode0_tr80_va10_ts10')
```

This provides:
- Number of clusters
- Average cluster size
- Number of singleton clusters
- Maximum cluster size
- A histogram of cluster sizes

## References

- MMSeqs2: https://github.com/soedinglab/MMseqs2
- Steinegger M, Söding J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35, 1026–1028. doi: 10.1038/nbt.3988
