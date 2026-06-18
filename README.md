# gecko_downstream_analysis

These scripts are for preparation of sgRNA/target read count data from raw .fastq files for genome-wide CRISPR screening experiments, with collection of some basic quality control statistics. The algorithm will only count sgRNAs which exactly match defined sgRNAs in a specific library.

### Requirements:
- Python 3
- pandas

- Sequence data in .fastq format from high-throughput sequencing (e.g. Illumina) of sgRNA protospacer regions, e.g. as described in Joung et al. (2017), Nature Protocols. This is designed for single-end sequencing - if paired end sequencing has been used, only R1 will be used here.
- sgRNA:target definition file, as a tab-delimited file with first column for 20-nt sgRNA sequence and second column for target (no header row)
Note that the script is currently set up for 20-nt sgRNAs only, but could be easily modified to allow for other lengths.

### Outputs:
- Tab-delimited file of count data for each cell pool, in a format suitable for downstream applications such as the MAGeCK pipeline:

  | seq | target | pool1_count | pool2_count | etc. |
  | --- | ---- | ---- | ---- | ---- |
  |CGAAGCGGGGCTTCTTACGT|Gene1|1221|657|...|
  |ACTTCGCCAGCTTACAAAAG|Gene1|156|12|...|

- A text file named 'pool_[pool name]_summary.txt' with basic summary statistics for read counts and library quality control parameters such as sgRNA coverage
- The script 'collate_sgrna_summary_stats.py' will also output a .csv with a single table collating the summary statistics for all cell pools (named with the above naming convention) in a single directory.

### Running instructions
To count sgRNA reads from fastq files, run the script 'count_gecko_sgrnas.py'. The following parameters need to be specified - these are currently hard-coded in global variables pending script refactoring, and so need to be modified in the script itself before running:
```
POOLS = {
         'pool1_name': '/path/to/pool1_fastqfile.fastq',
         'pool2_name': '/path/to/pool2_fastqfile.fastq'
        }
```
A dictionary of the names you want to assign to each cell pool in output files, and their associated fastq files.
Note - the script can support splitting of sequencing data for a single cell pool across separate files for up to 4 lanes, if the file names are distinguished by 'L001', 'L002' etc but are otherwise identical.

```
GUIDELIST = "/path/to/sgrna_seq_target_file.txt"
```
The sgRNA:target definition file in the format specified above

```
OUTPUTFILE = "/path/to/outputfile.txt"
```
Chosen file name for output count file

```
OUTPUT_NONMATCHING = False
```
If set to True, the script will output an additional file (prefixed 'errors.') containing the counts of all sequences (between the specified flanks sequences, see below) that don't match sgRNAs specified in the sgRNA:target definition file, including sequences that are the wrong length etc. This can be useful for troubleshooting, e.g. to rule out library contamination or in the case of unexpected poor quality control metrics. However it will substantially slow down the script and the resulting output can be huge, so this is not recommended for routine use.

```
KEY1 = "CGAAACACCG"
KEY2 = "GTTT"
```
These are the left and right flank sequences for the sgRNA protospacer sequence - the script is looking for 20nt sequences between these flanks, and will flag reads as 'wrong length' if the spacing is not 20nt. These flanks can be set at any length, according to the vector used, read length and tolerance for sequencing errors (longer sequences will lead to discarding more reads). The recommended KEY1 sequence CGAAACACCG above will work for many vectors. The KEY2 (right flank sequence) will depend particularly on the tracrRNA variant in the chosen vector - e.g. GTTTT will work well for vectors such as LentiCRISPRv2, while GTTT will work for vectors such as pkLV2.

```
SEQ_COLUMN = 1
GENE_COLUMN = 0
```
Columns (0-indexed) containing sequence / target data in the sgRNA:target file, to allow some flexibility to deviate from the format specified above.

```
MIN_QUAL = ','
```
ASCII character to denote the lowest acceptable Phred quality score in the fastq files. Protospacer sequences containing any base with a score lower than this will be discarded. It is usually recommended to set this score to a fairly permissive level (e.g. ',' as above) unless you have some highly-similar sgRNAs where there is a realistic possibility that read errors could lead to sgRNA misclassification.


After running counts, if you want to collate the quality control stats, run the script collate_sgrna_summary_stats.py as follows:

```
python collate_sgrna_summary_stats.py /path/to/target_directory /path/to/outputfile.csv
```
This script will read summary stats from every file in the target directory starting with the word 'pool' (so make sure this has not been used for any other files in the directory!). The output file will be in .csv format.

Statistics given for each pool are:
 - Total read count in .fastq files
 - Total and % perfectly matched reads (i.e. matching an sgRNA in the specified library)
 - Number and % of reads with 'wrong length' guides, i.e. not 20nt between specified flank sequences
 - Number and % of non-matching 20mer guides (i.e. guides of the correct length but not matching any sgRNA in the specified library - this can be a useful indicator of library contamination)
 - Guides found - total number of sgRNAs with at least 1 read in the cell pool
 - Expected library size - total number of sgRNAs in the supplied sgRNA:target library definition file
 - % guides found (i.e. guides in the library with at least 1 read)
 - Reads per sgRNA (Mean number of matched reads / sgRNA, i.e. total matched reads / library size - note this ignores reads without a perfectly matched sgRNA sequence, and is calculated based on total (expected) library size rather than the number of sgRNAs actually detected. Typically >500X is recommended for novel libraries, but >100X can be sufficient in many contexts)
 - Skew ratio. This is the ratio of the 90th:10th percentile read counts. This is calculated based on the full expected library size, so will not be calculable if library coverage is <90%.
   
## Stouffer metaanalysis of MAGeCK output

The script mageck_stouffer_metaanalysis.py enables meta-analysis of p-values from multiple independent CRISPR screening experiments, using results generated by the MAGeCK RRA method (Li et al. 2019). This is an alternative to the MAGeCK inbuilt methods of combining replicates by simple combination (combining read counts for technical replicates) or --paired analysis (which considers each sgRNA result for a gene as independent even when the sgRNAs are the same, i.e. 3 replicates with n sgRNAs per gene would be considered as 3n independent sgRNAs per gene). This method is particularly suitable where replicates share a common baseline or where data quality is variable across the experiments to be combined, and we have found that it can improve screen discriminatory performance compared to the inbuilt methods for some datasets.

This script uses the 'Stouffer' z-score method to combine p-values from different experiments. This method has been selected due to its ability to account for direction of effects, which is particularly appropriate for the one-sided tests in MAGeCK RRA, as well as the ability to introduce weightings. As a weighting measure, we use (1/standard deviation of non-targeting guide LFCs) as a measure of precision, which is conceptually similar to the sqrt(n) method commonly used.

The script requires:
 - Python 3, pandas, scipy, numpy and matplotlib
 - For each experiment to be combined:
   - A gene-level MAGeCK output file, tab-delimited, named <exp_name>.gene_summary.txt, with column headers including 'pos|lfc', 'pos|p-value', 'pos|fdr', 'neg|lfc', 'neg|p-value' and 'neg|fdr', with gene name as the first column.
   - An sgrna-level MAGeCK output file, tab-delimited, name <exp_name>.sgrna_summary.txt in the same folder, with column headers including 'Gene' and 'LFC'. This must include nontargeting guides (note - these may be absent if MAGeCK was run with the nontargeting guides specified as control sgRNAs for variance computations), identified by including the word 'nontargeting' in the gene name column.

To run the analysis:

Until I get around to coding the command line variables properly, the following parameters need to be hard-coded in the script:
- cwd: working directory (default "./")
- runs: a list of experiment names. These are the prefixes of the .gene_summary.txt and .sgrna_summary.txt files mentioned above, i.e. if your file is named 'experiment_vs_control.gene_summary.txt' just enter 'experiment_vs_control' in this list.
- DIRECTION: a list of directions to meta-analyse, i.e. whether to look at positive enrichment ('pos'), negative enrichment ('neg') or both
- OUTSTEM: prefix for the outputfile(s) - this will be appended with the direction of enrichment for the analysis + '_stouffer_metaanalysis.txt'

Then just run the script as
```
python mageck_stouffer_metaanalysis.py
```

The resulting output file will include collated results of individual screens, weighting factors, median LFC across experiments for each gene, the p-value for Stouffer meta-analysis, and the FDR (BH method) for the meta-analysis.
