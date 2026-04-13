# Beth-2
Beth-2 is a Python toolkit for analyzing SARS-CoV-2 within-host viral mutation datasets.

It supports preprocessing of mutation call data, site-wise mutation statistics, and robust evaluation of mutation prediction models against gold-standard answer sets.

## Requirements
Python 3.7+
All dependencies can be installed with:
```
conda env create -f environment.yml
conda activate gisaid
```

## Quick Start

### 1. Preprocessing
Prepare raw mutation data.
```
python Scripts/Preprocessing.py \
  --input_dir vcf/ \
  --output_dir csv/ \
  --metadata_path Data/results_read_run_tsv.tsv \
  --start_pos 21563 \
  --end_pos 25384 \
  --region IND
```
#### Parameters:
•	--input_dir: Mutation vcf files directory
 
•	--output_dir: Mutation csv files directory

### 2. Hypothesis Testing
Compute site-specific mutation significance and prepare data for evaluation.
```
python Hypothesis.py --input_dir /path/to/Data --output_dir ../output
```
#### Parameters:
•	--input_dir: Mutation csv files directory
  
•	--output_dir: Where result tables are written
 
•	--rho: Depth threshold
 
•	--theta: Allele frequency threshold
 
•	--ns: Sample size threshold
 
•	--regions: Region labels to process
 
 
### 3. Evaluation
Compare predicted mutations to answer set and compute metrics.
```
python Scripts/Evaluation.py \
  --data_file Output/hypothesis_results.csv \
  --answer_file Data/answer_set_pr50.csv \
  --output_dir Output/ 
```
#### Parameters:
•	--data_file: Output from Hypothesis step
 
•	--answer_file:  Mutation set with over 50% population prevalence
 
•	--output_dir: Output directory
 
•	--p_value, --OR: Statistical thresholds

## Update

An update of BETH-2-predicted positive Spike mutations based on data available through February 2026 is provided in [`Output/reviewer_update_2026-02/`](./Output/reviewer_update_2026-02).


## Citation

If you use Beth2 for your research, please cite:

Su S, et al. Predicting Key Viral Mutations from Within-Host Deep Sequencing Data. (Unpublished Manuscript, 2026)

