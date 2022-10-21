# CLEMENT
- Genomic decomposition and reconstruction of **non-tumor** diploid subclones (2022)
- CLonal decomposition via Expectation-Maximization algorithm established in Non-Tumor setting
- Support multiple diploid sample

## Overview of CLEMENT workflow and core algorithms
<br/>

![Figure1 overview](https://user-images.githubusercontent.com/56012432/195979886-cd29df09-8291-4150-9001-db7dde5e7567.png)
<br/>

## Installation
### Dependencies
- python 3.6.x
- matplotlib 3.5.2
- seaborn 0.11.2
- numpy 1.21.5
- pandas 1.3.4
- scikit-learn 1.0.2
- scipy 1.7.3
- palettable 3.3.0

### Install from github (requires Python 3.6.* or newer)
1. https://github.com/Yonsei-TGIL/CLEMENT.git
2. cd CLEMENT/scripts

## Version update
1.0.0 (Oct 20th, 2022)

## Input format
As now of 1.0.0, CLEMENT only supports standardized TSV input. Examples of input file is shown in "example" directory.
- 1st column:	mutation ID (CHR_POS is recommended)
- 2nd column: label, if possible. If user don't know the label, just set 0
- 3rd column: **Depth1,Alt1,Depth2,Alt2....,Depth_n,Alt_n**    * should be comma-separated, and no space permitted

## Running
1. usage1: python3 CLEMENT.py [OPTIONS]
2. usage2: conda install CLEMENT (on ready)

**options**

	$These options are regarding User's input and output format.$
		--INPUT_TSV	Input data whether TSV or VCF. The tool automatically detects the number of samples. (mandatory)
		--NPVAF_DIR Directory where temporary VAF information TSV file deposits. (mandatory)
		--CLEMENT_DIR Directory where the outputs of CLEMENT be saved. (mandatory)
		--IMAGE_FORMAT Format of image files produced by CLEMENT (default: jpg)  (options: jpg, pdf, svg)

	$These options are regarding downsizing User's input or not$
		--USE_ALL If user want to make use of all the datasets, set True. If user want to downsize the sample, set False and --RANDOM_PICK a integer (default: False)
		--RANDOM_PICK If user want to downsize the sample, set --USE_ALL False and --RANDOM_PICK integer (e.g., --RANDOM-PICK 500)
		--DEPTH_CUTOFF	Threshold of input data according to depth information. About 50% of total coverage is recommended. (default: 100)

	$These options are adjusting E-M algorithm parameter$
		--NUM_CLONE_TRIAL_START Minimum number of expected cluster_hards (initation of K) (default: 2)
		--NUM_CLONE_TRIAL_END Maximum number of expected cluster_hards (termination of K) (default: 7)
		--TRIAL_NO Trial number in each candidate cluster_hard number. DO NOT recommend over 15 (default: 5)
		--KMEANS_CLUSTERNO	Number of initial K-means cluster_hard (default: 8)
		--MIN_CLUSTER_SIZE	The minimum cluster size that is acceptable (default: 15)

	$Miscelleneous$
		--VERBOSE	Print process →  0: no record,  1: simplified record,  2: verbose record (default: 2)


## Example
	-python3 CLEMENT.py --INPUT_TSV "../examples/example1/input.tsv" --NPVAF_DIR "../examples/example1"   --CLEMENT_DIR "../examples/example1"   --NUM_CLONE_TRIAL_START 2 --NUM_CLONE_TRIAL_END 6 --DEPTH_CUTOFF 10 --VERBOSE 2 --TRIAL_NO 3  --RANDOM_SEED 1
	-python3 CLEMENT.py --INPUT_TSV "../examples/example2/input.tsv" --NPVAF_DIR "../examples/example2"   --CLEMENT_DIR "../examples/example2"   --NUM_CLONE_TRIAL_START 2 --NUM_CLONE_TRIAL_END 6 --DEPTH_CUTOFF 10 --VERBOSE 1 --TRIAL_NO 3  --RANDOM_SEED 1

## Contact
	goldpm1@yuhs.ac
