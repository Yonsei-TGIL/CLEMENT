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
1. git clone https://github.com/Yonsei-TGIL/CLEMENT.git   
    cd CLEMENT/source   
*or*   
2. pip3 install CLEMENT-DNA  (to be updated) 
*or*   
3. conda install CLEMENT-DNA (to be updated)

## Version update
1.0.1 (June 12th, 2023)

## Input format
As now of 1.0.1, CLEMENT only supports standardized TSV input. Examples of input file is shown in "example" directory.
- 1st column: mutation ID (CHR_POS is recommended)
- 2nd column: label  (answer), if possible. If user don't know the label (answer), just set 0
- 3rd column: **Depth1,Alt1,Depth2,Alt2....,Depth_n,Alt_n**    * should be comma-separated, and no space permitted
- 4th column: **BQ1,BQ2....,BQ_n**    * should be comma-separated, and no space permitted. If absent, CLEMENT set default BQ as 20.

## Running
python3 CLEMENT-DNA.py [OPTIONS]


**options**

	$These options are regarding User's input and output format.$
		--INPUT_TSV	Input data whether TSV. The tool automatically detects the number of samples. (mandatory)
		--CLEMENT_DIR Directory where the outputs of CLEMENT be saved. (mandatory)

	$These options are regarding downsizing User's input or not$
		--RANDOM_PICK Set this variable to user want to downsize the sample. If user don't want to downsize, set -1 (default).
	
	$These options are adjusting E-M algorithm parameter$
		--NUM_CLONE_TRIAL_START Minimum number of expected cluster_hards (initation of K) (default: 3)
		--NUM_CLONE_TRIAL_END Maximum number of expected cluster_hards (termination of K) (default: 5)
		--TRIAL_NO Trial number in each candidate cluster_hard number. DO NOT recommend over 15 (default: 5)
		--KMEANS_CLUSTERNO	Number of initial K-means cluster. Recommendation : 5~8 for one-sample, 8-15 for larger-sample (default: 8)
		--MIN_CLUSTER_SIZE	The minimum cluster size that is acceptable (default: 9)

	$Other options$
		--MODE	Selection of clustering method  "Hard" : hard clustering only,  "Both" : both hard and soft (fuzzy) clustering (default: "Both")
		--MAKEONE_STRICT  1:strict, 2:lenient. (default : 1)
		--TN_CONFIDENTIALITY  Confidentiality that negative being negative (TN). Recommendation : > 0.99. (default : 0.995)

	$Miscelleneous$
		--VERBOSE	Print process →  0: no record,  1: simplified record,  2: verbose record (default: 2)


**output**

	${CLEMENT_DIR}"/result"
		CLEMENT_decision	CLEMENT's best recommendation among hard and soft clustering.
		CLEMENT_hard_1st CLEMENT's best decomposition by hard clustering.
		CLEMENT_hard.gapstatistics.txt Selecting the optimal K in hard clustering based on gap* stastics.
		CLEMENT_soft_1st CLEMENT's best decomposition by soft (fuzzy) clustering.
		membership.txt Membership assignment of all variants to each clusters. 
		membership_count.txt Count matrix of the membership assignment to each clusters.
		mixture.txt Centroid of each clusters

## Example
	python3 CLEMENT-DNA.py --INPUT_TSV "../example/2.CellData/MRS_2D/M1-5_M1-6/M1-5_M1-6_input.txt" --CLEMENT_DIR "../example/2.CellData/MRS_2D/M1-5_M1-6"  --NUM_CLONE_TRIAL_START 2 --NUM_CLONE_TRIAL_END 6 --RANDOM_PICK 500

## Contact
	goldpm1@yuhs.ac
