# CLEMENT
- Genomic decomposition and reconstruction of **non-tumor** diploid subclones (2023)
- CLonal decomposition via Expectation-Maximization algorithm established in Non-Tumor setting
- Support multiple diploid sample
- Biallelic variants (Homo, 1/1) can degrade the performance of CLEMENT.

## Overview of CLEMENT workflow and core algorithms
<br/>

![CLEMENT_overview](https://github.com/Yonsei-TGIL/CLEMENT/assets/111939069/e8ff11b3-5fa8-4e2e-b045-47e4da90b01c)
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

### Install from github
1. git clone https://github.com/Yonsei-TGIL/CLEMENT.git   
    cd CLEMENT   
    pip3 install .   

2. pip3 install git+https://github.com/Yonsei-TGIL/CLEMENT.git    

### Install from PyPi
3. pip3 install CLEMENTDNA   

## Version update
1.0.5 (Jan 1st, 2024)

## Input format
As now of 1.0.4, CLEMENT only supports standardized TSV input. Examples of input file is shown in _"example"_ directory.
- 1st column: mutation ID (CHR_POS is recommended)
- 2nd column: label (answer), if possible. If user don't know the label (answer), just set 0
- 3rd column: **Depth1,Alt1,Depth2,Alt2....,Depth_n,Alt_n**    * should be comma-separated, and no space permitted
- 4th column: **BQ1,BQ2....,BQ_n**    * should be comma-separated, and no space permitted. If absent, CLEMENT set default BQ as 20.

## Running
### command line interface
	CLEMENT [OPTIONS]   


### options

	(Mandatory) These options are regarding User's input and output format
		--INPUT_TSV		Input data whether TSV. The tool automatically detects the number of samples
		--CLEMENT_DIR 		Directory where the outputs of CLEMENT be saved

	These options are regarding downsizing User's input or not
		--RANDOM_PICK 		Set this variable to user want to downsize the sample. If user don't want to downsize, set -1. (default : -1).
	
	These options are adjusting E-M algorithm parameter
		--NUM_CLONE_TRIAL_START 	Minimum number of expected cluster_hards (initation of K) 	(default: 3)
		--NUM_CLONE_TRIAL_END 		Maximum number of expected cluster_hards (termination of K)	 (default: 5)
		--TRIAL_NO 			Trial number in each candidate cluster_hard number. DO NOT recommend over 15 (default: 5)
    		--FP_PRIOR FP_PRIOR   		Prior of false positive (FP) being generated. Recommendation : <= 0.1. Default : 0.01
  		--TN_PRIOR TN_PRIOR   		Confidentiality that negative being negative (TN). Recommendation : > 0.99. Default : 0.99
		--KMEANS_CLUSTERNO		Number of initial K-means cluster. Recommendation : 5~8 for one-sample, 8-15 for larger-sample (default: 8)
		--MIN_CLUSTER_SIZE		The minimum cluster size that is acceptable. Recommendation : 1-3% of total variants number 	(default: 9)

	Other options
		--MODE			Selection of clustering method.
					"Hard": hard clustering only,  "Both": both hard and soft (fuzzy) clustering (default: "Both")
		--MAKEONE_STRICT  	1: strict, 2: lenient, 3: most lenient (default : 1)
  		--SCORING		True : comparing with the answer set, False : just visualization (default: False)
		

	Miscelleneous
		--FONT_FAMILY		Font family that displayed in the plots (default : "arial")
  		--VISUALIZATION		Whether produce image in every E-M step (default: True)
    		--IMAGE_FORMAT		Image format that displayed in the plots (default : jpg)
		--VERBOSE		0: no record,  1: simplified record,  2: verbose record (default: 2)


### output

**${CLEMENT_DIR}"/result"**
- **CLEMENT_decision**		_CLEMENT's best recommendation among hard and soft clustering._
- **CLEMENT_hard_1st**  	_CLEMENT's best decomposition by hard clustering._
- **CLEMENT_hard.gapstatistics.txt** 	_Selecting the optimal K in hard clustering based on gap* stastics._
- **CLEMENT_soft_1st** 	_CLEMENT's best decomposition by soft (fuzzy) clustering._
- **membership.txt** 	_Membership assignment of all variants to each clusters._
- **membership_count.txt** 	_Count matrix of the membership assignment to each clusters._
- **mixture.txt** 	_Centroid of each clusters_

## Example
	DIR=[YOUR_DIRECTORY]

	# Example 1
	CLEMENT \
		--INPUT_TSV ${DIR}"/example/2.CellData/MRS_2D/M1-8_M2-4/M1-8_M2-4_input.txt" \
		--CLEMENT_DIR ${DIR}"/example/2.CellData/MRS_2D/M1-8_M2-4"  \
		--NUM_CLONE_TRIAL_START 2 \
		--NUM_CLONE_TRIAL_END 6 \
		--RANDOM_PICK 500
	
		


![example1](https://github.com/Yonsei-TGIL/CLEMENT/assets/56012432/a5a6beb2-e0ac-44ad-8a5a-1b9aa4480010)
![example2](https://github.com/Yonsei-TGIL/CLEMENT/assets/56012432/3ee2c4a3-4627-40a3-80e6-666a981a6c20)
<br/>

## Contact
	goldpm1@yuhs.ac


