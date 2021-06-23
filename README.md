# ngsPool_assessment

# Introduction
  This project aims to compare the software ngsPool with commonly used software for Pooled-seq data such as Popoolation2, Snape and VarScan

  For a full description and installation of the software ngsPool, visit the project page: https://github.com/mfumagalli/ngsJulia.git

# Installation
- Popoolation2 (https://sourceforge.net/projects/popoolation2/)
- Snape (https://code.google.com/archive/p/snape-pooled/downloads)
- VarScan (https://sourceforge.net/projects/varscan/; mannual http://varscan.sourceforge.net/using-varscan.html#v2.3_mpileup2snp)

# Annotation of scripts
## 1. Scripts and functions
### 1.1 scripts for simulating mpileup
- simulMpileup.R 		      # original script quoted from the repository: SamueleSoraggi/HMMploidy
- simulMpileup_qqVector.R	# modified to simulate mpileup based on fixed Site frequency spectrum
- simulMpileup_qq.R		    # modified to simulate mpileup based on fixed Minor allele frequency (MAF)
-----------------------
### 1.2 functions for Downstream analysis
- mannual_output.jl # annotation of each column of output file of 4 programs (mannual for func_read_standard.jl)
- func_read_standard.jl # read the minor allele frequency (MAF) estimation from the output file of programs
- func_analysis_standard.jl # standardized functions for comparing by plotting
-----------------------
## 2. Script for estimating real data 
- real_data_pipeline # pipeline in bash script to get the estimation of MAF from real data (FASTQ file) 

(Ref: Fracassetti M, Griffin PC, Willi Y. 2015. Validation of Pooled Whole-Genome Re-Sequencing 
in Arabidopsis lyrata. PLoS One. 10:e0140462. doi: 10.1371/journal.pone.0140462.)

-----------------------
## 3. Comparing programs (linux command script writing + analysis in one file) using real and simulated data
### 3.1 Test with different simulated data
- compare_power_qqVector.jl  # calculated power and FPR of SNP calling using simulated data with fixed expected Minor allele frequency (using simulMpileup_qq.R)
- compare_power_sample.jl #calculated power and FPR of SNP calling using simulated data with sampled SFS (using simulMpileup.R: expected minor allele frequency is sampled from sampled Site Frequency Spectrum)

-----------------------
### 3.2 Exploring the setting of program
- compare3estimates_ngsPool.jl # compare 3 estimates of ngsPool: maximum likelihood (GSS), maximum likelihood (SFS), expected estimate
- compare_Snape.jl # exploring the setting of Snape (i.e. infuence of informative prior and flat prior MAF)
- compare_ngsPool_snape.jl # compare the SNP calling by two software at different singificance levels
- compare_snape_varscan.jl









