# ngsPool_assessment

## Introduction
  This project aims to compare the software ngsPool with commonly used software for Pooled-seq data such as Popoolation2, Snape and VarScan

  For a full description and installation of the software ngsPool, visit the project page: https://github.com/mfumagalli/ngsJulia.git

## Installation
- Popoolation2 (https://sourceforge.net/projects/popoolation2/)
- Snape (https://code.google.com/archive/p/snape-pooled/downloads)
- VarScan (https://sourceforge.net/projects/varscan/; mannual http://varscan.sourceforge.net/using-varscan.html#v2.3_mpileup2snp)

## Annotation of scripts

- func_read_standard.jl # read the minor allele frequency (MAF) estimation from the output file of software
- mannual_output.jl # annotation of each column of output file of 4 software
- func_analysis_standard.jl # standardized functions for comparing by plotting

- compare3estimates_ngsPool.jl # compare 3 estimates of ngsPool: maximum likelihood (GSS), maximum likelihood (SFS), expected estimate
- compare_Snape.jl # exploring the setting of Snape (i.e. infuence of informative prior and flat prior MAF)
- compare_ngsPool_snape.jl # compare the SNP calling by two software at different singificance levels
- compare_snape_varscan.jl

- compare_power_qqVector.jl  #analysises of which the true MAFs are fixed (specially focus on low MAFs)
- compare_power_sample.jl #analysises of which the true MAFs are nromally sampled (the sampling is simulating the SFS of real data)

- real_data_pipeline # pipeline in bash script to get the estimation of MAF from real data (FASTQ file) 
         







