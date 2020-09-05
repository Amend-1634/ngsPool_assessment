##ngsPool output######################################
#the first row:
 # chrom   position        reference       nonreference    major   minor   lrtSNP  lrtBia  lrtTria maf     freqMax freqE
#corresponding annotation:
   #chromosome, position,reference allele, nonreference allele, major allele, minor allele,
   #likelihood ration test statistic (lrt) of SNP, lrt of biallelic site, lrt of triallelic,
   #minor allele frequency (from Golden Section Search),
   #Maximum likelihood estimation (from SFS)
   #Expected estimation (from SFS)

##Popoolation2######################################
# (source: https://sourceforge.net/p/popoolation2/wiki/Tutorial/)
# Sample of a synchronized file:
#
# 2R  2302    N   0:7:0:0:0:0 0:7:0:0:0:0
# 2R  2303    N   0:8:0:0:0:0 0:8:0:0:0:0
# 2R  2304    N   0:0:9:0:0:0 0:0:9:0:0:0
# 2R  2305    N   1:0:9:0:0:0 0:0:9:1:0:0
# col1: reference contig
# col2: position within the refernce contig
# col3: reference character
# col4: allele frequencies of population number 1
# col5: allele frequencies of population number 2
# coln: allele frequencies of population number n
# The allele frequencies are in the format A:T:C:G:N:del, i.e: count of bases 'A',
# count of bases 'T',... and deletion count in the end (character '*' in the mpileup)

##Snape output######################################
#10 fields:
    #1. chr: genomic coordinates
    #2. position along the chronosome
    #3. ref nucleotides
    #4. num of minor (alternative) nucleotide
    #5. mean quality of the reference nucleotide
    #6. mean quality of the alternative nucleotide
    #7. first and second most frequent nucleotides in the pileup
    #8. 1-p(0), p(f) as the probability distribution function for the maf (f)
    #9. p(1)
    #10. mean value of f

    # what's num of minor (alternative) nucleotide for?
    # 5.6. are around 60 and 70, seems unrealistic(0.0001% error rate) compared to
        #(Kofler et al., 2011) Popoolation, should be around 20 to 30
        #what's the unit?
    # 1-p(0)-p(1) is the probability of SNP? Any threshold?
    # why 1-p(0) can be greater than 1? (no matter accumulative or density)
    # can 10. be compared with the expected mean in ngsPool?



##Varscan#############################################
# OUTPUT
# 	Tab-delimited SNP calls with the following columns:
# 	1.Chrom		chromosome name
# 	2.Position	position (1-based)
# 	3.Ref			reference allele at this position
# 	4.Var			variant allele observed
# 	5.PoolCall	Cross-sample call using all data (Cons:Cov:Reads1:Reads2:Freq:P-value)
# 			1)Cons - consensus genotype in IUPAC format
# 			2)Cov - total depth of coverage
# 			3)Reads1 - number of reads supporting reference
# 			4)Reads2 - number of reads supporting variant
# 			5)Freq - the variant allele frequency by read count
##5.5)lower false discovery rate in SNP calling
# 			6)P-value - FET p-value of observed reads vs expected non-variant
# 	6.StrandFilt	Information to look for strand bias using all reads (R1+:R1-:R2+:R2-:pval)
# 			R1+ = reference supporting reads on forward strand
# 			R1- = reference supporting reads on reverse strand
# 			R2+ = variant supporting reads on forward strand
# 			R2- = variant supporting reads on reverse strand
# 			pval = FET p-value for strand distribution, R1 versus R2
# 	7.SamplesRef	Number of samples called reference (wildtype)
# 	8.SamplesHet	Number of samples called heterozygous-variant
# 	9.SamplesHom	Number of samples called homozygous-variant
# 	10.SamplesNC	Number of samples not covered / not called
# 	11.SampleCalls	The calls for each sample in the mpileup, space-delimited
#     			Each sample has six values separated by colons:
# 			Cons - consensus genotype in IUPAC format
# 			Cov - total depth of coverage
# 			Reads1 - number of reads supporting reference
# 			Reads2 - number of reads supporting variant
# 			Freq - the variant allele frequency by read count
##11.5)collect as
# 			P-value - FET p-value of observed reads vs expected non-variant
                    #fisher's exact test?
