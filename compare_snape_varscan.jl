#A trial: compare Snape and VarScan estimation with real data (not included in report)

#Snape outpu
cd("C:\\simu_raw\\sig_snape")
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")

cd("C:\\Fraca_raw\\out")
using CSV, DataFrames
name="PopB_pooled_L4"

varscan_p=read_varscan(name) #maf and Fisher's exact test p values?
varscan_nallele= read_nallele_varscan2(name) #mumber of reference or alternative allele

snape_p, snape_nallele = read_snape(name)

snape_p #275
varscan_p #108
snape_varscan_p = innerjoin(snape_p, varscan_p, on = [:chrom, :position])
#74

snape_nallele #275
varscan_nallele #108
snape_varscan_nallele = innerjoin(snape_nallele, varscan_nallele, on = [:chrom, :position])
#74

#plot estimation by Software Snape (x axis) with Varscan estimation (y axis)
scatter(snape_varscan_p[:,3],snape_varscan_p[:,4], xlabel="snape_maf", ylabel="varscan_maf", label="cross-sample", legend=:bottomright) #differ
    scatter!(snape_varscan_p[:,3],snape_varscan_p[:,6], label="each-sample") #differ
    Plots.abline!(1, 0, line=:dash, label = "")#abline(slope,intercept)
#each-sample test are closer to snape result while cros sample are sparser

scatter(snape_varscan_nallele[:,4],snape_varscan_nallele[:,6],
xlabel="alternative-snape", ylabel="alt-varscan",legend=:false) #referece
#why varscan tend to have a much higher count of alterntive allele
#the quality and read number of snape is absurd...why?

# scatter(snape_varscan_nallele[:,3],snape_varscan_nallele[:,5]) #alternative allele
#snape ref actually categorical(AGCT)

#has two series: each-sample test estimation of VarScan and multiple-sample test of VarScan as y
###function_plot
function plot_snape_varscan(name)
    varscan_p=read_varscan(name)
    varscan_nallele= read_nallele_varscan2(name) #mumber of reference or alternative allele

    snape_p, snape_nallele = read_snape(name)

    print("snape  ", size(snape_p,1))
    print("varscan  ", size(varscan_p,1))
    snape_varscan_p = innerjoin(snape_p, varscan_p, on = [:Chrom, :Position])
    print("inner_join  ", size(snape_varscan_p,1))

    snape_varscan_nallele = innerjoin(snape_nallele, varscan_nallele, on = [:Chrom, :Position])

    p_maf = scatter(snape_varscan_p[:,3],snape_varscan_p[:,4], xlabel="snape_maf", ylabel="varscan_maf", label="cross-sample", legend=:bottomright) #differ
    scatter!(snape_varscan_p[:,3],snape_varscan_p[:,6], label="each-sample") #differ
    Plots.abline!(1, 0, line=:dash, label = "")#abline(slope,intercept)
    display(p_maf)

    p_nallele = scatter(snape_varscan_nallele[:,4],snape_varscan_nallele[:,6],
    xlabel="alternative-snape", ylabel="alt-varscan",legend=:false) #referece
    display(p_nallele)

    return p_maf, p_nallele
end

cd("C:\\Fraca_raw\\out")
using CSV, DataFrames
name="PopB_pooled_L4"
p_maf, p_nallele = plot_snape_varscan(name)


##whether 4 lanes of real data are the similar coverage for the same sites(8scaffold respectively for the 8 chromosomes)?
 #if so I can compare the estimation of maf with elevated depth
 #4 lanes seperately, random combination of 1 2 3 4
