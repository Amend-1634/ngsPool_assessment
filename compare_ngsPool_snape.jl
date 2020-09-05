
#compare the response of two software (ngsPool and Snape) to significance levels employed
using Statistics, StatsPlots


# Estimated MAF plotted with true MAF
cd("C:\\simu_raw\\sig")
include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
name="simu-1-50"
df_true = read_true(name)

name="simu-1-50-0.001"
df_snp=read_ngsPool_snp(name)
snape_p, snape_nallele = read_snape(name)

#xyplot
df_true[:,:true_maf]
df_snp[!,3:end]
snape_p[!,:maf_snape]

# @df mydata scatter(:a, :b)
all_sites=outerjoin(df_true, df_snp, snape_p, on=[:chrom, :position])

all_sites2=all_sites[1:1000,:] #401 the first point?

est2=convert(Matrix, all_sites2[:,[:GSS_ngsPool, :maf_snape]])
    scatter(all_sites2.true_maf, est2, xscale=:log10, title="Estimation of minor allele frequency",
        label=["GSS_ngsPool" "snape"], legend=:topleft, dpi=300)
    plot!(all_sites2.true_maf,all_sites2.true_maf, label="true")



##boxplot comparing the minor allele frequency of snape and ngsPool among significance level (x axis)
cd("C:\\simu_raw\\sig")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
nsamp=50 #no 100 as snape doesn't have data
	# dpt=1
	dpts = [1,5]
	sigs = [0.1, 0.05, 0.01, 0.005, 0.001]
	p = ngsPool_snape_sig_maf(nsamp, dpts, sigs)
	#with plots (suplots layout: row as depth, column as significance level)

#arrays of plots and matrixs
cd("C:\\simu_raw\\sig")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
	nsamp=50 #no 100 as snape doesn't have data
	dpts = [1,5]
	sigs = [0.1, 0.05, 0.01, 0.005, 0.001]
	# dpts = [1]
	# sigs = [0.1]
	xyplots, rmse_plots, mpe_plots, xy_all_p, rmse_all_p, mpe_all_p= all_3plots(nsamp, dpts, sigs)
	#xyplots: true MAF vs estimated MAF
	#rmse: root mean square error
	#mpe: mean percentage error
	#xy_all_p: the MAF matrix containing true MAF and estimtes
