##Exploring software Snape and compare ngsPool and Snape

    #also tested why snape are all empty when sample size=10

#-nchr \$nChroms => -nchr $nsamp
#maxcov back to be 100

#compare the significance threshold of ngsPool and Snape:
using Distributions #using HypothesisTests, Distributions
cd("C:\\simu_raw\\pipeline")
    f=open("sig_script", "w")
        pre_cmd=
        "cd /mnt/c/simu_raw/sig #simulation data
    # mkdir sig; cd sig
    other_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\"
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
        "
        write(f, pre_cmd, "\n")

        for dpt in [1,5], nsamp in [50, 100], sig in [0.1, 0.05, 0.01, 0.005, 0.001]
            lrtSnp=cquantile(Chisq(1), sig)
            pp_snape=1-sig

            cmd=
            "lrtSnp=$lrtSnp
    name=\"simu-$dpt-$nsamp\"
    nChroms=$(2*nsamp)
        # time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out \$name.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$name.mpileup.gz
        # julia \$ngsPool_dir\"ngsPool.jl\" --fin \$name.mpileup.gz --fout \$name-$sig.snp.gz --nChroms \$nChroms --lrtSnp $lrtSnp
        # gzip -d \$name.mpileup.gz

        # echo \"npstat time\"
        # maxcov=\$nChroms*$dpt #set as default 100
        # -n as haploid sample size (what about diploid?)
    #time \$other_dir\"npstat\" -n \$nChroms -l 1000 -mincov 4 -maxcov 100 -minqual 20 -nolowfreq 2 \$name.mpileup #default output samename.stats

    theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$name.mpileup.stats) #nucleotide diversity
    D=\$(awk '{sum+=\$13} END {print sum/NR}' \$name.mpileup.stats) #Prior genetic difference between reference genome and population
    echo \"Chromosomes \$name, theta is \$theta, D is \$D\"

    echo \"snape time\"
    time \$other_dir\"snape-pooled\" -nchr \$nChroms -theta \$theta -D \$D -fold unfolded -priortype informative < \$name.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='2'' > \$name-$sig.snape
    "
            println(cmd)
            println("")
            write(f, cmd, "\n")
        end
        close(f)

#only 5-100 combination has empty snape output


#######################################test_npstat
using Distributions #using HypothesisTests, Distributions
cd("C:\\simu_raw\\pipeline")
    f=open("test_npstat", "w")
    pre_cmd=
    "cd /mnt/c/simu_raw/sig #simulation data
    # mkdir sig; cd sig
    other_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\"
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    "
    write(f, pre_cmd, "\n")

    for dpt in [1,5], nsamp in [50, 100], sig in [0.1, 0.05, 0.01, 0.005, 0.001]
        lrtSnp=cquantile(Chisq(1), sig)
        pp_snape=1-sig

        cmd=
        "
    lrtSnp=$lrtSnp
    name=\"simu-$dpt-$nsamp\"
    nChroms=$(2*nsamp)
    maxcov=\$nChroms*$dpt

    # time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out test-$dpt-$nsamp.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$name.mpileup.gz
    # julia \$ngsPool_dir\"ngsPool.jl\" --fin \$name.mpileup.gz --fout \$name-$sig.snp.gz --nChroms \$nChroms --lrtSnp $lrtSnp
    # gzip -d \$name.mpileup.gz

    echo \"npstat time\"
    time \$other_dir\"npstat\" -n $(2*nsamp) -l 1000 -mincov 4 -maxcov \$maxcov -minqual 20 -nolowfreq 2 \$name.mpileup #default output samename.stats

    theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$name.mpileup.stats) #nucleotide diversity
    D=\$(awk '{sum+=\$13} END {print sum/NR}' \$name.mpileup.stats) #Prior genetic difference between reference genome and population
    echo \"Chromosomes \$name, theta is \$theta, D is \$D\"

    echo \"snape time\"
    time \$other_dir\"snape-pooled\" -nchr \$nChroms -theta \$theta -D \$D -fold unfolded -priortype informative < \$name.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='2'' > \$name-$sig.snape
    "
        println(cmd)
        println("")
        write(f, cmd, "\n")
    end
    close(f)

#when it finish, check
cd /mnt/c/simu_raw/sig
wc *snape*
#still zero #not due to max coverage

#pilot: as replicate number=1


##
cd("C:\\simu_raw\\sig")
include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_plot.jl")
for dpt in [1,5], nsamp in [50, 100], sig in [0.1, 0.05, 0.01, 0.005, 0.001]
    pre_name="simu-$dpt-$nsamp"
    $name-$sig.snp.gz
    name="$pre_name-$sig.snape"
    snape_p, snape_nallele = read_snape(name)


#############################test_snape
using Distributions #using HypothesisTests, Distributions
cd("C:\\simu_raw\\pipeline")
    f=open("test_snape", "w")
    pre_cmd=
    "cd /mnt/c/simu_raw #simulation data
    mkdir sig_snape; cd sig_snape
    other_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\"
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    "
    write(f, pre_cmd, "\n")

    # for dpt in [1,5], nsamp in [50, 100], sig in [0.1, 0.05, 0.01, 0.005, 0.001]
    for nsamp in [50,70,90, 110]
        sig=0.05
        dpt=1
        lrtSnp=cquantile(Chisq(1), sig)
        pp_snape=1-sig

        cmd=
        "
    lrtSnp=$lrtSnp
    name=\"simu-$dpt-$nsamp\"
    nChroms=$(2*nsamp) #total chromosome pooled
    maxcov=\$nChroms*$dpt

    time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out test-$dpt-$nsamp.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$name.mpileup.gz
    # julia \$ngsPool_dir\"ngsPool.jl\" --fin \$name.mpileup.gz --fout \$name-$sig.snp.gz --nChroms \$nChroms --lrtSnp $lrtSnp
    gzip -d \$name.mpileup.gz

    echo \"npstat time\"
    time \$other_dir\"npstat\" -n $(2*nsamp) -l 1000 -mincov 4 -maxcov \$maxcov -minqual 20 -nolowfreq 2 \$name.mpileup #default output samename.stats

    theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$name.mpileup.stats) #nucleotide diversity
    D=\$(awk '{sum+=\$13} END {print sum/NR}' \$name.mpileup.stats) #Prior genetic difference between reference genome and population
    echo \"Chromosomes \$name, theta is \$theta, D is \$D\"

    echo \"snape time\"
    time \$other_dir\"snape-pooled\" -nchr \$nChroms -theta \$theta -D \$D -fold unfolded -priortype informative < \$name.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='2'' > \$name-$sig.snape
    "
        println(cmd)
        println("")
        write(f, cmd, "\n")
    end
    close(f)

    #all except 50 are empty #[50,70,90, 110]

cd("C:\\simu_raw\\pipeline")
    f=open("test_snape2", "w")
    pre_cmd=
    "cd /mnt/c/simu_raw #simulation data
    #mkdir sig_snape;
    cd sig_snape
    other_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\"
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    "
    write(f, pre_cmd, "\n")

    # for dpt in [1,5], nsamp in [50, 100], sig in [0.1, 0.05, 0.01, 0.005, 0.001]
    for nsamp in [30, 40, 55, 60]
        sig=0.05
        dpt=1
        lrtSnp=cquantile(Chisq(1), sig)
        pp_snape=1-sig

        cmd=
        "lrtSnp=$lrtSnp
        name=\"simu-$dpt-$nsamp\"
        nChroms=$(2*nsamp)
        time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out \$name.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$name.mpileup.gz
        # julia \$ngsPool_dir\"ngsPool.jl\" --fin \$name.mpileup.gz --fout \$name-$sig.snp.gz --nChroms \$nChroms --lrtSnp $lrtSnp
        gzip -d \$name.mpileup.gz

        # echo \"npstat time\"
        # maxcov=\$nChroms*$dpt #set as default 100
        # -n as haploid sample size (what about diploid?)
        time \$other_dir\"npstat\" -n $(2*nsamp) -l 1000 -mincov 4 -maxcov 100 -minqual 20 -nolowfreq 2 \$name.mpileup #default output samename.stats

        theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$name.mpileup.stats) #nucleotide diversity
        D=\$(awk '{sum+=\$13} END {print sum/NR}' \$name.mpileup.stats) #Prior genetic difference between reference genome and population
        echo \"Chromosomes \$name, theta is \$theta, D is \$D\"

        echo \"snape time\"
        time \$other_dir\"snape-pooled\" -nchr $nsamp -theta \$theta -D \$D -fold unfolded -priortype informative < \$name.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='2'' > \$name-$sig.snape
        "
        println(cmd)
        println("")
        write(f, cmd, "\n")
    end
    close(f) #[30, 40, 55, 60]

#after change  -nchr \$nChroms => -nchr $nsamp; before beyond 1-50 => no data
 #no it's beyond 1-60
# 60   660  3208 simu-1-60-0.05.snape
# 0     0     0 simu-1-70-0.05.snape

##################does informative and flat prior MAF make a difference in the estimation of Snape
#if so, what significance level of informative/flat Snape are comparable with what sig of ngsPool

using Distributions #using HypothesisTests, Distributions
cd("C:\\simu_raw\\pipeline")
    f=open("test_snape_flat", "w")
    pre_cmd=
    "cd /mnt/c/simu_raw #simulation data
    # mkdir sig_snape;
    cd sig_snape
    other_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\"
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    "
    write(f, pre_cmd, "\n")

    # for dpt in [1,5], nsamp in [50, 100], sig in [0.1, 0.05, 0.01, 0.005, 0.001]
    for nsamp in [20,30,40,50,55,60]

        sig=0.05
        dpt=1
        lrtSnp=cquantile(Chisq(1), sig)
        pp_snape=1-sig
        name="simu-$dpt-$nsamp"

        cmd=
        "
    lrtSnp=$lrtSnp
    nChroms=$(2*nsamp) #total chromosome pooled
    maxcov=\$nChroms*$dpt

    time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out $name.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > $name.mpileup.gz
    # julia \$ngsPool_dir\"ngsPool.jl\" --fin $name.mpileup.gz --fout $name-$sig.snp.gz --nChroms \$nChroms --lrtSnp $lrtSnp
    gzip -d $name.mpileup.gz

    echo \"npstat time\"
    time \$other_dir\"npstat\" -n $(2*nsamp) -l 1000 -mincov 4 -maxcov \$maxcov -minqual 20 -nolowfreq 2 $name.mpileup #default output samename.stats

    theta=\$(awk '{sum+=\$6} END {print sum/NR}' $name.mpileup.stats) #nucleotide diversity
    D=\$(awk '{sum+=\$13} END {print sum/NR}' $name.mpileup.stats) #Prior genetic difference between reference genome and population
    echo \"Chromosomes $name, theta is \$theta, D is \$D\"

    echo \"snape time\"
    #informative prior
    time \$other_dir\"snape-pooled\" -nchr \$nChroms -theta \$theta -D \$D -fold unfolded -priortype informative < $name.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='2'' > $name-$sig.snape
    #flat prior
    time \$other_dir\"snape-pooled\" -nchr \$nChroms -theta \$theta -D \$D -fold unfolded -priortype flat < $name.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='2'' > $name-$sig-flat.snape
    "
        println(cmd)
        println("")
        write(f, cmd, "\n")
    end
    close(f)

#rows of flat is almost half of the informative ones (the upper boundary of sample size is still 55)
cd("C:\\simu_raw\\sig_snape")
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")

name="simu-1-30"
df_true = read_true(name)

name="simu-1-30-0.05"
df_snape = read_snape(name)[1]
df_snape_flat = read_snape("$name-flat")[1]

# df_ngs=read_ngsPool_snp(name)
using DataFrames
# df_all=outerjoin(df_true,df_snape,df_snape_flat, on = (:chrom,:position))
# df_all=outerjoin(df_true,df_snape,df_snape_flat, on = [:chrom,:position])
df_all=outerjoin(df_true,df_snape,df_snape_flat, on = [:chrom,:position],makeunique=true)
# names(df_all)
filt_df_all=filter([:maf_snape, :maf_snape_1] => (maf_snape, maf_snape_1) -> !ismissing(maf_snape) || !ismissing(maf_snape_1), df_all)
gr()
Plots.plot(Matrix(filt_df_all[:, [:maf_snape, :maf_snape_1]]), label=["informative" "flat"],
    title="Snape-depth1x-samp30-significancance0.05",dpi=300,grid=false,
    xlabel="sites", ylabel="Minor allele frequency")
using Plots
Plots.savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\snape-info-flat.png")

#flat calls even less sites
    #  13   143   686 simu-1-20-0.05-flat.snape
    #  23   253  1191 simu-1-20-0.05.snape
    #  11   121   592 simu-1-30-0.05-flat.snape
    #  17   187   888 simu-1-30-0.05.snape
    #  18   198   963 simu-1-40-0.05-flat.snape
    #  31   341  1614 simu-1-40-0.05.snape
    #  34   374  1805 simu-1-50-0.05-flat.snape
    #  56   616  2952 simu-1-50-0.05.snape

################################################################################
# analysis -- COMPARE NGSPOOL AND SNAPE
################################################################################


include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
using Statistics, StatsPlots
# "$name.txt" #population allele frequency simulated
# "$name.snp.gz" #estimation of ngsPool with SNP calling
# "$name.sync" #Popoolation2
# "$name.snape"
# "$name.varscan"

#compare significance threshold



cd("C:\\simu_raw\\sig")
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

# println(names(all_sites))
# plot(df_true[!,:true_maf])

# @df all_sites scatter(:true_maf, [:GSS_ngsPool,:freqMax_ngsPool,:freqE_ngsPool, :maf_snape])
# est=convert(Matrix, all_sites2[:,[:GSS_ngsPool,:freqMax_ngsPool,:freqE_ngsPool, :maf_snape]])
    est2=convert(Matrix, all_sites2[:,[:GSS_ngsPool, :maf_snape]])
    scatter(all_sites2.true_maf, est2, xscale=:log10, title="Estimation of minor allele frequency",
        label=["GSS_ngsPool" "snape"], legend=:topleft, dpi=300)
        plot!(all_sites2.true_maf,all_sites2.true_maf, label="true")
# savefig("H\\julia\\ngsJulia\\ngsPool\\geneticdrift_sd.jpg")####hang on
#1e-3:1e-1:0.5



include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_plot.jl")
# all_sites

# using StatsPlots
boxplot(rmse_mtx,["1", "2", "3", "4"]) #rmse distribution

plot!(label=["1", "2", "3", "4"])
boxplot(rmse_mtx, ["1", "2", "3", "4"]) #rmse distribution
boxplot(rmse_mtx) #rmse distribution
scatter!(rmse_ave_vec, color=:black) #mean

boxplot_missing(mpe_mtx) #all are very close to 0, snape seem to be with less outlier


 #the snape calls are too little? only call at high confidence?

# scatter!(mpe_ave_vec)
mpe_ave_vec #all missing


#RMSE

######################
##boxplot comparing the minor allele frequency of snape and ngsPool among significance level (x axis)
cd("C:\\simu_raw\\sig")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
    nsamp=50 #no 100 as snape doesn't have data
	# dpt=1
	dpts = [1,5]
	sigs = [0.1, 0.05, 0.01, 0.005, 0.001]
	p = ngsPool_snape_sig_maf(nsamp, dpts, sigs)
#no difference between significance...
#problem in plotting but have no time to check (2020/8/14)

#collect detective minor allele frequency

# for dpt in [1,5], sig in [0.1, 0.05, 0.01, 0.005, 0.001]

cd("C:\\simu_raw\\sig")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
	nsamp=50 #no 100 as snape doesn't have data
	dpts = [1,5]
	sigs = [0.1, 0.05, 0.01, 0.005, 0.001]
	# dpts = [1]
	# sigs = [0.1]
	xyplots, rmse_plots, mpe_plots, xy_all_p, rmse_all_p, mpe_all_p= all_3plots(nsamp, dpts, sigs)

    #after setting all missing point (sites not called), the advantage in RMSE of ngsPool is flattened.

nsamp=50 #no 100 as snape doesn't have data
	dpts = [1,5]
	sigs = [0.1, 0.05, 0.01, 0.005, 0.001]
	xyplots = sig_all_xyplots(nsamp, dpts, sigs)
xyplots


gr()
rmse_all_p
# rmse_all_p'
# plot(rmse_all_p..., layout=(5,2))
# plot(rmse_all_p..., layout=(2,5))
# plot(rmse_all_p...)
mpe_all_p

# using PyPlot
# plot(xyplots..., layout=(2,5),titlefontsize=11, size=(1400,600), margin=1mm,
# 	xlabel="", ylabel="")
# 	# annotate!([(700,300,text("Significance level", 48))])
# 	annotate!([(.6, 1.5, text(L"L = M\frac{KP}{1+KP}",48))])
#
#
# annotate!(2050, 1/n_s, text("1/n_s",f, :left))

# plot!(xlabel=false, ylabel=false)
# plot!(xlabel="significance level")
xy_all_p
rmse_all_p
mpe_all_p
