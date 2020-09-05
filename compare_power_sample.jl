##Script for power and FPR

##using randomly sampled true MAFs (simulMpileup.R)

#settings annotation:
    #sigs: significance levels
    #dpts: depths
    #nsamps: sample sizes

##############################
#Interest: very low minor allele frequency
#settings: sigs=collect(0.001:0.001:0.1); dpts=[3,5]; nsamps=[30, 50]

# 1) Generating linux script
 #significance from 0.001 to 0.1 #but found that this is far from composing the ROC curve
 #script "power_FPR" # directory "power"
cd("C:\\simu_raw\\pipeline") #qq-0
    using Statistics, Distributions
    f=open("power_FPR", "w") #no preset true MAF
        pre_cmd=
    "
    cd /mnt/c/simu_raw/ #fixed population allele frequency
    mkdir power; cd power
    Popoo_dir=\"/mnt/H/julia/ngsJulia/ngsPool/popoolation2_1201/\"
    Fraca_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\" #snape and varscan not written but used by Fraca
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    echo \"\" > time_collection #open a file to collect time information
    "
    write(f, pre_cmd, "\n")
    for dpt in [1,3,5], nsamp in [30, 50]
        # sig in [0.1, 0.05, 0.01, 0.005, 0.001] #more points needed for ROC curve
    # for dpt in [1], nsamp in [30], qq in [1e-4] #test
        # sig=0.05

        pre_cmd="
        mkdir power_rnd; cd power_rnd
        out=\"simu-$dpt-$nsamp\"
        echo \"\$out mpileup##################################\" >> time_collection
        (time Rscript \$ngsPool_dir\"simulMpileup.R\" --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection
        "
        write(f, pre_cmd, "\n")

        for sig in 0.001:0.001:0.1
            lrtSnp=cquantile(Chisq(1), sig)
            pp_snape=1-sig

            cmd ="
            out2=\"simu-$dpt-$nsamp-$sig\" #*: here the last element is sig instead of qq(preset true MAF)
            gzip -d \$out.mpileup.gz
            #1.Popoolation
            echo \"Popoolation java\" >> time_collection
            (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out2.sync) 2>> time_collection

            #parameters for snape or varscan
            #chr_pool=8 #Arabidopsis lyrata has 8 chromosomes
            n_chr=$(2*nsamp) #total chromosomes pooled = diploid * num of individuals
            min=4 #read depth should be at least 4
            max=100 #default maximum read depth
            min_qual=20 #base quality
            l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
            min_all=2 #minimum allele count for Snape, VarScan and NPStat
            # pp_snape=0.9 #posterior probability threshold for Snape

            #2. Snape
            #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
            #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
            # -n as haploid sample size (what about diploid?)

            echo \"npstat time\" >> time_collection
            (time \$Fraca_dir\"npstat\" -n \$n_chr -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out2.stats) 2>> time_collection
            #will output \$out\".stats\" in default

            #NR: number of records variable (awk NR gives the total number of records being processed or line number)
            theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out2\".stats\") #h.mpileup.stats\" as input file #nucleotide diversity
            D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out2\".stats\")
             #Prior genetic difference between reference genome and population
            #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print

            echo \"snape time\" >> time_collection
            (time \$Fraca_dir\"snape-pooled\" -nchr \$n_chr -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='\$min_all'' > \$out2.snape) 2>> time_collection

            #3. VarScan
            echo \"VarScan time\" >> time_collection
            (time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value $pp_snape > \$out2.varscan) 2>> time_collection
            #--p-value as threshold to call variant

            #4. ngsPool
            gzip \$out.mpileup
            echo \"julia time\" >> time_collection
            (time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out2.gz --nChroms \$n_chr --lrtSnp $lrtSnp) 2>> time_collection
            "
            # println(cmd)
            # println("")
            write(f, cmd, "\n")
        end

        #################################set qq as 0
        pre_cmd="
        cd ../
        mkdir power_0; cd power_0
        out=\"simu-$dpt-$nsamp\"
        echo \"\$out mpileup##################################\" >> time_collection
        (time Rscript \$ngsPool_dir\"simulMpileup_qq.R\" --qq 0 --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection
        "
        write(f, pre_cmd, "\n")

        for sig in 0.001:0.001:0.1
            lrtSnp=cquantile(Chisq(1), sig)
            pp_snape=1-sig

            cmd ="
            out2=\"simu-$dpt-$nsamp-$sig\" #*: here the last element is sig instead of qq(preset true MAF)
            gzip -d \$out.mpileup.gz
            #1.Popoolation
            echo \"Popoolation java\" >> time_collection
            (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out2.sync) 2>> time_collection

            #parameters for snape or varscan
            #chr_pool=8 #Arabidopsis lyrata has 8 chromosomes
            n_chr=$(2*nsamp) #total chromosomes pooled = diploid * num of individuals
            min=4 #read depth should be at least 4
            max=100 #default maximum read depth
            min_qual=20 #base quality
            l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
            min_all=2 #minimum allele count for Snape, VarScan and NPStat
            # pp_snape=0.9 #posterior probability threshold for Snape

            #2. Snape
            #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
            #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
            # -n as haploid sample size (what about diploid?)

            echo \"npstat time\" >> time_collection
            (time \$Fraca_dir\"npstat\" -n \$n_chr -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out2.stats) 2>> time_collection
            #will output \$out\".stats\" in default

            #NR: number of records variable (awk NR gives the total number of records being processed or line number)
            theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out2\".stats\") #h.mpileup.stats\" as input file #nucleotide diversity
            D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out2\".stats\")
             #Prior genetic difference between reference genome and population
            #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print

            echo \"snape time\" >> time_collection
            (time \$Fraca_dir\"snape-pooled\" -nchr \$n_chr -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='\$min_all'' > \$out2.snape) 2>> time_collection

            #3. VarScan
            echo \"VarScan time\" >> time_collection
            (time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value $sig > \$out2.varscan) 2>> time_collection
                                                                                                                        #pp_snape amended to be sig after generating data(might be the reason why no much difference made among sigs in varscan)
            #--p-value as threshold to call variant

            #4. ngsPool
            gzip \$out.mpileup
            echo \"julia time\" >> time_collection
            (time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out2.gz --nChroms \$n_chr --lrtSnp $lrtSnp) 2>> time_collection
            "
            # println(cmd)
            # println("")
            write(f, cmd, "\n")
        end

    end
    close(f)

# 2) data analysis--ROC curve
sigs=collect(0.001:0.001:0.1); dpts=[3,5]; nsamps=[30, 50] #settings

include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
    cd("C:\\simu_raw\\power\\power_rnd") #power
power = metric_2mtx(dpts, nsamps, sigs)
power[1,1]

cd("C:\\simu_raw\\power\\power_0") #false positive rate
FPR = metric_2mtx(dpts, nsamps, sigs)
FPR[1,1]
#all element along the column has the same value

#the matrix named as "FPR" and "power":
    #outer matrix: the number of depth as row number, the number of sample size as column number
    #each element of the outer matrix => inner matrix
    #inner matrix: row as significance level, column as corresponding software

Plots.scatter(FPR[1,1], power[1,1], size=(600,500), ylim=(0,1), xlim=(0,1))
    Plots.abline!(1,0)

Plots.scatter(FPR[2,1], power[2,1], size=(600,500), ylim=(0,1), xlim=(0,1))
    Plots.abline!(1,0)

#the array collecting all ROC curve plots under difference combination of sample size(nsamp) and depth(dpt)
n_dpts=length(dpts)
n_nsamps=length(nsamps)
for d in 1:n_dpts, n in 1:n_nsamps
    p=plot(power[d,n], FDR[d,n])
end



##################################
####with wider spectrum of significance level#################################
#significance from 0.1 to 0.99
#settings: depth: 1, sample size: [30, 50]
# significance level: 0.01:0.01:0.99


##bash script to simulate and estimate the MAF of called sites #using simulMpileup.R
    #script "power_FPR_wider", directory "power_wider"
cd("C:\\simu_raw\\pipeline") #qq-0
    using Statistics, Distributions
    f=open("power_FPR_wider", "w") #no preset true MAF
        pre_cmd=
    "
    cd /mnt/c/simu_raw/ #fixed population allele frequency
    mkdir power_wider; cd power_wider
    Popoo_dir=\"/mnt/H/julia/ngsJulia/ngsPool/popoolation2_1201/\"
    Fraca_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\" #snape and varscan not written but used by Fraca
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    echo \"\" > time_collection #open a file to collect time information
    "
    write(f, pre_cmd, "\n")
    # for dpt in [3,5], nsamp in [30, 50]
    for dpt in [1], nsamp in [30, 50]
        # sig in [0.1, 0.05, 0.01, 0.005, 0.001] #more points needed for ROC curve
    # for dpt in [1], nsamp in [30], qq in [1e-4] #test
        # sig=0.05

        pre_cmd="
        mkdir power_rnd; cd power_rnd
        out=\"simu-$dpt-$nsamp\"
        echo \"\$out mpileup##################################\" >> time_collection
        (time Rscript \$ngsPool_dir\"simulMpileup.R\" --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection

        gzip -d \$out.mpileup.gz
        #1.Popoolation
        echo \"Popoolation java\" >> time_collection
        (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out.sync) 2>> time_collection
        "
        write(f, pre_cmd, "\n")

        for sig in 0.01:0.01:0.99
            lrtSnp=cquantile(Chisq(1), sig)
            pp_snape=1-sig

            cmd ="
            out2=\"simu-$dpt-$nsamp-$sig\" #*: here the last element is sig instead of qq(preset true MAF)

            #parameters for snape or varscan
            #chr_pool=8 #Arabidopsis lyrata has 8 chromosomes
            n_chr=$(2*nsamp) #total chromosomes pooled = diploid * num of individuals
            min=4 #read depth should be at least 4
            max=100 #default maximum read depth
            min_qual=20 #base quality
            l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
            min_all=2 #minimum allele count for Snape, VarScan and NPStat
            # pp_snape=0.9 #posterior probability threshold for Snape

            #2. Snape
            #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
            #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
            # -n as haploid sample size (what about diploid?)

            echo \"npstat time\" >> time_collection
            (time \$Fraca_dir\"npstat\" -n \$n_chr -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out2.stats) 2>> time_collection
            #will output \$out\".stats\" in default

            #NR: number of records variable (awk NR gives the total number of records being processed or line number)
            theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out2\".stats\") #h.mpileup.stats\" as input file #nucleotide diversity
            D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out2\".stats\")
             #Prior genetic difference between reference genome and population
            #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print

            echo \"snape time\" >> time_collection
            (time \$Fraca_dir\"snape-pooled\" -nchr \$n_chr -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='\$min_all'' > \$out2.snape) 2>> time_collection

            #3. VarScan
            echo \"VarScan time\" >> time_collection
            (time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value $sig > \$out2.varscan) 2>> time_collection
            #--p-value as threshold to call variant
            "
            write(f, cmd, "\n")
        end

        write(f, "gzip \$out.mpileup", "\n")

        #input as mpileup.gz
        for sig in 0.01:0.01:0.99
            lrtSnp=cquantile(Chisq(1), sig)
            # pp_snape=1-sig
            cmd ="
            out2=\"simu-$dpt-$nsamp-$sig\"
            #4. ngsPool
            echo \"julia time\" >> time_collection
            (time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out2.gz --nChroms \$n_chr --lrtSnp $lrtSnp) 2>> time_collection
            "
            # println(cmd)
            # println("")
            write(f, cmd, "\n")
        end

        #################################set qq as 0
        pre_cmd="
        cd ../
        mkdir power_0; cd power_0
        out=\"simu-$dpt-$nsamp\"
        echo \"\$out mpileup##################################\" >> time_collection
        (time Rscript \$ngsPool_dir\"simulMpileup_qq.R\" --qq 0 --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection

        gzip -d \$out.mpileup.gz
        #1.Popoolation
        echo \"Popoolation java\" >> time_collection
        (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out.sync) 2>> time_collection
        "
        write(f, pre_cmd, "\n")

        for sig in 0.01:0.01:0.99
            lrtSnp=cquantile(Chisq(1), sig)
            pp_snape=1-sig

            cmd ="
            out2=\"simu-$dpt-$nsamp-$sig\" #*: here the last element is sig instead of qq(preset true MAF)

            #parameters for snape or varscan
            #chr_pool=8 #Arabidopsis lyrata has 8 chromosomes
            n_chr=$(2*nsamp) #total chromosomes pooled = diploid * num of individuals
            min=4 #read depth should be at least 4
            max=100 #default maximum read depth
            min_qual=20 #base quality
            l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
            min_all=2 #minimum allele count for Snape, VarScan and NPStat
            # pp_snape=0.9 #posterior probability threshold for Snape

            #2. Snape
            #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
            #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
            # -n as haploid sample size (what about diploid?)

            echo \"npstat time\" >> time_collection
            (time \$Fraca_dir\"npstat\" -n \$n_chr -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out2.stats) 2>> time_collection
            #will output \$out\".stats\" in default

            #NR: number of records variable (awk NR gives the total number of records being processed or line number)
            theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out2\".stats\") #h.mpileup.stats\" as input file #nucleotide diversity
            D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out2\".stats\")
             #Prior genetic difference between reference genome and population
            #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print

            echo \"snape time\" >> time_collection
            (time \$Fraca_dir\"snape-pooled\" -nchr \$n_chr -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='\$min_all'' > \$out2.snape) 2>> time_collection

            #3. VarScan
            echo \"VarScan time\" >> time_collection
            (time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value $sig > \$out2.varscan) 2>> time_collection
            #--p-value as threshold to call variant
            "
            write(f, cmd, "\n")
        end

        write(f, "gzip \$out.mpileup", "\n")

        #input as mpileup.gz
        for sig in 0.01:0.01:0.99
            lrtSnp=cquantile(Chisq(1), sig)
            # pp_snape=1-sig
            cmd ="
            out2=\"simu-$dpt-$nsamp-$sig\"
            #4. ngsPool
            echo \"julia time\" >> time_collection
            (time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out2.gz --nChroms \$n_chr --lrtSnp $lrtSnp) 2>> time_collection
            "
            # println(cmd)
            # println("")
            write(f, cmd, "\n")
        end
        write(f, "cd ../", "\n")

    end
    close(f)

##data analysis
sigs=collect(0.01:0.01:0.99);
    dpts=[3, 5]; nsamps=[30, 50]

    cd("C:\\simu_raw\\power_wider\\power_rnd") #power
    power2 = metric_2mtx2(dpts, nsamps, sigs)
    # power2[1,1]

    cd("C:\\simu_raw\\power_wider\\power_0") #false positive rate
    FPR2 = metric_2mtx2(dpts, nsamps, sigs)

power2[1,1]#0.531
power2[1,2]#0.626
power2[2,1]#0.597
power2[2,2]#0.724

FPR2[1,1]#0.21
FPR2[1,2]#0.275
FPR2[2,1]#0.305
FPR2[2,2]#0.451

#single plot with one depth (1) and sample size (30)
scatter(FPR2[1,1], power2[1,1], size=(600,500), ylim=(0,1), xlim=(0,1),
    label=["ngsPool" "Snape" "Popoolation" "VarScan"], legend=:bottom,
    fg_legend = :transparent, background_color_legend=:transparent,markerstrokewidth=0,
    grid=false, legendfontsize=10, markersize=5, marker=:hex, dpi=300)
    Plots.abline!(1,0, label="")
    plot!(ylabel="sample size 30\npower(true positive rate)")
    plot!(xlabel="False positive rate\ndepth=1")
savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\ROC_dpt1_sample30.png")


#the array of ROC curve plots of combinations with dpts=[3, 5]; nsamps=[30, 50]
allPlots = let allPlots=Array{Plots.Plot{Plots.GRBackend},2}(undef,2,2) #uni-nitialized plot array
    k=1
    for i in 1:2, j in 1:2

            if i==2 && j==2
                p=scatter(FPR2[i,j], power2[i,j], size=(600,500), ylim=(0,1), xlim=(0,1),
                    label=["ngsPool" "Snape" "Popoolation" "VarScan"],
                    legend=:best, markerstrokewidth=0, xlims=(0,0.7),
                    fg_legend = :transparent, background_color_legend=:transparent,
                    grid=false, legendfontsize=10, markersize=5, marker=:hex,
                    # xlabel="False positive rate", ylabel="power(true positive rate)"
                    )
            else
                p=scatter(FPR2[i,j], power2[i,j], size=(600,500), ylim=(0,1), xlim=(0,1),
                    label="",markerstrokewidth=0,
                    fg_legend = :transparent, background_color_legend=:transparent,
                    grid=false, legendfontsize=8, markersize=5, marker=:hex,
                    # xlabel="False positive rate", ylabel="power(true positive rate)"
                    )
            end
            if i==1 && j==1
                plot!(ylabel="sample size 30\npower(true positive rate)")
            elseif i==1 && j==2
                plot!(ylabel="sample size 50\npower(true positive rate)")
                plot!(xlabel="False positive rate\ndepth=3")
            elseif i==2 && j==2
                plot!(xlabel="False positive rate\ndepth=5")
            end

            plot!(title="($k)")
            Plots.abline!(1,0, label="")
        allPlots[i,j]=p
        k+=1
    end
    allPlots
end
# dpts=[3, 5]; nsamps=[30, 50]

# allPlots'
display(plot(allPlots...,layout=(2,2), dpi=300))
savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\ROC_dpt3_5_sample_30_50.png")


############ combine significance level both 0.001:0.001:0.1 and 0.1:0.01:0.99
dpts=[1]; nsamps=[30, 50];sigs=collect(0.001:0.001:0.1)
    cd("C:\\simu_raw\\power_wider\\power_rnd") #power
    power22 = metric_2mtx2(dpts, nsamps, sigs)

cd("C:\\simu_raw\\power_wider\\power_0") #false positive rate
    FPR22 = metric_2mtx2(dpts, nsamps, sigs)

plot(power2[1,1]) #very small change compared to ngsPool
#vertically concacenate the matrixes with significance level 0.001:0.1 to the matrixes with 0.1:0.99
all_power, all_FPR = let all_power=Array{Array{Float64,2},2}(undef,2,2), all_FPR=Array{Array{Float64,2},2}(undef,2,2)
    for i in 1:2, j in 1:2
        # println(size(all_power))
        # println(size(all_FPR))
        all_power[i,j]=vcat(power[i,j], power2[i,j])
        all_FPR[i,j]=vcat(FPR[i,j], FPR2[i,j])
    end
    all_power, all_FPR
end

power2=all_power
FPR2=all_FPR
#to reuse the block of codes of collecting plots with array


allPlots = let allPlots=Array{Plots.Plot{Plots.GRBackend},2}(undef,2,2) #uni-nitialized plot array
    k=1
    for i in 1:2, j in 1:2

            if i==2 && j==2
                p=scatter(FPR2[i,j], power2[i,j], size=(600,500), ylim=(0,1), xlim=(0,1),
                    label=["ngsPool" "Snape" "Popoolation" "VarScan"],
                    legend=:best, markerstrokewidth=0, xlims=(0,0.7),
                    fg_legend = :transparent, background_color_legend=:transparent,
                    grid=false, legendfontsize=10, markersize=5, marker=:hex,
                    # xlabel="False positive rate", ylabel="power(true positive rate)"
                    )
            else
                p=scatter(FPR2[i,j], power2[i,j], size=(600,500), ylim=(0,1), xlim=(0,1),
                    label="",markerstrokewidth=0,
                    fg_legend = :transparent, background_color_legend=:transparent,
                    grid=false, legendfontsize=8, markersize=5, marker=:hex,
                    # xlabel="False positive rate", ylabel="power(true positive rate)"
                    )
            end
            if i==1 && j==1
                plot!(ylabel="sample size 30\npower(true positive rate)")
            elseif i==1 && j==2
                plot!(ylabel="sample size 50\npower(true positive rate)")
                plot!(xlabel="False positive rate\ndepth=3")
            elseif i==2 && j==2
                plot!(xlabel="False positive rate\ndepth=5")
            end

            plot!(title="($k)")
            Plots.abline!(1,0, label="")
        allPlots[i,j]=p
        k+=1
    end
    allPlots
end
# dpts=[3, 5]; nsamps=[30, 50]

# allPlots'
display(plot(allPlots...,layout=(2,2), dpi=300))
savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\ROC_dpt1_sample_30_50.png")

# The significance levels contain one sequence from 0.001 to 0.1 in the unit of 0.001
# and another from 0.01 to 0.99 with the spacing of 0.01. The depths and sample sizes
# presented in this figure are {3, 5} and {30, 50}, which are annotated accordingly
# to the row and column of subplots. The diagonal line represents guesses made randomly.
# The further the ROC curve goes beyond the line, the better quality is suggested.

# The point of Snape tends to be slightly lower than the ROC curve of ngsPool.
# This also applies between Snape and ngsPool at relatively higher depth and
# sample size in Figure 4(4). In other cases, ngsPool estimates are associated
# with a higher power at a same false positive rate. The ROC sketch of VarScan
# can be seen close to the y axis. The extremely low false positive rates can
# be caused by its rare SNP callings.


################################################################################
# sample size: 30, 50, 100
# depth per individual: 1,3,5,10,20
# sigs = [0.001, 0.01, 0.1, (0.2:0.1:0.9...,)...]

cd("C:\\simu_raw\\pipeline") #qq-0
    using Statistics, Distributions
    f=open("power_FPR_p2", "w") #no preset true MAF
        pre_cmd=
    "
    cd /mnt/c/simu_raw/ #fixed population allele frequency
    mkdir power_p2; cd power_p2
    Popoo_dir=\"/mnt/H/julia/ngsJulia/ngsPool/popoolation2_1201/\"
    Fraca_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\" #snape and varscan not written but used by Fraca
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    echo \"\" > time_collection #open a file to collect time information
    "
    write(f, pre_cmd, "\n")
    # for dpt in [1,3,5], nsamp in [30, 50]
    for dpt in [1,3,5,10,20], nsamp in [30, 50, 100]
        # sig in [0.1, 0.05, 0.01, 0.005, 0.001] #more points needed for ROC curve
    # for dpt in [1], nsamp in [30], qq in [1e-4] #test
        # sig=0.05

        pre_cmd="
        mkdir power_rnd; cd power_rnd
        out=\"simu-$dpt-$nsamp\"
        echo \"\$out mpileup##################################\" >> time_collection
        (time Rscript \$ngsPool_dir\"simulMpileup.R\" --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection
        "
        write(f, pre_cmd, "\n")

        for sig in [0.001, 0.01, 0.1, (0.2:0.1:0.9...,)...]
            lrtSnp=cquantile(Chisq(1), sig)
            pp_snape=1-sig

            cmd ="
            out2=\"simu-$dpt-$nsamp-$sig\" #*: here the last element is sig instead of qq(preset true MAF)
            gzip -d \$out.mpileup.gz
            #1.Popoolation
            echo \"Popoolation java\" >> time_collection
            (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out2.sync) 2>> time_collection

            #parameters for snape or varscan
            #chr_pool=8 #Arabidopsis lyrata has 8 chromosomes
            n_chr=$(2*nsamp) #total chromosomes pooled = diploid * num of individuals
            min=4 #read depth should be at least 4
            max=100 #default maximum read depth
            min_qual=20 #base quality
            l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
            min_all=2 #minimum allele count for Snape, VarScan and NPStat
            # pp_snape=0.9 #posterior probability threshold for Snape

            #2. Snape
            #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
            #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
            # -n as haploid sample size (what about diploid?)

            echo \"npstat time\" >> time_collection
            (time \$Fraca_dir\"npstat\" -n \$n_chr -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out2.stats) 2>> time_collection
            #will output \$out\".stats\" in default

            #NR: number of records variable (awk NR gives the total number of records being processed or line number)
            theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out2\".stats\") #h.mpileup.stats\" as input file #nucleotide diversity
            D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out2\".stats\")
             #Prior genetic difference between reference genome and population
            #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print

            echo \"snape time\" >> time_collection
            (time \$Fraca_dir\"snape-pooled\" -nchr \$n_chr -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='\$min_all'' > \$out2.snape) 2>> time_collection

            #3. VarScan
            echo \"VarScan time\" >> time_collection
            (time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value $pp_snape > \$out2.varscan) 2>> time_collection
            #--p-value as threshold to call variant

            #4. ngsPool
            gzip \$out.mpileup
            echo \"julia time\" >> time_collection
            (time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out2.gz --nChroms \$n_chr --lrtSnp $lrtSnp) 2>> time_collection
            "
            # println(cmd)
            # println("")
            write(f, cmd, "\n")
        end

        #################################set qq as 0
        pre_cmd="
        cd ../
        mkdir power_0; cd power_0
        out=\"simu-$dpt-$nsamp\"
        echo \"\$out mpileup##################################\" >> time_collection
        (time Rscript \$ngsPool_dir\"simulMpileup_qq.R\" --qq 0 --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection
        "
        write(f, pre_cmd, "\n")

        for sig in [0.001, 0.01, 0.1, (0.2:0.1:0.9...,)...]
            lrtSnp=cquantile(Chisq(1), sig)
            pp_snape=1-sig

            cmd ="
            out2=\"simu-$dpt-$nsamp-$sig\" #*: here the last element is sig instead of qq(preset true MAF)
            gzip -d \$out.mpileup.gz
            #1.Popoolation
            echo \"Popoolation java\" >> time_collection
            (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out2.sync) 2>> time_collection

            #parameters for snape or varscan
            #chr_pool=8 #Arabidopsis lyrata has 8 chromosomes
            n_chr=$(2*nsamp) #total chromosomes pooled = diploid * num of individuals
            min=4 #read depth should be at least 4
            max=100 #default maximum read depth
            min_qual=20 #base quality
            l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
            min_all=2 #minimum allele count for Snape, VarScan and NPStat
            # pp_snape=0.9 #posterior probability threshold for Snape

            #2. Snape
            #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
            #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
            # -n as haploid sample size (what about diploid?)

            echo \"npstat time\" >> time_collection
            (time \$Fraca_dir\"npstat\" -n \$n_chr -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out2.stats) 2>> time_collection
            #will output \$out\".stats\" in default

            #NR: number of records variable (awk NR gives the total number of records being processed or line number)
            theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out2\".stats\") #h.mpileup.stats\" as input file #nucleotide diversity
            D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out2\".stats\")
             #Prior genetic difference between reference genome and population
            #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print

            echo \"snape time\" >> time_collection
            (time \$Fraca_dir\"snape-pooled\" -nchr \$n_chr -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='\$min_all'' > \$out2.snape) 2>> time_collection

            #3. VarScan
            echo \"VarScan time\" >> time_collection
            (time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value $sig > \$out2.varscan) 2>> time_collection
                                                                                                                        #pp_snape amended to be sig after generating data(might be the reason why no much difference made among sigs in varscan)
            #--p-value as threshold to call variant

            #4. ngsPool
            gzip \$out.mpileup
            echo \"julia time\" >> time_collection
            (time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out2.gz --nChroms \$n_chr --lrtSnp $lrtSnp) 2>> time_collection
            "
            # println(cmd)
            # println("")
            write(f, cmd, "\n")
        end
        write(f, "cd ../", "\n")
    end
    close(f)

#

#####with wider spectrum of significance level#################################
sigs=[0.001, 0.01, 0.1, (0.2:0.1:0.9...,)...]
    dpts=[1,3,5,10,20]
    nsamps=[30, 50, 100]

    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")

    cd("C:\\simu_raw\\power_p2\\power_rnd") #power
    power2 = metric_2mtx2(dpts, nsamps, sigs)
    # power2[1,1]

    cd("C:\\simu_raw\\power_p2\\power_0") #false positive rate
    FPR2 = metric_2mtx2(dpts, nsamps, sigs)

allPlots = let allPlots=Array{Plots.Plot{Plots.GRBackend},2}(undef,5,3) #uni-nitialized plot array
    k=1
    for i in 1:5, j in 1:3

            if i==5 && j==3
                p=scatter(FPR2[i,j], power2[i,j], size=(600,500), ylim=(-0.1,1.1), xlim=(-0.1,1),
                    label=["ngspool" "snape" "P2" "varscan"],
                    legend=:bottomleft, markerstrokewidth=0, xlims=(0,0.7),
                    fg_legend = :transparent, background_color_legend=:transparent,
                     legendfontsize=11, markersize=5, marker=:hex
                    # xlabel="False positive rate", ylabel="power(true positive rate)"
                    )
            else
                p=scatter(FPR2[i,j], power2[i,j], size=(300,300), ylim=(-0.1,1.1), xlim=(-0.1,1),
                    label="",markerstrokewidth=0,
                    fg_legend = :transparent, background_color_legend=:transparent,
                    legendfontsize=8, markersize=5, marker=:hex,
                    # xlabel="False positive rate", ylabel="power(true positive rate)"
                    )
            end

            if i==1
                plot!(ylabel="sample size $(nsamps[j])\npower(true positive rate)")
            end

            if j==3
                plot!(xlabel="False positive rate\ndepth $(dpts[i])")
            end

            Plots.abline!(1,0, label="")
        allPlots[i,j]=p
        k+=1
    end
    allPlots
end
allPlots[1,1]
# dpts=[3, 5]; nsamps=[30, 50]

# allPlots'
using Plots.PlotMeasures
    display(plot(allPlots...,layout=(3,5), dpi=300, size=(1500, 1000), margin=13mm))

power2_stack=zeros(0,4)
for i in 1:5, j in 1:3


savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\ROC_p2.png")

##compare script_power2 and script_power##########################
#using 1/N as input MAF frequency
include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
# dpts=[1,3,5,10,20]
    sigs=collect(0.01:0.01:0.99);
    dpts=[1]; nsamps=[30]
    # nsamps=[30, 50, 100]

    # the name of .sync with significance level
    cd("C:\\simu_raw\\power_wider_qqVector\\power_rnd") #power
    power2 = metric_2mtx2(dpts, nsamps, sigs)
    # power2[1,1]

    # cd("C:\\simu_raw\\power_wider_qqVector\\power_0") #false positive rate
    # FPR2 = metric_2mtx2(dpts, nsamps, sigs)

#using FIXED 10^n as input MAF
    cd("C:\\simu_raw\\power_wider\\power_rnd") #power
    power22 = metric_2mtx2(dpts, nsamps, sigs)
    # power2[1,1]

    cd("C:\\simu_raw\\power_wider\\power_0") #false positive rate
    FPR22 = metric_2mtx2(dpts, nsamps, sigs)

#only compare the ngsPool
function stack2(a::Array)
    n_row=size(a,1)
    n_col=size(a,2)
    inner_col=size(a[1,1],2)
    new=zeros(0, inner_col)
    for i in 1:n_row, j in 1:n_col
        new = vcat(a[i, j], new)
    end
    return new
end


power2_stack=stack2(power2)
# FPR2_stack=stack2(FPR2)

power22_stack=stack2(power22)
FPR22_stack=stack2(FPR22)

#useless block
    # plot(power2_stack[:,1], power22_stack[:,1], ylim=(0,0.7), xlim=(0,0.7),
    #     legend=:bottom, label="", xlabel="power (MAF input probability as 1/n)", ylabel="power (MAF input probability as 10^(-4:-1) )") #similar but not the same
    #     Plots.abline!(1,0, line=:dash, label="")
    #
    # scatter(FPR2_stack[:,1], power2_stack[:,1], label="10^sequence",
    #     xlabel="False positive rate", ylabel="Power", grid=false, legend=:bottomright,
    #     dpi=300, size=(500,500), ylim=(0,0.8), xlim=(0,0.8))
    #     scatter!(FPR22_stack[:,1], power22_stack[:,1], alpha=0.4, markersize=1, label="1/sequence", color=:grey)
    #     Plots.abline!(1,0,line=:dash, label="")
    # #they are exactly the same #so input MAF ferquency doesn't have a influence in the power
    #
    #
    # savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\input MAF seqeuncy.png")
    #
    #
    # scatter(FPR2_stack[:,2], power2_stack[:,2]) #, label="10^sequence",
    #     scatter!(FPR22_stack[:,2], power22_stack[:,2], alpha=0.3)#, alpha=0.4, markersize=1, label="1/sequence", color=:grey)
    #     #totally overlap
    #
    # scatter(FPR2_stack[:,3], power2_stack[:,3]) #, label="10^sequence",
    #     scatter!(FPR22_stack[:,3], power22_stack[:,3], alpha=0.3)#, alpha=0.4, markersize=1, label="1/sequence", color=:grey)
    #
    # scatter(FPR2_stack[:,4], power2_stack[:,4]) #, label="10^sequence",
    #     scatter!(FPR22_stack[:,4], power22_stack[:,4], alpha=0.3)#, alpha=0.4, markersize=1, label="1/sequence", color=:grey)
    #     #also overlap
