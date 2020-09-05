#Script for power and FPR

#using fixed MAFs (simulMpileup_qqVector.R)



#################### using simulMpileup_qqVector.R ###########################

#significance from 0.1 to 0.99 ###using simulMpileup_qqVector instead of simulMpileup (different input MAF)
 #script "power_FPR_wider_qqVector"; directory: "power_wider_qqVector"
cd("C:\\simu_raw\\pipeline") #qq-0
    using Statistics, Distributions
    f=open("power_FPR_wider_qqVector", "w") #no preset true MAF
        pre_cmd=
    "
    cd /mnt/c/simu_raw/ #fixed population allele frequency
    mkdir power_wider_qqVector; cd power_wider_qqVector
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
        (time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection

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

        write(f, "gzip \$out.mpileup", "\n") #input as mpileup.gz

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
        # pre_cmd="
        # cd ../
        # mkdir power_0; cd power_0
        # out=\"simu-$dpt-$nsamp\"
        # echo \"\$out mpileup##################################\" >> time_collection
        # (time Rscript \$ngsPool_dir\"simulMpileup_qq.R\" --qq 0 --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection
        #
        # gzip -d \$out.mpileup.gz
        # #1.Popoolation
        # echo \"Popoolation java\" >> time_collection
        # (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out.sync) 2>> time_collection
        # "
        # write(f, pre_cmd, "\n")
        #
        # for sig in 0.01:0.01:0.99
        #     lrtSnp=cquantile(Chisq(1), sig)
        #     pp_snape=1-sig
        #
        #     cmd ="
        #     out2=\"simu-$dpt-$nsamp-$sig\" #*: here the last element is sig instead of qq(preset true MAF)
        #
        #     #parameters for snape or varscan
        #     #chr_pool=8 #Arabidopsis lyrata has 8 chromosomes
        #     n_chr=$(2*nsamp) #total chromosomes pooled = diploid * num of individuals
        #     min=4 #read depth should be at least 4
        #     max=100 #default maximum read depth
        #     min_qual=20 #base quality
        #     l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
        #     min_all=2 #minimum allele count for Snape, VarScan and NPStat
        #     # pp_snape=0.9 #posterior probability threshold for Snape
        #
        #     #2. Snape
        #     #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
        #     #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
        #     # -n as haploid sample size (what about diploid?)
        #
        #     echo \"npstat time\" >> time_collection
        #     (time \$Fraca_dir\"npstat\" -n \$n_chr -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out2.stats) 2>> time_collection
        #     #will output \$out\".stats\" in default
        #
        #     #NR: number of records variable (awk NR gives the total number of records being processed or line number)
        #     theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out2\".stats\") #h.mpileup.stats\" as input file #nucleotide diversity
        #     D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out2\".stats\")
        #      #Prior genetic difference between reference genome and population
        #     #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print
        #
        #     echo \"snape time\" >> time_collection
        #     (time \$Fraca_dir\"snape-pooled\" -nchr \$n_chr -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='\$min_all'' > \$out2.snape) 2>> time_collection
        #
        #     #3. VarScan
        #     echo \"VarScan time\" >> time_collection
        #     (time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value $sig > \$out2.varscan) 2>> time_collection
        #     #--p-value as threshold to call variant
        #     "
        #     write(f, cmd, "\n")
        # end
        #
        # write(f, "gzip \$out.mpileup", "\n")
        #
        # #input as mpileup.gz
        # for sig in 0.01:0.01:0.99
        #     lrtSnp=cquantile(Chisq(1), sig)
        #     # pp_snape=1-sig
        #     cmd ="
        #     out2=\"simu-$dpt-$nsamp-$sig\"
        #     #4. ngsPool
        #     echo \"julia time\" >> time_collection
        #     (time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out2.gz --nChroms \$n_chr --lrtSnp $lrtSnp) 2>> time_collection
        #     "
        #     # println(cmd)
        #     # println("")
        #     write(f, cmd, "\n")
        # end
        # write(f, "cd ../", "\n")

    end
    close(f)

#significance from 0.001:0.001:0.1 ###using simulMpileup_qqVector instead of simulMpileup (different input MAF)
    #just ngsPool
    # script "power_FPR_qqVector2"; directory "power_wider_qqVector"
cd("C:\\simu_raw\\pipeline") #qq-0
    using Statistics, Distributions
    f=open("power_FPR_qqVector2", "w") #no preset true MAF
        pre_cmd=
    "
    cd /mnt/c/simu_raw/ #fixed population allele frequency
    mkdir power_wider_qqVector; cd power_wider_qqVector
    Popoo_dir=\"/mnt/H/julia/ngsJulia/ngsPool/popoolation2_1201/\"
    Fraca_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\" #snape and varscan not written but used by Fraca
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    echo \"\" > time_collection #open a file to collect time information
    "
    write(f, pre_cmd, "\n")
    # for dpt in [3,5], nsamp in [30, 50]
    for dpt in [1], nsamp in [30,50]
        # sig in [0.1, 0.05, 0.01, 0.005, 0.001] #more points needed for ROC curve
    # for dpt in [1], nsamp in [30], qq in [1e-4] #test
        # sig=0.05

        pre_cmd="
        mkdir power_rnd; cd power_rnd
        # out=\"simu-$dpt-$nsamp\"
        # echo \"\$out mpileup##################################\" >> time_collection
        # (time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection
        #
        # gzip -d \$out.mpileup.gz
        # #1.Popoolation
        # echo \"Popoolation java\" >> time_collection
        # (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out.sync) 2>> time_collection
        "
        write(f, pre_cmd, "\n")
        #
        # for sig in 0.001:0.001:0.1
        #     lrtSnp=cquantile(Chisq(1), sig)
        #     pp_snape=1-sig
        #
        #     cmd ="
        #     out2=\"simu-$dpt-$nsamp-$sig\" #*: here the last element is sig instead of qq(preset true MAF)
        #
        #     #parameters for snape or varscan
        #     #chr_pool=8 #Arabidopsis lyrata has 8 chromosomes
        #     n_chr=$(2*nsamp) #total chromosomes pooled = diploid * num of individuals
        #     min=4 #read depth should be at least 4
        #     max=100 #default maximum read depth
        #     min_qual=20 #base quality
        #     l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
        #     min_all=2 #minimum allele count for Snape, VarScan and NPStat
        #     # pp_snape=0.9 #posterior probability threshold for Snape
        #
        #     #2. Snape
        #     #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
        #     #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
        #     # -n as haploid sample size (what about diploid?)
        #
        #     echo \"npstat time\" >> time_collection
        #     (time \$Fraca_dir\"npstat\" -n \$n_chr -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out2.stats) 2>> time_collection
        #     #will output \$out\".stats\" in default
        #
        #     #NR: number of records variable (awk NR gives the total number of records being processed or line number)
        #     theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out2\".stats\") #h.mpileup.stats\" as input file #nucleotide diversity
        #     D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out2\".stats\")
        #      #Prior genetic difference between reference genome and population
        #     #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print
        #
        #     echo \"snape time\" >> time_collection
        #     (time \$Fraca_dir\"snape-pooled\" -nchr \$n_chr -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '$pp_snape'' | awk '\$5 >='\$min_all'' > \$out2.snape) 2>> time_collection
        #
        #     #3. VarScan
        #     echo \"VarScan time\" >> time_collection
        #     (time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value $sig > \$out2.varscan) 2>> time_collection
        #     #--p-value as threshold to call variant
        #     "
        #     write(f, cmd, "\n")
        # end
        #
        # write(f, "gzip \$out.mpileup", "\n") #input as mpileup.gz

        for sig in 0.001:0.001:0.1
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

    end
    close(f)

#######################
# sig in [0.001, 0.01, 0.1] ---------"power_FPR_narrower" ###using simulMpileup_qqVector instead of simulMpileup (different input MAF)
# rep in 1:100 #replicates
 #script "power_FPR_narrower"; directory: "power_narrower"
cd("C:\\simu_raw\\pipeline")
    using Statistics, Distributions
    f=open("power_FPR_narrower", "w") #no preset true MAF
        pre_cmd=
    "
    cd /mnt/c/simu_raw/ #fixed population allele frequency
    mkdir power_narrower; cd power_narrower
    Popoo_dir=\"/mnt/H/julia/ngsJulia/ngsPool/popoolation2_1201/\"
    Fraca_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\" #snape and varscan not written but used by Fraca
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    echo \"\" > time_collection #open a file to collect time information
    # mkdir power_rnd;
    cd power_rnd
    "
    write(f, pre_cmd, "\n")
    # for dpt in [3,5], nsamp in [30, 50] #simulated
    # # for dpt in [1], nsamp in [30] #simulated
    for dpt in [1], nsamp in [50], rep in 1:100
        pre_cmd="
        out=\"simu-$dpt-$nsamp-$rep\"
        echo \"\$out mpileup##################################\" >> time_collection
        (time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out \$out.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual 20 --ksfs 1 --ne 10000 --pool | gzip > \$out.mpileup.gz) 2>> time_collection

        gzip -d \$out.mpileup.gz
        #1.Popoolation
        echo \"Popoolation java\" >> time_collection
        (time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out.sync) 2>> time_collection
        "
        write(f, pre_cmd, "\n")

        for sig in [0.001, 0.01, 0.1]
            lrtSnp=cquantile(Chisq(1), sig)
            pp_snape=1-sig

            cmd ="
            out2=\"simu-$dpt-$nsamp-$rep-$sig\" #*: here the last element is sig instead of qq(preset true MAF)

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
        for sig in [0.001, 0.01, 0.1]
            lrtSnp=cquantile(Chisq(1), sig)
            # pp_snape=1-sig
            cmd ="
            out=\"simu-$dpt-$nsamp-$rep\"
            out2=\"simu-$dpt-$nsamp-$rep-$sig\"
            #4. ngsPool
            echo \"julia time\" >> time_collection
            (time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out2.gz --nChroms \$n_chr --lrtSnp $lrtSnp) 2>> time_collection
            "
            # println(cmd)
            # println("")
            write(f, cmd, "\n")
        end
    end
    close(f)

#so little difference in the number of called sites is introduced if fixing probability vector of MAFs (qqVector)!!!
#replicates has no meaning after fixing qqVector(though MAF varies)
cd("C:\\simu_raw\\power_narrower\\power_rnd") #power
p2_FPR = let p2_FPR=zeros(100)
    # for dpt in dpts, nsamp in nsamps, rep in reps
    # i=0
    # for d in 1:2, s in 1:2
        # i+=1
        dpt=1
        nsamp=50
        for r in 1:100
            # name="simu-$(dpts[d])-$(nsamps[s])-$(reps[r])"
            name="simu-$dpt-$nsamp-$r"
            df_p2 = read_p2(name)
            df_p2 = filter(:popoo => x -> !iszero(x), df_p2)
            p2_FPR[r]=nrow(df_p2)
        end
    # end
    p2_FPR
end
plot(p2_FPR) # called sites number are all 444

################################################################################
#why only ngsPool have different reaction accoring to significance level?
    #the snp calling result visually only vary when depth is very low (dpt=1)
    # so test under depth 1 and sample size 50 (according to Figure 4)

dpt=1; nsamp=50;sigs=[0.001, 0.1]; replicates=1:100
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
    cd("C:\\simu_raw\\power_narrower\\power_rnd") #power
    power3 = metric_mtx3(dpt, nsamp, sigs, replicates)

cd("C:\\simu_raw\\power_narrower\\power_0") #false positive rate
FPR3 = metric_mtx3(dpt, nsamp, sigs, replicates)

scatter(FPR3, power3)
dfs=[df_ngspool, df_snape, df_p2,  df_varscan]

plot(power3) #only significance level made a difference in the number of sites called
    #all the replicate should have no impact???!

plot(FPR3)

#######################################
#compare 4 software in a wider range of depth and sample size

#qual in [20, 30], dpt in [1,3,5,10,20], nsamp in [30, 50, 100]
 #mean base quality (for simulating the base quality in mpileup file) =20 in sys20/;
 #mean base quality = 30 in sys30/ directory

# 1) produce mpileup files
cd("C:\\simu_raw\\pipeline")
    f=open("sys_mpileup", "w")
    pre_cmd=
    "cd /mnt/c/simu_raw/
    mkdir sys
    cd sys
    mkdir sys20; mkdir sys30
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    "
    write(f, pre_cmd, "\n")

    for qual in [20, 30], dpt in [1,3,5,10,20], nsamp in [30, 50, 100]
        cmd=
        "
    # name=\"sys$qual/simu-$dpt-$nsamp\" >> names_sys #not used
    time Rscript \$ngsPool_dir\"simulMpileup_qqVector.R\" --out \$name.txt --copy 2x$nsamp --sites 1000 --depth $dpt --qual $qual --ksfs 1 --ne 10000 --pool | gzip > \$name.mpileup.gz
    "
	#qqVector: 1e-4:length=1000, 1e-1
        println(cmd)
        println("")
        write(f, cmd, "\n")
    end
    close(f)
#each simulation takes around 20s

# 2) inference by programs
#all significance level are set as 0.05 here
cd("C:\\simu_raw\\pipeline")
    f=open("sys_script", "w")
        pre_cmd=
    "
    cd /mnt/c/simu_raw/sys/
    # other_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\" #other program
    Popoo_dir=\"/mnt/H/julia/ngsJulia/ngsPool/popoolation2_1201/\"
    Fraca_dir=\"/mnt/H/julia/ngsJulia/ngsPool/pooled_WGS/\" #snape and varscan not written but used by Fraca
    ngsPool_dir=\"/mnt/H/julia/ngsJulia/ngsPool/\"
    "
    write(f, pre_cmd, "\n")

    for qual in [20, 30], dpt in [1,3,5,10,20], nsamp in [30, 50, 100]
        #qual=20; dpt=1;nsamp=30

	        # sig=0.05
	        # lrtSnp=cquantile(Chisq(1), sig)
	        # pp_snape=1-sig
		#threshold set for three programs (follows their default setting):
		 #varscan: 0.05
		 #ngsPool: 0.01
		 #snape: 0.1
		 #Popoolaion: no (said to be bad at SNP calling, anything do with no threshold?)

        cmd="
        out=\"sys$qual/simu-$dpt-$nsamp\"
        gzip -d \$out.mpileup.gz
        #1.Popoolation
    echo \"Popoolation############################java time\"
        time java -ea -Xmx7g -jar \$Popoo_dir\"mpileup2sync.jar\" --input \$out.mpileup --fastq-type sanger --min-qual 20 --threads 8 --output \$out.sync

        #parameters for snape or varscan
        chr_pool=8 #Arabidopsis lyrata has 8 chromosomes (without information specified from Fraca, assume they're all pooled)
        min=4 #read depth should be at least 4
        max=100 #default maximum read depth
        min_qual=20 #base quality
        l_npstat=1000 #base pair(unit, bp) of the windows for NPStat
        min_all=2 #minimum allele count for Snape, VarScan and NPStat
        pp_snape=0.9 #posterior probability threshold for Snape

        #2. Snape
        #The nucleotide diversity and the genetic differentiation from the reference genome that are needed
        #to set prior probabilities in the Bayesian model of Snape were calculated by NPStat [37].
        # -n as haploid sample size (what about diploid?)
    echo \"npstat##############################time\"
        time \$Fraca_dir\"npstat\" -n $(2*nsamp) -l \$l_npstat -mincov \$min -maxcov \$max -minqual \$min_qual -nolowfreq \$min_all \$out.mpileup > \$out.stats
        #will output \$out\".stats\" in default

        #NR: number of records variable (awk NR gives the total number of records being processed or line number)
        theta=\$(awk '{sum+=\$6} END {print sum/NR}' \$out\".mpileup.stats\") #h.mpileup.stats\" as input file #nucleotide diversity
        D=\$(awk '{sum+=\$13} END {print sum/NR}' \$out\".mpileup.stats\")
         #Prior genetic difference between reference genome and population
        #echo \"#Chromosomes \$out, theta is \$theta, D is \$D\" > \$out\"_stat\" #won't print

    echo \"snape##############################time\"
        time \$Fraca_dir\"snape-pooled\" -nchr $(2*nsamp) -theta \$theta -D \$D -fold unfolded -priortype informative < \$out.mpileup | awk '\$9 > '0.05'' | awk '\$5 >='\$min_all'' > \$out.snape

        #3. VarScan
        echo \"VarScan time\"
        time java -Xmx2g -jar \$Fraca_dir\"VarScan.v2.3.9.jar\" mpileup2snp \$out\".mpileup\" --min-coverage \$min --min-avg-qual \$min_qual --min-reads2 \$min_all --p-value 0.9 > \$out.varscan
        #--p-value as threshold of posterior probability to call variant

        #4. ngsPool
        gzip \$out.mpileup
    echo \"julia##############################time\"
        time julia \$ngsPool_dir\"ngsPool.jl\" --fin \$out.mpileup.gz --fout \$out.gz --nChroms $(2*nsamp) --lrtSnp 7.879438576622414 #0.01 as the significance level
        "
        println(cmd)
        println("")
        write(f, cmd, "\n")
    end
    close(f)

    #elements in out_names are without ""

# 3) data analysis
#sig=0.05
qual in [20, 30], dpt in [1,3,5,10,20], nsamp in [30, 50, 100]
qual=20

# include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_plot.jl")
qual=20
	dpts=[1,3,5,10,20]
	nsamps=[30, 50, 100]


cd("C:\\simu_raw\\sys\\sys20")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
	# #what I am interested just some combinations
	# dpt=1;nsamp=100 #snape no varscan and snape
	# # dpt=3;nsamp=30 #
	# #
	# # dpt=1;nsamp=30
	# # dpt=20;nsamp=100 #all fits so well, popoolation vertical has some pattern between 1e-2.5, and 1e-3
	# xyplots = sys_xyplot(dpt,nsamp)

	qual=20
	dpts=[1,3,5,10,20]
	nsamps=[30, 50, 100]
	rmse_plots, mpe_plots = sys_all_boxplots(nsamps, dpts) #LoadError: syntax: invalid named tuple element "SSAValue(61137)" ##????


################################################################################
#the probability of a potential MAF plotted with the potential MAF (the corrected version)

using Plots.PlotMeasures
    Ne=10000
    poly=1 ./ collect(1:(Ne-1)) #probability for polymorphic sites
    poly=poly/sum(poly) #normalize it
    all_prob=[0; poly; 0] #prabability including non-polymorphic and fixed sites-all sites
    MAFs_range=collect(0:Ne) ./ Ne #the ragne of potential MAF

    scatter(MAFs_range, all_prob, yscale=:log10, ylims=(1e-5, 0.11), label="sampling", dpi=300,margin=4mm,
    xlabel="MAF for sampling", ylabel="Probability of sampling MAF", markersize=2.5, markershape=:hex)

    #around dozen => therefore might not make a difference even with 100 replicates
    scatter!(10 .^ range(-1,-4,length=1000), fill(1/1000, 1000), label="10^(-1:-4)", markersize=2.5, markershape=:cross)

savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\why_ROC_initially_points_correct.png")
