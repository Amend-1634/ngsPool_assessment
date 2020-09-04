using Statistics, Distributions, Plots, Plots.PlotMeasures, DataFrames, StatsPlots

#annotation:
	#dpts: depths
	#sigs: significance levels
	#nsamps: sample sizes
	#xyplot: x as true MAFs, y as estimated MAFs
#the commented section following each function is the testing/example of this function


#input the vector of true MAFs and estimated MAFs (the latter with missing values)
#output the root mean square error vector without the entry of missing values
function RMSE_missing(true_maf, maf)
	#extract only sites with data
	nrow_use=findall(x->!ismissing(x), maf)
	maf=maf[nrow_use,:]; true_maf=true_maf[nrow_use,:]
	# rmse_vec = sqrt(sum( (maf .- true_maf).^2 ) ./ length(maf) )
	rmse_vec = sqrt.( (maf .- true_maf).^2 ) ./ length(maf)  #sum(replicate) but not among sites
	rmse_ave = mean(rmse_vec)
	return rmse_vec, rmse_ave, nrow_use
end

#calculate Mean percentage error
    #delete the entry with missing value or true_maf=0 (this might not be the standard MPE methods)
function MPE_missing(true_maf, maf) #input two vectors
    #extract only sites with data
    nrow_use=[findall(x->!ismissing(x), maf); findall(x->!iszero(x), true_maf)]

    maf=maf[nrow_use,:]; true_maf=true_maf[nrow_use,:]
    mpe_vec = (maf .- true_maf) ./ true_maf
    mpe_ave = mean(mpe_vec)
    return mpe_vec, mpe_ave, nrow_use
end


####### New idea of the RMSE and MPE (set all sites not called as 0)############
	#input the outerjoined output dataframe #missing value already replaced by 0
	 #(columns as chrom, position, true, estimation of each software(each col as a software))
function rmse_mpe(all_sites) #input as dataframe
    nrow_all=nrow(all_sites) #size() not working for df
	ncol_all=ncol(all_sites)-3 #
    rmse_mtx=Array{Union{Missing, Float64}}(undef, nrow_all, ncol_all)
    mpe_mtx=Array{Union{Missing, Float64}}(undef, nrow_all, ncol_all)
    rmse_ave_vec=[]
    mpe_ave_vec=[]

    for i in 4:ncol(all_sites)#1:3 chrom position true_maf estimated_maf...
        rmse_vec, rmse_ave, rmse_nrow_use = RMSE_missing(all_sites.true_maf, all_sites[:,i])
        mpe_vec, mpe_ave, mpe_nrow_use = MPE_missing(all_sites.true_maf, all_sites[:,i])

        rmse_mtx[rmse_nrow_use, [i-3]] .= rmse_vec
        mpe_mtx[mpe_nrow_use, [i-3]] .= mpe_vec

        push!(rmse_ave_vec, rmse_ave)
        push!(mpe_ave_vec, mpe_ave)
    end
    return rmse_mtx, mpe_mtx, rmse_ave_vec, mpe_ave_vec
end

#boxplot of four columns in a matrix (GSS, maximum likelihood estimation and expected estimation in ngsPool and Snape estimation)
	#allowing missing value of each column
function boxplot_missing(m)
    bp=boxplot(xticks=(1:4,["GSS_ngsPool" "freqMax_ngsPool" "freqE_ngsPool" "Snape"]),
	legend=false, dpi=300, fg_legend = :transparent, grid=false, xrotation=-15)

    for i in 1:size(m,2)
		bp=boxplot!(filter(!ismissing,m[:,i]))
		# label=["GSS_ngsPool" "freqMax_ngsPool" "freqE_ngsPool" "Snape"][i],
    end
    # display(bp)
    return bp
end

#input name and output plot
#output three plots: xyplot, 2 boxplots of rmse and mpe
function name_to_plot_sig(nsamp, dpt, sig)
	name="simu-$dpt-$nsamp"
	df_true = read_true(name)

	name="$name-$sig"
	df_snp=read_ngsPool_snp(name)
	df_snape_p, df_snape_nallele = read_snape(name)
	# println("df_snape_p")
	# println(head(df_snape_p))
	# println(typeof(df_snape_p))

	all_sites=outerjoin(df_true, df_snp, df_snape_p, on=[:chrom, :position])
	# all_sites = coalesce.(all_sites, 0) #all plots generated with this added is named with _nomissing
	# println(size(all_sites))

	all_sites2=all_sites[1:1000,:] #used to cut the x
	est2=convert(Matrix, all_sites2[:,[:GSS_ngsPool, :maf_snape]])
	xyplot = Plots.scatter(all_sites2.true_maf, est2, xscale=:log10, grid=false,
		title="depth $dpt x  p threshold $sig",
		label=["GSS_ngsPool" "snape"], legend=:topleft, fg_legend = :transparent,
		dpi=300, ylabel="Minor allele frequency(maf)", xlabel="Population maf") #,legend=false
		plot!(all_sites2.true_maf, all_sites2.true_maf, label="true")

		filt_gss = filter(x->!ismissing(x), all_sites2[:,:GSS_ngsPool])
		min_gss = minimum((filt_gss...,))
		# println(min_gss)
		filt_snape = filter(x->!ismissing(x), all_sites2[:,:maf_snape])
		min_snape = minimum((filt_snape...,))
		# println(min_snape)
		Plots.abline!(0, min_gss, line=:dash, color=:blue,label="") #abline(slope,intercept)
		Plots.abline!(0, min_snape, line=:dash, color=:orange,label="")
		plot!(title="(a)")

	rmse_mtx, mpe_mtx, rmse_ave_vec, mpe_ave_vec = rmse_mpe(all_sites)

	# println("rmse_mtx")
	# println(typeof(rmse_mtx))
	rmse_bp = boxplot_missing(rmse_mtx)
	plot!(ylabel="Root mean square error")
	plot!(title="(b)")
	mpe_bp = boxplot_missing(mpe_mtx)
	plot!(ylabel="Mean percentage error (%)")
	plot!(title="(c)")

	# Plots.scalefontsizes(15)
	display(Plots.plot(xyplot, rmse_bp, mpe_bp, layout=(1,3), size=(1200,500),
		margin=12mm))#tickfontsize=20, legendfontsize=20, guidefontsize=20))
	savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\3p-nsamp$nsamp-dpt$dpt-sig$sig.png")
	# savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\3p-nsamp$nsamp-dpt$dpt-sig$sig-nomissing.png")

	return xyplot, rmse_bp, mpe_bp
end

# name_to_plot_sig(50,1,0.05)

#output array of plots (so we can index the plot in interest out of it instead
 #of having a fixed layout of plot)
function all_3plots(nsamp, dpts, sigs) #3plots: xyplot, 2 boxplots of rmse and mpe
	n_sigs = length(sigs); n_dpts = length(dpts)
	xyplots = Array{Plots.Plot{Plots.GRBackend},2}(undef,n_sigs, n_dpts)
	rmse_plots = Array{Plots.Plot{Plots.GRBackend},2}(undef,n_sigs, n_dpts)
	mpe_plots = Array{Plots.Plot{Plots.GRBackend},2}(undef,n_sigs, n_dpts)

	for d in 1:n_dpts, s in 1:n_sigs
		sig=sigs[s]
		dpt=dpts[d]
		# lrtSnp=cquantile(Chisq(1), sig)
		# pp_snape=1-sig

		xyplot, rmse_bp, mpe_bp = name_to_plot_sig(nsamp, dpt, sig)
		xyplots[s,d]=xyplot
		rmse_plots[s,d]=rmse_bp
		mpe_plots[s,d]=mpe_bp
	end

	xy_all_p = plot(xyplots..., legend=false) #, layout=(n_sigs,n_dpts)

	rmse_all_p = plot(rmse_plots...) #, layout=(n_sigs,n_dpts)
	mpe_all_p = plot(mpe_plots...) #, layout=(n_sigs,n_dpts)

	return xyplots, rmse_plots, mpe_plots, xy_all_p, rmse_all_p, rmse_all_p
end

#all xyplots
function sig_all_xyplots(nsamp, dpts, sigs)
	n_sigs = length(sigs); n_dpts = length(dpts)
	xyplots = Array{Plots.Plot{Plots.GRBackend},2}(undef, n_sigs, n_dpts)

	for s in 1:n_sigs, d in 1:n_dpts
		name="simu-$(dpts[d])-$nsamp"
		df_true = read_true(name)

		name="$name-$(sigs[s])"
		df_snp=read_ngsPool_snp(name)
		df_snape_p, df_snape_nallele = read_snape(name)
		all_sites=outerjoin(df_true, df_snp, df_snape_p, on=[:chrom, :position])
		# println("all_sites")
		# println(typeof(all_sites))

		all_sites2=all_sites[1:1000,:] #used to cut the x

		est2=convert(Matrix, all_sites2[:,[:GSS_ngsPool, :maf_snape]])
		xyplot = Plots.scatter(all_sites2.true_maf, est2, xscale=:log10, grid=false,
			label="", xticksfontsize=12, yticksfontsize=12,
			# title=" $(dpts[d]) x  p  $(sigs[s])", titlefontsize=13,
			# label=["GSS_ngsPool" "snape"], legend=:topleft, fg_legend = :transparent,
			dpi=300) #,legend=false
			plot!(all_sites2.true_maf, all_sites2.true_maf, label="") # ="true"

			filt_gss = filter(x->!ismissing(x), all_sites2[:,:GSS_ngsPool])
			min_gss = minimum((filt_gss...,))
			# println(min_gss)
			filt_snape = filter(x->!ismissing(x), all_sites2[:,:maf_snape])
			min_snape = minimum((filt_snape...,))
			# println(min_snape)
			Plots.abline!(0, min_gss, line=:dash, color=:blue,label="") #abline(slope,intercept)
			Plots.abline!(0, min_snape, line=:dash, color=:orange,label="")

			if d==2 #the end row
				plot!(xlabel= "True maf \n significance $(sigs[s])")
			end

			if s==1
				plot!(ylabel="depth $(dpts[d])x \n Estimated maf", titlefontsize=12)
			end
			xyplots[s,d]=xyplot

			# if s==5 && d==2
			# 	plot!(label=["GSS_ngsPool" "snape" "true"], legend=:topleft)
			# end
	end
	display(plot(xyplots..., layout=(n_dpts, n_sigs), size=(1500,800),
		left_margin=13mm, bottom_margin=10mm))
	Plots.savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\compare-sig.png")

	return xyplots
end

# nsamp=50 #no 100 as snape doesn't have data
# dpts = [1,5]
# sigs = [0.1, 0.05, 0.01, 0.005, 0.001]
# xyplots = sig_all_plot(nsamp, dpts, sigs)


######compare significance level => all boxplot (RMSE, MPE)
function sig_all_boxplots(nsamp, dpts, sigs) #significance threshold
	n_sigs = length(sigs); n_dpts = length(dpts)
	rmse_plots = Array{Plots.Plot{Plots.GRBackend},2}(undef,n_sigs, n_dpts)
	mpe_plots = Array{Plots.Plot{Plots.GRBackend},2}(undef,n_sigs, n_dpts)

	for s in 1:n_sigs, d in 1:n_dpts
		name="simu-$(dpts[d])-$nsamp"
		df_true = read_true(name)

		name="$name-$(sigs[s])"
		df_snp=read_ngsPool_snp(name)
		df_snape_p, df_snape_nallele = read_snape(name)
		all_sites=outerjoin(df_true, df_snp, df_snape_p, on=[:chrom, :position])

		rmse_mtx, mpe_mtx, rmse_ave_vec, mpe_ave_vec = rmse_mpe(all_sites)

		rmse_bp = boxplot_missing(rmse_mtx)
		# plot!(ylabel="Root mean square error")
		if d==2 #the end row
			plot!(xlabel= "Program \n significance $(sigs[s])")
		end
		if s==1
			plot!(ylabel="depth $(dpts[d])x \n Root mean square error", titlefontsize=12)
		end

		mpe_bp = boxplot_missing(mpe_mtx)
		# plot!(ylabel="Mean percentage error")
		if d==2 #the end row
			plot!(xlabel= "Program \n significance $(sigs[s])")
		end

		if s==1
			plot!(ylabel="depth $(dpts[d])x \n Mean percentage error", titlefontsize=12)
		end

		rmse_plots[s,d] = rmse_bp
		mpe_plots[s,d] = mpe_bp
	end
	display(plot(rmse_plots..., layout=(n_dpts, n_sigs), size=(1500,800),
		left_margin=13mm, bottom_margin=10mm))

	display(plot(mpe_plots..., layout=(n_dpts, n_sigs), size=(1500,800),
		left_margin=13mm, bottom_margin=10mm))
	# savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\compare-sig.png")

	return rmse_plots, mpe_plots
end

# nsamp=50 #no 100 as snape doesn't have data
# dpts = [1,5]
# sigs = [0.1, 0.05, 0.01, 0.005, 0.001]
# rmse_plots, mpe_plots = sig_all_boxplots(nsamp, dpts, sigs)
#
# plot(rmse_plots..., layout=(2,5), size=(1500, 800))


##boxplot comparing the minor allele frequency of snape and ngsPool among significance level (x axis)
function ngsPool_snape_sig_maf(nsamp, dpts, sigs) #compare significance level
	# for d in 1:n_dpts,
	for dpt in dpts
		n_sigs=length(sigs)
		# gss_mtx=Array{Union{Missing, Float64}}(undef, 1000, n_sigs) #1000 sites
		# snape_mtx=Array{Union{Missing, Float64}}(undef, 1000, n_sigs)
		gss_df = DataFrame(Array{Union{Missing, Float64}}(undef, 1000, n_sigs))
		rename!(gss_df, string.(sigs))
		snape_df = DataFrame(Array{Union{Missing, Float64}}(undef, 1000, n_sigs))
		rename!(snape_df, string.(sigs))
		true_maf=[]

		for	s in 1:n_sigs

			# nsamp=50; dpt=1; sig=0.05
			name="simu-$dpt-$nsamp"
			df_true = read_true(name)

			name="$name-$(sigs[s])"
			df_snp=read_ngsPool_snp(name)
			df_snape_p, df_snape_nallele = read_snape(name)
			all_sites=outerjoin(df_true, df_snp, df_snape_p, on=[:chrom, :position])
			# if s == 1
			# 	true_maf = all_sites.true_maf
			# end

			gss_df[:,s] .= all_sites.GSS_ngsPool
			# println(all_sites.maf_snape)
			# println(size(all_sites.maf_snape))
			# println(size(snape_df[:,s]))

			snape_df[:,s] .= all_sites.maf_snape
			#matrix x-all sites, y-significance level
		end

		long_gss = stack(gss_df)
		# long_gss[!,"program"] = fill(1, size(long_gss,1))
		long_gss[!,"program"] = fill("GSS ", size(long_gss,1))
		long_snape = stack(snape_df)
		# long_snape[!,"program"] = fill(2, size(long_snape,1))
		long_snape[!,"program"] = fill("Snape", size(long_snape,1))
		long_df=vcat(long_gss, long_snape)
		filter!(:value => value -> !ismissing(value), long_df)
		variable=String.(long_df.variable)
		long_df.variable=variable

		p=StatsPlots.groupedboxplot(long_df.variable, long_df.value,
		group=long_df.program, grid=false,
		label=["GSS_ngsPool" "Snape"], legend=:outertopleft,alpha=0.9,
		fg_legend = :transparent, dpi=300, size=(500, 300),
		xlabel="Significance level", ylabel="Minor allele frequency")
		display(p)
		savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\boxplot-sig-nsamp$nsamp-dpt$dpt.png")

	end #dpt
end

# cd("C:\\simu_raw\\sys\\sys20")
# 	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
# 	include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
#
# name="simu-1-30"

#compare the estimations (each site a data point in plot) of four software
 #output one xyplot both x and y in log10 scale
function sys_xyplot(dpt,nsamp)

		name="simu-$dpt-$nsamp"
		df_true = read_true(name)
		true_maf=df_true.true_maf

		# name="$name-$(sigs[s])"
		df_snp=read_ngsPool_snp2(name)
		df_snape_p, df_snape_nallele = read_snape(name)
		df_p2 =read_p2(name)
		df_varscan=read_varscan(name)[:,[1,2,3,5]] #4,6 as FET p value
		if isempty(df_snape_p) #snape doesn't cope with high amount of individuals well--out of bounds errors
			all_sites=outerjoin(df_true, df_snp,df_p2,df_varscan,on=[:chrom, :position])
			# all_sites = all_sites[:,3:end]
			#reorder for plotting
			all_sites=all_sites[:,[:true_maf, :popoo, :GSS_ngsPool, :cross_samp_maf, :each_samp_maf]]
			#rename for legend
			rename!(all_sites, [:true_MAF, :Popoolation, :GSS_ngsPool, :cross_sample_VarScan, :each_sample_VarScan])
			all_sites2 = stack(all_sites, [:Popoolation, :GSS_ngsPool, :cross_sample_VarScan, :each_sample_VarScan],
			[:true_MAF])
		else
			all_sites=outerjoin(df_true, df_snp,df_p2, df_snape_p, df_varscan,on=[:chrom, :position])
			# all_sites = all_sites[:,3:end]
			#reorder for plotting
			all_sites=all_sites[:,[:true_maf, :popoo, :GSS_ngsPool, :maf_snape, :cross_samp_maf, :each_samp_maf]]
			#rename for legend
			rename!(all_sites, [:true_MAF, :Popoolation, :GSS_ngsPool, :Snape, :cross_sample_VarScan, :each_sample_VarScan])
			all_sites2 = stack(all_sites, [:Popoolation, :GSS_ngsPool, :Snape, :cross_sample_VarScan, :each_sample_VarScan],
			[:true_MAF])
		end

		xyplot = @df all_sites2 StatsPlots.scatter(
			:true_MAF,
			:value,
			group = :variable, dpi=300,  #background_color=:transparent, foreground_color=:black,
			# title = "",
			xlabel = "True MAF",
			ylabel = "Estimated MAF",#leftmargin=10mm,
			m = (0.5, [:cross :hex :star7 :square :circle], 3), grid=false,
			# bg = RGB(0.2, 0.2, 0.2), #black background
			fg_legend = :transparent, background_color_legend=:transparent,
			xscale=:log10, size=(600, 500),
			# legend=:topleft,  ylim=(1e-3,0.3),
			yscale=:log10, ylim=(1e-3,0.5),
			legend=:bottomleft
			# legend=:topleft
			)
			# Plots.plot!(yformatter = y->string(round.(10 .^ (-3:0.25:-0.5), digits=3))) #not working
			# plot!(margin=15mm)
			plot!(true_maf, true_maf, xscale=:log10, color=:grey,label="true MAF")
			vline!([0.05], line=:dot,label="", linewidth=3,alpha=0.6)
		display(xyplot)
		Plots.savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\1000sites_5programs_dpt($dpt)_samp$nsamp.png")
		# return xyplot
end


# #3d plot (not the best display way)
# xyplot=Plots.scatter()
# for i in 1:size(est2, 2) #column number
# 	Plots.plot!(true_maf, fill(i,size(est2, 1)), est2[:,i]) #x, z(i here as the profile value), y
# 		# Plots.plot!(true_maf, fill(i,size(est2, 1)), est2[:,i], fill=(0, .5,:orange)) #x, z(i here as the profile value), y
# 		# Plots.scatter!(true_maf, fill(i,size(est2, 1)), est2[:,i], leg=true)
# 	Plots.plot!(true_maf, fill(i,size(est2, 1)), true_maf, alpha=0.5) #x, z(i here as the profile value), y
#
# 	# surface!(true_maf, est2[:,i], fill(i,size(est2, 1))) #no it connect the points
# end


####function related with replicates relate with the script "script_qq"#########
	#fixed true MAF and 1000 replicates

#get the mean from 1000 replicates
function mean_replicates(dpt, nsamp)
    #all the replicates => a dataframe containing only the mean maf and variant called number

    qq_all=[0; 10 .^ range(-4,-1,length=20)] #all population maf for visualisation

    #average mean maf out of 1000 replicates at 21 MAF points
    all_maf=DataFrame(true_maf=qq_all, ngspool=zeros(21), p2=zeros(21),
        crs_varscan=zeros(21), each_varscan=zeros(21),
        snape=zeros(21))
    a = Array{Union{Missing, Float64},2}(undef,21, 6) #uninitialized array
    a .= missing; a[:,1]=qq_all
    # println(a)
    all_maf=DataFrame(a, [:true_maf, :ngspool,:snape,:p2,:crs_varscan,:each_varscan])
    #num of variant called out of 1000 replicates
    all_called=similar(all_maf); all_called .= 0; all_called.true_maf = qq_all


    i=1
    for qq in qq_all
        if qq == 0
            qq="0"
        end

        name="simu-$dpt-$nsamp-$qq"
        # df_true = read_true(name) #[:chrom, :position, :true_maf]
        df_ngs=dropmissing(read_ngsPool_snp2(name)) #[:chrom, :position, :GSS_ngsPool]

        df_p2 = dropmissing(read_p2(name))  #[:chrom, :position, :popoo]

        df_varscan = dropmissing(read_varscan(name), [:cross_samp_p, :each_samp_maf]) #[:chrom,:position,:cross_samp_maf, :cross_samp_p, :each_samp_maf, :each_samp_p])

        df_snape_p=read_snape(name)[1]
        if !isempty(df_snape_p)
            df_snape = dropmissing(df_snape_p, [:maf_snape]) #[:chrom,:position,:maf_snape]
            all_maf.snape[i] = mean(df_snape.maf_snape)
            all_called.snape[i] = length(df_snape.maf_snape)
        end

        all_maf.ngspool[i] = mean(df_ngs.GSS_ngsPool)
        all_called.ngspool[i] = length(df_ngs.GSS_ngsPool)

        all_maf.p2[i] = mean(df_p2.popoo)
        all_called.p2[i] = length(df_p2.popoo)

        all_maf.crs_varscan[i] = mean(df_varscan.cross_samp_maf)
        all_called.crs_varscan[i] = length(df_varscan.cross_samp_maf)

        all_maf.each_varscan[i] = mean(df_varscan.each_samp_maf)
        all_called.each_varscan[i] = length(df_varscan.each_samp_maf)

        i+=1
    end
    # replace!() #NaN and missing are the same

    all_metrics= all_called[:,2:end] ./ 1000
    all_metrics.true_maf =all_called.true_maf
	#first row true MAF=0, all else are not

    # all_FPR = all_metrics[1, :]
    # all_power = all_metrics[2:end, :]
    return all_maf, all_metrics
end

######### boplot: x as true maf, y as rmse or mpe, series as software
#input one depth and one sample size => boxplot comparing four software
function RMSE_MPE_replicates(dpt, nsamp)
	rmse_df = DataFrame(program="", qq=0.0, rmse=0.0)

    qq_all=10 .^ range(-4,-1,length=20) #all population maf for visualisation

    i=1
    for qq in qq_all
		if occursin("0.", name)
			name="simu-$dpt-$nsamp-$qq"
		else
			name="simu-$dpt-$nsamp"
		end
		df_p2 = read_p2(name)  #[:chrom, :position, :popoo]

        name="simu-$dpt-$nsamp-$qq"
        # df_true = read_true(name) #[:chrom, :position, :true_maf]
        df_ngs=read_ngsPool_snp2(name) #[:chrom, :position, :GSS_ngsPool]


        df_varscan = read_varscan(name) #[:chrom,:position,:cross_samp_maf, :cross_samp_p, :each_samp_maf, :each_samp_p])

        df_snape_p=read_snape(name)[1] #[:chrom,:position,:maf_snape]

        j=2#start from ngsPool
        for df_i in [df_ngs, df_snape_p, df_p2, df_varscan]
            if !isempty(df_i)
                n_row = nrow(df_i)
                true_maf=fill(qq, n_row)
                maf = df_i[:,end]
                rmse_vec = RMSE_missing(true_maf, maf)[1][:]
                mpe_vec = MPE_missing(true_maf, maf)[1][:]
				n_row = length(rmse_vec)
				name=names(df_i)[end]
				if df_i == df_varscan
					name="each-sample_varscan"
				end
				rmse_df_i = DataFrame(program=fill(name, n_row),
					qq=fill(qq, n_row), rmse=rmse_vec)

				# println(names(df_i)[end])

				rmse_df = vcat(rmse_df,rmse_df_i)
            end
            j+=1
        end #df_i

		# df_varscan
		n_row = nrow(df_varscan)
		true_maf=fill(qq, n_row)
		maf = df_varscan[:,end-2]
		rmse_vec = RMSE_missing(true_maf, maf)[1][:]
		mpe_vec = MPE_missing(true_maf, maf)[1][:]
		n_row = length(rmse_vec)
		# println(names(df_varscan)[end-2])
		rmse_df_i = DataFrame(program=fill("cross-sample_varscan", n_row),
			qq=fill(qq, n_row), rmse=rmse_vec)
		rmse_df = vcat(rmse_df,rmse_df_i)

        i+=1
    end #qq
	rmse_df = rmse_df[2:end, :] #delete the first row (meaningless)

	println(head(rmse_df))
	dropmissing!(rmse_df, [:rmse])
	println("after dropping, NAs?")
	a=findall(isnan, rmse_df.rmse)
	println(a)
	if !isempty(a)
		println(nrow(rmse_df),ncol(rmse_df))
		rmse_df=rmse_df[setdiff(1:end, a), :]
	end


	# println("rmse_df  ", findall(x->x<0, completecases(rmse_df)))
	bp = groupedboxplot(rmse_df.qq, rmse_df.rmse, group = rmse_df.program, #xscale=:log10,
		bar_width = 0.01, xlabel="True MAF", ylabel="Root mean square error (RMSE)",
		fg_legend = :transparent, background_color_legend=:transparent,label="",
		#legend=:topright, #		fg_legend = :transparent,
 		legendfontsize=9,grid=false) #, xscale=:log10)
	display(bp)
	# savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\qq-rmse-nsamp$nsamp-dpt$dpt.png")

    return bp #boxplot
end

############################## script: power_FPR ###############################
#power and fasle positive rate


#return the matrix of the proportion of called sites
	#outer matrix: the number of depth as row number, the number of sample size as column number
	#each element of the outer matrix => inner matrix
	#inner matrix: row as significance level, column as corresponding software
#depending on the input MAF probability (if all 0 => not polymorphic site => the output is FPR)
	#(if not => the matrix of power
	# dpt-d
	# nsamp-n
	# significance-s
	# program-p
function metric_2mtx(dpts, nsamps, sigs)
    n_dpts=length(dpts)
    n_nsamps=length(nsamps)
    n_sigs=length(sigs)

    power = Array{Array{Float64,2},2}(undef,n_dpts, n_nsamps)
    for d in 1:n_dpts, n in 1:n_nsamps
        power_i = zeros(n_sigs, 4) #row as 100 significance, col as 4programs(varscan occupy two)
        for s in 1:n_sigs
            dpt=dpts[d]
            nsamp=nsamps[n]
            sig=sigs[s]

            # dpt=1; nsamp=30; sig=0.001
			if occursin("0.", name)
				name="simu-$dpt-$nsamp-$sig"
			else
				name="simu-$dpt-$nsamp"
			end

            df_p2 = read_p2(name)
			name="simu-$dpt-$nsamp-$sig"
            df_ngspool = read_ngsPool_snp2(name)

            df_p2 = filter(:popoo => x -> !iszero(x), df_p2)

            df_snape = read_snape(name)[1]
            # df_snape = filter(:maf_snape => x -> !iszero(x), df_snape)

            df_varscan = read_varscan(name)[:,[1,2,3,5]]
            #cross-sample and each-sample call the sites samely
            # findall(ismissing, convert(Matrix, df_varscan)) #0
            # df_varscan = filter(:cross_samp_maf => x -> !iszero(x), df_varscan)
            # df_varscan = filter(:each_samp_maf => x -> !iszero(x), df_varscan)

            #checked: all sites are (polymorphic) called variant sites

            dfs=[df_ngspool, df_snape, df_p2,  df_varscan]
            for p in 1:4
                power_i[s, p] = nrow(dfs[p])
            end
        end#significance
        power[d, n] = power_i
    end
    power .= power ./ 1000 # proportion of sites called out of 1000 sites
    return power
end

##simulation linux script: power_FPR_wider
	#only difference from the function metric_2mtx() above is the naming of Popoolation2 output (*.sync)
	 #significance level not annotated in *.sync as Popoolaiton2 does not require significance level
function metric_2mtx2(dpts, nsamps, sigs)
    n_dpts=length(dpts)
    n_nsamps=length(nsamps)
    n_sigs=length(sigs)

    power = Array{Array{Float64,2},2}(undef,n_dpts, n_nsamps)
    for d in 1:n_dpts, n in 1:n_nsamps
        power_i = zeros(n_sigs, 4) #row as 100 significance, col as 4programs(varscan occupy two)
        dpt=dpts[d]
        nsamp=nsamps[n]
        for s in 1:n_sigs
            sig=sigs[s]

            # dpt=1; nsamp=30; sig=0.001

			name="simu-$dpt-$nsamp"
			# name="simu-$dpt-$nsamp-$sig"

			df_p2 = read_p2(name)
			df_p2 = filter(:popoo => x -> !iszero(x), df_p2)

			name="simu-$dpt-$nsamp-$sig"
			df_ngspool = read_ngsPool_snp2(name)

            df_snape = read_snape(name)[1]
            # df_snape = filter(:maf_snape => x -> !iszero(x), df_snape)

            df_varscan = read_varscan(name)[:,[1,2,3,5]]
                #cross-sample and each-sample call the sites samely
            # findall(ismissing, convert(Matrix, df_varscan)) #0
            # df_varscan = filter(:cross_samp_maf => x -> !iszero(x), df_varscan)
            # df_varscan = filter(:each_samp_maf => x -> !iszero(x), df_varscan)
			cro=count(x->!ismissing(x) && !iszero(x), df_varscan[:,"cross_samp_maf"])
			each=count(x->!ismissing(x) && !iszero(x), df_varscan[:,"cross_samp_maf"])
			if cro != each
				println("not equal")
				println("cross $cro")
				println("each $each")
			end
            #checked: all sites are (polymorphic) called variant sites

            dfs=[df_ngspool, df_snape, df_p2,  df_varscan]
            # println(nrow.(dfs))

            for p in 1:4 #4 programs
                power_i[s, p] = nrow(dfs[p])
            end
        end#significance
        power[d, n] = power_i
    end
    power .= power ./ 1000 # proportion of sites called out of 1000 sites
    return power
end #for files in power_wider (sig 0.1:0.99)

##linux script: power_FPR_wider
	#input: one depth and one sample size, 100 significance levels, 4 programs
	#output: the matrix of the proporiton of called variants (it depends on the true MAF to have this matrix as power or false positive rate)
function metric_mtx3(dpt, nsamp, sigs, replicates)
    # n_dpts=length(dpts)
    # n_nsamps=length(nsamps)
	n_sigs=length(sigs)
    n_reps=length(replicates)

    power = Array{Array{Float64,2},1}(undef,n_sigs)
    for s in 1:n_sigs
        power_i = zeros(n_reps, 4) #row as 100 significance, col as 4programs(varscan occupy two)
        # dpt=dpts[d]
        # nsamp=nsamps[n]
        for rep in 1:n_reps
            sig=sigs[s]

            # dpt=1; nsamp=30; sig=0.001

            name="simu-$dpt-$nsamp-$rep"
            df_p2 = read_p2(name)
            df_p2 = filter(:popoo => x -> !iszero(x), df_p2)
            name="simu-$dpt-$nsamp-$rep-$sig"
            # println(name)
            df_ngspool = read_ngsPool_snp2(name)

            df_snape = read_snape(name)[1]
            # df_snape = filter(:maf_snape => x -> !iszero(x), df_snape)

            df_varscan = read_varscan(name)[:,[1,2,3,5]]
                #cross-sample and each-sample call the sites samely
            # findall(ismissing, convert(Matrix, df_varscan)) #0
            # df_varscan = filter(:cross_samp_maf => x -> !iszero(x), df_varscan)
            # df_varscan = filter(:each_samp_maf => x -> !iszero(x), df_varscan)

            #checked: all sites are (polymorphic) called variant sites

            dfs=[df_ngspool, df_snape, df_p2,  df_varscan]
            # println(nrow.(dfs))

            for p in 1:4 #4 programs
                power_i[rep, p] = nrow(dfs[p])
            end
        end#significance
        power[s] = power_i
    end
    power .= power ./ 1000 # proportion of sites called out of 1000 sites
	# return power

		power_mtx=zeros(0,4)
		for s in 1:n_sigs
			power_mtx= vcat(power_mtx, power[s])
		end
		return power_mtx #what matter is the corresponding of power and FPR, all else not

end #for files in power_wider (sig 0.1:0.99)
