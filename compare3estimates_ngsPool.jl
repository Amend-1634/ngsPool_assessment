#func_read_plot2(2020/8/12)
using CSV, CodecZlib, Plots, DelimitedFiles, Plots.PlotMeasures, DataFrames

#form: :chrom(string), :position(Int), :algorithm_software(float)]
#name: the full name

#read the gzipped file containing estimation of ngsPool with SNP calling
function read_ngsPool_snp(filename) #all are SNP calling
    maf_snp = open(GzipDecompressorStream, filename, "r") do stream #estimated
        CSV.read(stream) #output to df
    end #1000 rows--sites, maf as dataframe

    # df_snp = maf_snp[!, [:chrom,:position, :maf,:freqMax, :freqE]]

    df_snp = hcat(maf_snp[!,1], maf_snp[!, [:position, :maf,:freqMax, :freqE]])

    rename!(df_snp, [:chrom, :position, :GSS_ngsPool,:freqMax_ngsPool, :freqE_ngsPool])
     #Golden-section search (GSS); maf with maximum SFS likelihood, expected maf with SFS likelihood
    return df_snp #dataframe
end

cd("C:\\simu_raw\\sig")
filename="simu-1-100-0.1.snp.gz"
df_snp=read_ngsPool_snp(filename)

#read the txt file containing true MAFs
function read_true(filename)
    if occursin("txt", filename)
        daf = readdlm(filename, '\t', String, '\n') #true #daf as matrix
        daf_1 = daf[:, 1]
        daf_2 = parse.(Int64, daf[:, 2])
        daf_end = parse.(Float64, daf[:, 5])#wrong, true maf is in 5th column, the end column is the simple reads proportion of reads(no meaning)
        df_true = DataFrame(chrom = daf_1, position = daf_2, true_maf = daf_end)
    elseif occursin("gz", filename)
        daf = open(GzipDecompressorStream, filename, "r") do stream #estimated
            CSV.read(stream, header=false) #output to df
        end #1000 rows--sites, maf as dataframe
        df_true = daf[!,[1,2,5]]
        df_true = rename!(df_true, [:chrom, :position, :true_maf])
    end
    return df_true
end

cd("C:\\simu_raw\\sig")
filename="test-1-100.txt"
df_true = read_true(filename)

df_snp2=innerjoin(df_snp, df_true, on=[:chrom, :position])
df_snp3=sort(df_snp2, :true_maf)


mtx_snp3=convert(Matrix,df_snp3[3:end])
scatter(mtx_snp3[:,end], mtx_snp3[:,1:end-1])
plot(mtx_snp3)
#differ but very close

#read the estimation (.sync file) of Popoolation2
function read_p2(filename) #Popoolation2
    daf = readdlm(filename, '\t', String, '\n')
    daf_count = daf[:, 4]
    nsite = length(daf_count)
    mafs = zeros(nsite)

    for i in 1:length(daf_count)
        count = parse.(Int, split(daf_count[i], ":"))
        maf = sort(count)[end-1]/sum(count) #minor allele frequency
        mafs[i] = maf
    end
    daf_2 = parse.(Int, daf[:,2])
    df_final =
    return df_final
end
##################################
#using standardized funcions to read
include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_read_standard.jl")
    include("H:\\julia\\ngsJulia\\ngsPool\\compare\\func_analysis_standard.jl")
    cd("C:\\simu_raw\\power\\power_rnd") #power
    name="simu-1-50-0.005" #pick a random one
    df_snp = read_ngsPool_snp(name)
    name="simu-1-50" #pick a random one
    df_true = read_true(name)
    df_join=innerjoin(df_true, df_snp, on = [:chrom, :position])
names(df_join)

y=convert(Matrix, df_join[:, 4:6])
y2=map(x -> x > 0.5 ? 1-x : x, y)
x=map(x -> x > 0.5 ? 1-x : x, df_join.true_maf)

findall(x->x>0.5, y2)

#root mean square error
RMSE_mtx = let RMSE_mtx=zeros(nrow(df_join),3)
    for i in 1:nrow(df_join), j in 1:3
        RMSE_mtx[i,j]=RMSE.(df_join.true_maf[i], df_join[i, 3+j])
    end
    RMSE_mtx
end

#mean percentage error
MPE_mtx = let MPE_mtx=zeros(nrow(df_join),3)
    for i in 1:nrow(df_join), j in 1:3
        MPE_mtx[i,j]=MPE.(df_join.true_maf[i], df_join[i, 3+j])
    end
    MPE_mtx
end

RMSE_mtx
MPE_mtx

#four plots comparing the three estimations (GSS, maximum likelihood, expected estimation) of ngsPool 
using Plots.PlotMeasures
p=plot(dpi=300, layout=(2,2), size=(1200, 1000), grid=false, margin=9mm, guidefontsize=13,  tickfontsize=12, legendfontsize=12,
    fg_legend = :transparent, background_color_legend=:transparent)
    boxplot!(RMSE_mtx,label="",xtick=(1:3, ["GSS", "MLE", "Expected"]), ylabel="Root mean square error", subplot=3, title="(c)")
    boxplot!(MPE_mtx,label="",xtick=(4:7, ["GSS", "MLE", "Expected"]), ylabel="Mean percentage error", subplot=4, title="(d)")

    scatter!(x, y2, dpi=300, grid=false, legend=:topleft,xlabel="true MAF",alpha=0.7,ylabel="Estimated MAF",
        label=["Goden section search" "Maximum likelihood" "Expected"], subplot=1, markersize=3, title="(a)")
    scatter!(y2[:,1], y2[:,2:3], markersize=3, alpha=0.3, xlabel="GSS MAF estimate",ylabel="Estimated MAF", title="(b)",
        label=["Maximum likelihood" "Expected"], legend=:topleft, subplot=2)

savefig("H:\\julia\\ngsJulia\\ngsPool\\compare\\plot\\3estimates_ngsPool.png")

# names(df_join)
