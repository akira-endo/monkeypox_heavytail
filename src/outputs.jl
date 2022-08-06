using CSV
using Distributions
using DataFrames
using LinearAlgebra
using Measures
using Plots
using Random
using SpecialFunctions
using StatsPlots
using StatsBase
using ProgressBars
include("heavytail.jl")
include("outbreaklikelihood.jl")
gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=10,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))
    
# Util function
Base.:-(x::Tuple{Float64,Float64,Float64},y::Int)=x.-y

# Parameters
## Weibull parameters
if !isdefined(@__MODULE__,:κ_msm) κ_msm=0.6 end # Pareto exponent
if !isdefined(@__MODULE__,:α_msm) α_msm=0.15 end
#if !isdefined(@__MODULE__,:κ_other) κ_other=1.66 end # Pareto exponent
#if !isdefined(@__MODULE__,:κ_other) α_other = 0.04 end
## Epi parameters
if !isdefined(@__MODULE__,:inf_period) inf_period=21 end
## Demography parameters
if !isdefined(@__MODULE__,:p_msm) p_msm = 0.02 end# proportion of sexually active MSM
#if !isdefined(@__MODULE__,:p_other) p_other = 0.86 end # proportion of sexually active non-MSM
if !isdefined(@__MODULE__,:p_hem) p_hem = 0.86 end # proportion of sexually active heterosexual men
if !isdefined(@__MODULE__,:p_hew) p_hew = 0.89 end # proportion of sexually active heterosexual women

## 2022 outbreak parameters
#currentcases = 584 - 16 # as of 29 May 2022. total - discarded from Global.health
if !isdefined(@__MODULE__,:currentcases) currentcases = 776 - 48 end # as of 31 May 2022. total - discarded from Global.health

## crontrol variables
if !isdefined(@__MODULE__,:printtable) printtable = true end
if !isdefined(@__MODULE__,:outbreak_iters) outbreak_iters = 100000 end
if !isdefined(@__MODULE__,:hist_upp) hist_upp = :auto end

# Degree distribution over 1 year
deg1y_msm = truncated(Weibull2(α_msm, κ_msm), lower=1)
#deg1y_other = truncated(Weibull2(α_other, κ_other), lower=1)
deg1y_hem = truncated(Weibull2(α_hem, κ_hem), lower=1)
deg1y_hew = truncated(Weibull2(α_hew, κ_hew), lower=1)

# Preparation for draws from ~ x*pdf(x)
ED_3w_msm = truncated(Weibull(params(deg1y_msm)[1], params(deg1y_msm)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull 
ED_3w_hem = truncated(Weibull(params(deg1y_hem)[1], params(deg1y_hem)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull
ED_3w_hew = truncated(Weibull(params(deg1y_hew)[1], params(deg1y_hew)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull

# Pre-draws from excess degree +1 distribution for efficiency
Random.seed!(2022)
drawsize = 1000000
edeg3w_msm = rand(ED_3w_msm, drawsize)
#edeg3w_other = rand(ED_3w_other, drawsize)
edeg3w_hem = rand(ED_3w_hem, drawsize)
edeg3w_hew = rand(ED_3w_hew, drawsize)

# Degree distributions for different types
deg3w_msm = truncated(Weibull(params(deg1y_msm)[1], params(deg1y_msm)[2]*inf_period/365), lower=inf_period/365)
#deg3w_other = truncated(Weibull(params(deg1y_other)[1], params(deg1y_other)[2]*inf_period/365), lower=inf_period/365)
deg3w_hem = truncated(Weibull(params(deg1y_hem)[1], params(deg1y_hem)[2]*inf_period/365), lower=inf_period/365)
deg3w_hew = truncated(Weibull(params(deg1y_hew)[1], params(deg1y_hew)[2]*inf_period/365), lower=inf_period/365)

contDirac = TruncatedNormal(0,1e-10,0,0.1)
deg3w_msmavg = MixtureModel([contDirac,deg3w_msm],[1-p_msm ,p_msm])
deg3w_hemavg = MixtureModel([contDirac,deg3w_hem],[1-p_hem ,p_hem])
deg3w_hewavg = MixtureModel([contDirac,deg3w_hew],[1-p_hew ,p_hew])


### Outputs
if !isdefined(@__MODULE__,:R0only) R0only = false end
if !isdefined(@__MODULE__,:R0finite) R0finite = false end
if !isdefined(@__MODULE__,:R0deplete) R0deplete = false end
if !isdefined(@__MODULE__,:R0assortative) R0assortative = false end
if R0only
    if any(R0finite)
        if any(.!R0finite) # generate new: takes ~ 1 hour
            ssrange=round.(Int,10 .^((10:50)./10))
            @time R0finite=[quantile([mean(max.(0,ISR_Weibull_finite(ED_3w_msm,ss).-1)) for i in 1:1000],[0.025,0.5,0.975]) for ss in ssrange]
            
            Rfvalues=reduce(hcat,R0finite)|>permutedims
            R0finite_output=[DataFrame([ssrange],[:primarysize]) DataFrame(Rfvalues,string.([0.025,0.5,0.975]))]
        else # load existing data
            R0finite_df=CSV.read("output/R0finite_quantile.csv",DataFrame)
            ssrange=R0finite_df.primarysize
            Rfvalues=R0finite_df[:,2:end]|>Matrix
        end

        plot(ssrange,Rfvalues[:,2].*[0.1 0.5],ribbons=((Rfvalues[:,2].-Rfvalues[:,1]).*[0.1 0.5],(Rfvalues[:,3].-Rfvalues[:,2]).*[0.1 0.5]),yaxis=:log,xaxis=:log,xticks=10 .^(1:5),
yticks=10.0 .^ (-1:3),label="SAR: ".*string.([0.1 0.5]),ylim=(0.1,1000),ylabel="R₀ over MSM sexual network",xlabel="size of primary cases")
hline!([1],linestyle=:dot,color=:black,label="")|>display

    elseif any(R0deplete)
        if !isdefined(@__MODULE__,:popsize) popsize = 500000 end
        if !isdefined(@__MODULE__,:lowlim) lowlim = "" end
        
        if any(.!R0deplete) # generate new: takes ~ 1 hour
            
            #=
            pop=rand(deg3w_msm,popsize)
            @time samples=[sample(pop,Weights(pop),min(100010,popsize),replace=false) for i in 1:1000]
            pop=nothing
            samplemat = reduce(hcat, samples)                
            R0mat=max.(0,samplemat.-1)
            sizerange=11:min(100000,popsize-10)
            movavgs=[col[i-10:i+10]|>mean for i in sizerange, col in eachcol(R0mat)]
            qs=quantile.(eachrow(max.(1e-8,movavgs)), Ref([0.025,0.5,0.975]))
            Rdvalues=reduce(hcat,qs)|>permutedims
            =# 
            rng=MersenneTwister(popsize)
            if lowlim==""
                deg=rand(deg3w_msm,popsize) # degree of susceptible pop
            else
                deg=rand(truncated(deg3w_msm.untruncated,quantile(deg3w_msm,lowlim),Inf),popsize)
            end
            pop=copy(deg)
            Rmat=Matrix{eltype(deg)}(undef, popsize,1000) # initialize result vector
            sumdeg=sum(deg)
            @time begin
            for j in 1:1000
                pop.=deg
                susprop=1.0
                for i in 1:min(100010,popsize)
                    sampleid=sample(rng,1:popsize, Weights(deg))
                    susprop-=pop[sampleid]/sumdeg
                    Rmat[i,j]=max(0,pop[sampleid]-1)*susprop/(1-pop[sampleid]/sumdeg)
                    pop[sampleid]=0
                end
            end
            end
            sizerange=11:min(100000,popsize-10)
            movavgs=[col[i-10:i+10]|>mean for i in sizerange, col in eachcol(Rmat)]
            qs=quantile.(eachrow(max.(1e-8,movavgs)), Ref([0.025,0.5,0.975]))
            Rdvalues=reduce(hcat,qs)|>permutedims    
 
            R0deplete_output=[DataFrame([sizerange],[:outbreaksize]) DataFrame(Rdvalues,string.([0.025,0.5,0.975]))]
            pop=nothing;deg=nothing;qs=nothing
        else # load existing data
            verstring=string(lowlim)*"_"*string(popsize÷1000)
            R0deplete_df=CSV.read("output/R0sdeplete"*verstring*"k.csv",DataFrame)
            sizerange=R0deplete_df.outbreaksize
            Rdvalues=R0deplete_df[:,2:end]|>Matrix
        end
        qlines=[Rdvalues[:,2] Rdvalues[:,2].-Rdvalues[:,1] Rdvalues[:,3].-Rdvalues[:,2]]
                                
        plot(sizerange,qlines[:,1].*[0.1 0.5],ribbon=(qlines[:,2].*[0.1 0.5],qlines[:,3].*[0.1 0.5]),
    xaxis=:log,yaxis=:log,xlim=(10,100000),ylim=(0.01,100),xticks=10 .^(1:5),label="SAR: ".*string.([0.1 0.5]),
    xlabel="cumulative cases", ylabel="effective R over sexual network")
hline!([1],linestyle=:dot,color=:black,label="")|>display

    else
        if R0assortative                                                    
            R0HS=[√(R0(truncated(deg3w_hem.untruncated,quantile(deg3w_hem,lowlim),Inf))*
                R0(truncated(deg3w_hew.untruncated,quantile(deg3w_hew,lowlim),Inf))) for lowlim in 0:0.001:0.999]
            plot(0:0.001:0.999,R0HS.*[0.1 0.5],
    yaxis=:log,ylim=(0.1,10),label="SAR: ".*string.([0.1 0.5]),legend=:topleft,
    xlabel="threshold for core group (quantile)", ylabel="R₀ among HS core group")
    hline!([1],linestyle=:dot,color=:black,label="")|>display
        else
            sarrange=exp.([-10:0.1:-0.5;-0.4:0.01:0])
            R0values=[R0(deg3w_msm) √(R0(deg3w_hem)*R0(deg3w_hew))].*sarrange;
        end
    end
else


# Degree over 1 year
Random.seed!(2022)
D_msm=[floor.(Int,rand_TWeibull(deg1y_msm,10000)) for i in 1:100]
D_hem=[floor.(Int,rand_TWeibull(deg1y_hem,10000)) for i in 1:100]
D_hew=[floor.(Int,rand_TWeibull(deg1y_hew,10000)) for i in 1:100]
freqs=[sum.([(==).(1:29);>(29)],Ref(x))./10000 for x in (D_msm[1], D_hem[1], D_hew[1])]

labels=["MSM" "heterosexual men" "heterosexual women"]
groupedbar(reduce(hcat,freqs),label=labels, linealpha=0, 
    xlabel = "number of partners over 1 year", ylabel="relative frequency")|>display
    
loglogplot=plot([sort.(D_msm),sort.(D_hem),sort.(D_hew)], [(length(D_msm[1]):-1:1)./length(D_msm[1]), (length(D_hem[1]):-1:1)./length(D_hem[1]), (length(D_hew[1]):-1:1)./length(D_hew[1])],xaxis=:log,yaxis=:log,label="",
    xlabel = "number of partners over 1 year", ylabel="relative frequency (cumulative)", color=repeat([1 2 3],inner=(1,100)),linealpha=0.2)
if isdefined(@__MODULE__,:partnersdata)
    freqs=[sum.([(==).(1:29);>(29)],Ref(vcat(fill.(1:length(x),x)...)))./sum(x) for x in eachcol(partnersdata)]
    groupedbar(reduce(hcat,freqs),label=labels, linealpha=0, 
    xlabel = "number of partners over 1 year", ylabel="relative frequency")|>display
    plot(reduce(hcat,freqs),label=labels,markersize=3, 
    xlabel = "number of partners over 1 year", ylabel="relative frequency")|>display


    y=[begin
            nonzero=col.>0
            1 .-((cumsum(col[nonzero]))[1:end-1]./sum(col))
        end for col in eachcol(partnersdata)]
    x=[(1:length(col))[col.>0][2:end] for col in eachcol(partnersdata)]
    plot!(loglogplot, x,y,linealpha=0,markershape=:circ,markersize=3,color=(1:3)',markerstrokewidth=1,label=labels)
end
display(loglogplot)


# Likelihood of an outbreak of current size or greater
if printtable
for SAR in [0.05, 0.1, 0.2, 0.5]
    outbreaklikelihood1(edeg3w_msm, edeg3w_hem, edeg3w_hew, currentcases, SAR, iter=outbreak_iters)
    print("\n\n")
    outbreaklikelihood2(edeg3w_msm, edeg3w_hem, edeg3w_hew, deg3w_msm, deg3w_msmavg, deg3w_hemavg,deg3w_hewavg, currentcases, SAR, iter=outbreak_iters)
    print("\n\n\n")
end
end

## Main plots

# The likelihood of outbreak vs SAR, sexually-associated MSM
Random.seed!(2022)
sarrange = [0;0.005;0.01:0.01:0.05;0.075;0.1:0.05:1]
seedrange = 1:5
currobsize=[simulatetillsize(edeg3w_msm,currentcases,sar;seed=s,iter=outbreak_iters) for sar in sarrange, s in seedrange]

plot(sarrange,currobsize,xlim=(0,1),
label=string.(seedrange').*[" initial case" fill("",1,4)],legend=:bottomright, ylim=(0,1),
xlabel="SAR per sexually-associated contact", ylabel="risk of outbreak (≥ $currentcases cases)", title="MSM index cases: sexually-associated exposure",
linewidth=1.5,color=(5:-1:1)'.*40,palette=:romaO) |> display

# The likelihood of outbreak vs SAR, random MSM
Random.seed!(2022)
sarrange = [0.02:0.02:0.08;0.1:0.1:1]
seedrange = 1:5
currobsize_rMSM=[simulatetillsize(edeg3w_msm,currentcases,sar;seed=s, init_deg=deg3w_msm,iter=outbreak_iters) for sar in sarrange, s in seedrange]

plot(sarrange,currobsize_rMSM,xlim=(0,1),
label=string.(seedrange').*[" initial case" fill("",1,4)],legend=(0.12,0.9), ylim=(0,1),
xlabel="SAR per sexually-associated contact", ylabel="risk of outbreak (≥ $currentcases cases)",title="MSM index cases: other exposure",
linewidth=1.5,color=(5:-1:1)'.*40,palette=:romaO) |> display

# The likelihood of outbreak vs SAR, random gen pop

Random.seed!(2022)
sarrange = [0.01:0.01:0.04;0.05;0.1:0.1:1]
seedrange = [5,10,50,100,500]
currobsize_ravgMSM=[simulatetillsize(edeg3w_msm,currentcases+s, sar;seed=s, init_deg=deg3w_msmavg,iter=outbreak_iters) for sar in sarrange, s in seedrange]

plot(sarrange,currobsize_ravgMSM,xlim=(0,1),
label=string.(seedrange').*[" initial cases" fill("",1,4)],legend=(0.12,0.9), ylim=(0,1),
xlabel="SAR per sexually-associated contact", ylabel="risk of outbreak (≥ $currentcases cases)",title="General index cases: other exposure",
linewidth=1.5,color=(5:-1:1)'.*50,palette=:romaO) |> display

# Outbreak sizes in non-MSM
Random.seed!(2022)
major_threshold = 20000
initdegs = repeat([[nothing,nothing], [deg3w_hemavg, deg3w_hewavg]],inner=10)
sarrange = repeat(0.1:0.1:1,2)
majobrisk=[simulatetillsize(edeg3w_hem,edeg3w_hew,major_threshold,sar;seed=1000, init_deg1=initdeg[1],init_deg2=initdeg[2], total=true,iter=1000) for (initdeg, sar) in zip(initdegs, sarrange)]

majobquant = quantile.(getfield.(majobrisk,:totalrecord),Ref((0.025,0.5,0.975))).-1000
plot(0.1:0.1:1,reshape(getindex.(majobquant,2),10,2),ribbon=(reshape(getindex.(majobquant,2).-getindex.(majobquant,1),10,2),reshape(getindex.(majobquant,3).-getindex.(majobquant,2),10,2)) ,
markershape=:none,markerstrokewidth=0,color=[1 2],xlim=:auto, ylim=(10,:auto), yaxis=:log, yticks=10 .^(1:4),
label="1000".*[" initial cases: sexually-associated" " initial cases: other"],legend=:topleft,
xlabel="SAR per sexually-associated contact", ylabel="Outbreak size",title="Outbreak size in non-MSM sexual network",
linewidth=1.5,linealpha=0.8) 
majobquart = quantile.(getfield.(majobrisk,:totalrecord),Ref((0.25,0.5,0.75))).-1000
plot!(0.1:0.1:1,reshape(getindex.(majobquart,2),10,2),ribbon=(reshape(getindex.(majobquart,2).-getindex.(majobquart,1),10,2),reshape(getindex.(majobquart,3).-getindex.(majobquart,2),10,2)),
color=(1:2)',label="") |> display


# R0 vs SAR
sarrange=exp.([-10:0.1:-0.5;-0.4:0.01:0])
R0values=[R0(deg3w_msm) √(R0(deg3w_hem)*R0(deg3w_hew))].*sarrange
R0valuesplot=plot(sarrange,R0values,yaxis=:log,ylim=(0.01,100),linewidth=1.5,
xlab="SAR per sexually-associated contact", ylab="R₀ over sexual network",
label=["MSM" "non-MSM"],legend=:topleft)
hline!(R0valuesplot,[1],linestyle=:dot,color=:black,label="") |> display

R0valuesplot=plot(sarrange,max.(0,min.(1,1 .- 1 ./R0values)),ylim=(0,1),linewidth=1.5,
xlab="SAR per sexually-associated contact", ylab="relative reduction in R₀ for control",
label=["MSM" "non-MSM"],legend=:topright)|>display

# Degree distributions with tail control
Random.seed!(2022)
αrange=0.04:0.01:100
κrange=0.6
qs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=1),0.99) for α in αrange, κ in κrange]
findq(q,qs)=findmin(x->abs(x-q),qs)
αs=αrange[last.(findq.([15,10],Ref(qs)))]
nonzeroPoisson(tw::Truncated{<:Weibull},n::Integer)=truncated.(Poisson.(ISR(tw,n,x->1-exp(-x))),lower=1)
origdeg=rand.(nonzeroPoisson(deg3w_msm,100000))#floor.(Int,rand(truncated(deg3w_msm.untruncated,lower=1),100))
orig_1p=nonzeroPoissonquantile(deg3w_msm,0.99)#round(quantile(origdeg,0.99),digits=0)
roundorig_1p=round(orig_1p,digits=1)
                                                                degs=[rand.(nonzeroPoisson(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=inf_period/365),100000)) for α in αs, κ in κrange]
freqs=[sum.([(==).(1:29);>(29)],Ref(x))./100000 for x in [[origdeg];degs]]
plot(bar.(freqs, linealpha=0, 
     ylabel="",legend=:none)...,
    xticks=[:none :none :native],
    ylim=(-0.01+0.009,hist_upp),xlim=(-0.4+4.8,31),layout=(3,1),size=(400,300),yticks=0:0.01:1,
label=["baseline, 1%tile: $roundorig_1p" "1%tile: 10" "1%tile: 15"],legend=:topright,
xlabel = ["" "" "number of partners over 21 days"],bottom_margin=[-2mm -2mm 4mm]) |> display

# R0 vs 1%ile
Random.seed!(2022)
αrange=[0.04:0.005:0.09;0.1:0.005:10]#BigFloat.(10 .^(-3:0.1:0))#
sarrange=[0.05,0.1,0.2,0.5]
κrange=repeat([κ_msm],length(sarrange))
R0s=repeat(sarrange,inner=length(κrange)÷length(sarrange))'.*max.(1e-15,[R0(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365),lower=inf_period/365),logvalue=false) for α in αrange, κ in κrange])
xs= [nonzeroPoissonquantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=inf_period/365),0.99,outbreak_iters) for α in αrange, κ in κrange]
R0capplot=plot(xs,R0s,yaxis=:log,ylim=(0.01,100),xlim=(0,30), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="R₀ over MSM sexual network",
legend=:topleft, label="SAR: ".*string.(sarrange'),
palette=:vikO,color=(1:4)'.*50,linewidth=1.5)
hline!(R0capplot,[1],linestyle=:dot,color=:black,linealpha=1,label="")
vline!(R0capplot,[orig_1p],color=:green,label="",linealpha=0.5,linestyle=:dash) |> display

# Reduction in R0 for control
R0capplot=plot(xs,max.(0,min.(1,1 .- 1 ./R0s)),ylim=(0,1),xlim=(0,30), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="relative reduction in R₀ for control",
legend=:topleft, label="SAR: ".*string.(sarrange'),
palette=:vikO,color=(1:4)'.*50,linewidth=1.5)
vline!(R0capplot,[orig_1p],color=:green,label="",linealpha=0.5,linestyle=:dash)|>display

# Eigenvector analysis
R_msm=R0(deg3w_msm)#meanm1sampling(edeg3w_msm)
R_hem=R0(deg3w_hem)
R_hew=R0(deg3w_hew)
r_hem=mean(deg3w_hemavg)
r_hew=mean(deg3w_hewavg)#mean(rand(deg3w_otheravg,1000000))

Rlist=[R_msm,R_hem,R_hew,r_hem,r_hew]
relRrange=0.:0.1:0.5
prange=[0,0.01,0.05,0.1]
@time eigs=[eigenNGM(1,R,p_m2o,Rlist...) for R in relRrange, p_m2o in prange];
evectors=getindex.(eigs,:vector)
plot(relRrange,last.(evectors),
xlabel="reproduction number ratio",ylabel="non-sexually-associated case ratio",
label=["MSM to non-MSM mixing: " "" "" ""].*string.(prange'),legend=:topleft) |> display

## Sensitivity analysis of the degree distribution shape
# MSM
αrange=BigFloat.(10 .^(-3:0.05:0))#0.04:0.01:10
κrange=repeat([0.77,0.66,0.88],2)
R0s1=repeat([0.1,0.5],inner=length(κrange)÷2)'.*max.(1e-15,[R0(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365),lower=inf_period/365),logvalue=false) for α in αrange, κ in κrange])
#xs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=1),0.99,100) for α in αrange, κ in κrange]
xs= [nonzeroPoissonquantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=inf_period/365),0.99,outbreak_iters) for α in αrange, κ in κrange]

labels="κ: ".*(x->string(x[1])*"–"*string(x[2]))(κrange[2:3])
R0capplot1=plot(xs,R0s1,yaxis=:log,ylim=(0.01,1000),xlim=(0,50), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="R₀ over sexual network",
legend=:bottomright, label=[labels "" ""],
palette=:vikO,color=(repeat([1,4],inner=3))'.*50,linewidth=[1.5 0.5 0.5])
hline!(R0capplot1,[1],linestyle=:dot,linealpha=1,color=:black,label="")
vline!(R0capplot1,[orig_1p],
    color=:green,label="",linealpha=0.5,linestyle=:dash) |> display 

# non-MSM

αrange=BigFloat.(10 .^(-3:0.05:0))
κrange=repeat([1.68,1.63,1.73,2.14,2.07,2.21],2)
@time R0s2=repeat([0.1,0.5],inner=length(κrange)÷2)'.*max.(1e-15,[R0(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365),lower=inf_period/365)) for α in αrange, κ in κrange])
xs= [nonzeroPoissonquantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=inf_period/365),0.99,outbreak_iters) for α in αrange, κ in κrange]
labels="κ: ".*(x->string(x[1])*"–"*string(x[2])).([[κrange[2:3]] [κrange[5:6]]])
R0capplot2=plot(xs,R0s2,yaxis=:log,ylim=(0.001,10),xlim=(0,5), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="one-way HS reproduction number",
legend=:bottomright, label=[""*labels[1] "" "" ""*labels[2] "" ""],
palette=:vikO,color=([repeat([1,2],inner=3);repeat([4,3],inner=3)])'.*50,linewidth=[1.5 0.5 0.5])
hline!(R0capplot2,[1],linestyle=:dot,linealpha=1,color=:black,label="")
vline!(R0capplot2,[nonzeroPoissonquantile(deg3w_hem,0.99)],
    color=:green,label="",linealpha=0.5,linestyle=:dash)
vline!(R0capplot2,[nonzeroPoissonquantile(deg3w_hew,0.99)],
    color=:limegreen,label="",linealpha=0.5,linestyle=:dashdot) |> display
    
end