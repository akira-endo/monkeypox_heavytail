using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributions
using LinearAlgebra
using Measures
using Plots
using Random
using SpecialFunctions
using StatsPlots
include("src/heavytail.jl")
include("src/outputs.jl")
gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=10,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))
    
# Util function
Base.:-(x::Tuple{Float64,Float64,Float64},y::Int)=x.-y

# Parameters
## Weibull parameters
κ_msm=0.6 # Pareto exponent
α_msm=0.15
κ_other=1.66 # Pareto exponent
α_other=0.04

## Epi parameters
inf_period=21

## Demography parameters
p_msm = 0.01 # proportion of sexually active MSM
p_other = 0.86 # proportion of sexually active non-MSM

## 2022 outbreak parameters
currentcases = 584 - 16 # as of 29 May 2022. total - discarded from Global.health

# Degree distribution over 1 year
deg1y_msm = truncated(Weibull2(α_msm, κ_msm), lower=1)
deg1y_other = truncated(Weibull2(α_other, κ_other), lower=1)

# Preparation for draws from ~ x*pdf(x)
ED_3w_msm = truncated(Weibull(params(deg1y_msm)[1], params(deg1y_msm)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull 
ED_3w_other = truncated(Weibull(params(deg1y_other)[1], params(deg1y_other)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull

# Pre-draws from excess degree +1 distribution for efficiency
drawsize = 1000000
edeg3w_msm = rand(ED_3w_msm, drawsize)
edeg3w_other = rand(ED_3w_other, drawsize)

# Degree distributions for different types
deg3w_msm = truncated(Weibull(params(deg1y_msm)[1], params(deg1y_msm)[2]*inf_period/365), lower=inf_period/365)
deg3w_other = truncated(Weibull(params(deg1y_other)[1], params(deg1y_other)[2]*inf_period/365), lower=inf_period/365)
contDirac = truncated(Normal(0,1e-10),lower=0,upper=0.1)
deg3w_msmavg = MixtureModel([contDirac,deg3w_msm],[1-p_msm ,p_msm])
deg3w_otheravg = MixtureModel([contDirac,deg3w_other],[1-p_other ,p_other])


### Outputs

# Degree over 1 year
Random.seed!(2022)
D_msm=floor.(Int,rand_TWeibull(deg1y_msm,500))
D_other=floor.(Int,rand_TWeibull(deg1y_other,500))
freqs=[sum.((==).(1:50),Ref(x))./500 for x in (D_msm, D_other)]

labels=["MSM" "non-MSM"]
groupedbar(reduce(hcat,freqs),label=labels, linealpha=0, 
    xlabel = "simulated number of partners over 1 year", ylabel="relative frequency")|>display
plot(sort.([D_msm,D_other]), (length(D_msm):-1:1)./length(D_msm),xaxis=:log,yaxis=:log,label=labels,
    xlabel = "simulated number of partners over 1 year", ylabel="relative frequency (cumulative)")


# Likelihood of an outbreak of current size or greater
for SAR in [0.05, 0.1, 0.2, 0.5]
    output1(edeg3w_msm, edeg3w_other, currentcases, SAR)
    print("\n")
    output2(edeg3w_msm, edeg3w_other, deg3w_msm, deg3w_msmavg, deg3w_otheravg, currentcases, SAR)
    print("\n")
end


## Main plots

# The likelihood of outbreak vs SAR, sexually-associated MSM
Random.seed!(2022)
sarrange = [0;0.005;0.01:0.01:0.05;0.075;0.1:0.05:1]
seedrange = 1:5
currobsize=[simulatetillsize(edeg3w_msm,currentcases,sar;seed=s,iter=100000) for sar in sarrange, s in seedrange]

plot(sarrange,currobsize,xlim=(0,1),
label=string.(seedrange').*[" initial case" fill("",1,4)],legend=:bottomright#=(0.12,0.9)=#, ylim=(0,1),
xlabel="SAR per sexually-associated contact", ylabel="risk of outbreak (≥ $currentcases cases)", title="MSM index cases: sexually-associated exposure",
linewidth=1.5,color=(5:-1:1)'.*40,palette=:romaO) |> display

# The likelihood of outbreak vs SAR, random MSM
Random.seed!(2022)
sarrange = [0.02:0.02:0.08;0.1:0.1:1]
seedrange = 1:5
currobsize_rMSM=[simulatetillsize(edeg3w_msm,currentcases,sar;seed=s, init_deg=deg3w_msm) for sar in sarrange, s in seedrange]

plot(sarrange,currobsize_rMSM,xlim=(0,1),
label=string.(seedrange').*[" initial case" fill("",1,4)],legend=(0.12,0.9), ylim=(0,1),
xlabel="SAR per sexually-associated contact", ylabel="risk of outbreak (≥ $currentcases cases)",title="MSM index cases: other exposure",
linewidth=1.5,color=(5:-1:1)'.*40,palette=:romaO) |> display

# The likelihood of outbreak vs SAR, random gen pop

Random.seed!(2022)
sarrange = [0.01:0.01:0.04;0.05;0.1:0.1:1]
seedrange = [5,10,50,100,500]
currobsize_ravgMSM=[simulatetillsize(edeg3w_msm,currentcases, sar;seed=s, init_deg=deg3w_msmavg) for sar in sarrange, s in seedrange]

plot(sarrange,currobsize_ravgMSM_highseed,xlim=(0,1),
label=string.(seedrange').*[" initial cases" fill("",1,4)],legend=(0.12,0.9), ylim=(0,1),
xlabel="SAR per sexually-associated contact", ylabel="risk of outbreak (≥ $currentcases cases)",title="General index cases: other exposure",
linewidth=1.5,color=(5:-1:1)'.*50,palette=:romaO) |> display

# Outbreak sizes in non-MSM
Random.seed!(2022)
major_threshold = 10000
initdegs = repeat([nothing, deg3w_otheravg],inner=10)
sarrange = repeat(0.1:0.1:1,2)
majobrisk=[simulatetillsize(edeg3w_other,major_threshold,sar;seed=1000, init_deg=initdeg, total=true,iter=1000) for (initdeg, sar) in zip(initdegs, sarrange)]

majobquant = quantile.(getfield.(majobrisk,:totalrecord),Ref((0.025,0.5,0.975))).-1000
plot(0.1:0.1:1,reshape(getindex.(majobquant,2),10,2),ribbon=(reshape(getindex.(majobquant,2).-getindex.(majobquant,1),10,2),reshape(getindex.(majobquant,3).-getindex.(majobquant,2),10,2)) ,
markershape=:circ,markerstrokewidth=0,color=[1 2],xlim=:auto, ylim=:auto,
label="1000".*[" initial cases: sexually-associated" " initial cases: other"],legend=:topleft,
xlabel="SAR per sexually-associated contact", ylabel="Outbreak size",title="Outbreak size in non-MSM sexual network",
linewidth=1.5) |> display

# R0 vs SAR
sarrange=exp.(-10:0.1:0)
R0values=[meanm1sampling(edeg3w_msm) meanm1sampling(edeg3w_other)].*sarrange;
R0valuesplot=plot(sarrange,R0values,yaxis=:log,ylim=(0.01,100),linewidth=1.5,
xlab="SAR per sexually-associated contact", ylab="R₀ over sexual network",
label=["MSM" "non-MSM"],legend=:topleft)
hline!(R0valuesplot,[1],linestyle=:dot,color=:black,label="") |> display

# Degree distributions with superspreading control
Random.seed!(2022)
αrange=0.04:0.01:100
κrange=0.6
qs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=inf_period/365),0.99) for α in αrange, κ in κrange]
findq(q,qs)=findmin(x->abs(x-q),qs)
αs=αrange[last.(findq.([7.5,5],Ref(qs)))]
origdeg=floor.(Int,rand(deg3w_msm,1000))

degs=[floor.(Int,rand(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=inf_period/365),1000)) for α in αs, κ in κrange]
freqs=[sum.([(==).(1:49);>(49)],Ref(x))./1000 for x in [[origdeg];degs]]
plot(bar.(freqs, linealpha=0, 
     ylabel="",legend=:none)...,
    xticks=[:none :none :native],
    ylim=(-0.01,0.06),xlim=(-0.4,51),layout=(3,1),size=(400,300),
label=["baseline, 1%tile: 9.5" "1%tile: 5" "1%tile: 7.5"],legend=:topright,
xlabel = ["" "" "simulated number of partners over 21 days"],bottom_margin=[-2mm -2mm 4mm]) |> display

# R0 vs 1%ile
Random.seed!(2022)
αrange=0.04:0.01:100
sarrange=[0.05,0.1,0.2,0.5]
κrange=repeat([0.6],length(sarrange))
R0s=repeat(sarrange,inner=length(κrange)÷length(sarrange))'.*max.(1e-15,[R0(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365),lower=inf_period/365),logvalue=false) for α in αrange, κ in κrange])
xs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=inf_period/365),0.99) for α in αrange, κ in κrange]
R0capplot=plot(xs,R0s,yaxis=:log,ylim=(0.01,100),xlim=(0,15), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="R₀ over MSM sexual network",
legend=:topleft, label="SAR: ".*string.(sarrange'),
palette=:vikO,color=(1:4)'.*50,linewidth=1.5)
hline!(R0capplot,[1],linestyle=:dot,color=:black,linealpha=1,label="")
vline!(R0capplot,[quantile(deg3w_msm,0.99)],color=:green,label="",linealpha=0.5,linestyle=:dash) |> display

#
R_msm=meanm1sampling(edeg3w_msm)
R_other=meanm1sampling(edeg3w_other)
r_other=mean(rand(deg3w_otheravg,1000000))

Rlist=[R_msm,R_other,r_other]
relRrange=0.:0.1:0.5
prange=[0,0.01,0.05]
eigs=[eigenNGM(1,R,p_m2o,Rlist...) for R in relRrange, p_m2o in prange];
evectors=getindex.(eigs,:vector)
plot(relRrange,last.(evectors),
xlabel="reproduction number ratio",ylabel="non-sexually-associated case ratio",
label=["sexual bridging ratio: " "" ""].*string.(prange'),legend=:topleft) |> display