using Distributions
using LinearAlgebra
using Measures
using Plots
using Random
using SpecialFunctions
using StatsPlots
include("heavytail.jl")
include("outbreaklikelihood.jl")
gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=10,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300))
    
# Util function
Base.:-(x::Tuple{Float64,Float64,Float64},y::Int)=x.-y

# Parameters
## Weibull parameters
if !isdefined(@__MODULE__,:κ_msm) κ_msm=0.6 end # Pareto exponent
if !isdefined(@__MODULE__,:α_msm) α_msm=0.15 end
if !isdefined(@__MODULE__,:κ_other) κ_other=1.66 end # Pareto exponent
if !isdefined(@__MODULE__,:κ_other) α_other=0.04 end
## Epi parameters
if !isdefined(@__MODULE__,:inf_period) inf_period=21 end
## Demography parameters
if !isdefined(@__MODULE__,:p_msm) p_msm = 0.01 end# proportion of sexually active MSM
if !isdefined(@__MODULE__,:p_other) p_other = 0.86 end # proportion of sexually active non-MSM

## 2022 outbreak parameters
#currentcases = 584 - 16 # as of 29 May 2022. total - discarded from Global.health
if !isdefined(@__MODULE__,:currentcases) currentcases = 776 - 48 end # as of 31 May 2022. total - discarded from Global.health

## crontrol variables
if !isdefined(@__MODULE__,:printtable) printtable = true end
if !isdefined(@__MODULE__,:outbreak_iters) outbreak_iters = 100000 end
if !isdefined(@__MODULE__,:hist_upp) hist_upp = :auto end

# Degree distribution over 1 year
deg1y_msm = truncated(Weibull2(α_msm, κ_msm), lower=1)
deg1y_other = truncated(Weibull2(α_other, κ_other), lower=1)

# Preparation for draws from ~ x*pdf(x)
ED_3w_msm = truncated(Weibull(params(deg1y_msm)[1], params(deg1y_msm)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull 
ED_3w_other = truncated(Weibull(params(deg1y_other)[1], params(deg1y_other)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull

# Pre-draws from excess degree +1 distribution for efficiency
Random.seed!(2022)
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
if !isdefined(@__MODULE__,:R0only) R0only = false end
if R0only
    sarrange=exp.([-10:0.1:-0.5;-0.4:0.01:0])
    R0values=[R0(deg3w_msm) R0(deg3w_other)].*sarrange;
else


# Degree over 1 year
Random.seed!(2022)
D_msm=[floor.(Int,rand_TWeibull(deg1y_msm,10000)) for i in 1:100]
D_other=[floor.(Int,rand_TWeibull(deg1y_other,10000)) for i in 1:100]
freqs=[sum.([(==).(1:49);>(49)],Ref(x))./10000 for x in (D_msm[1], D_other[1])]

labels=["MSM" "non-MSM"]
groupedbar(reduce(hcat,freqs),label=labels, linealpha=0, 
    xlabel = "number of partners over 1 year", ylabel="relative frequency")|>display
#capnan(x,cap)=ifelse.(x.<cap,NaN,x)
loglogplot=plot([sort.(D_msm),sort.(D_other)], [(length(D_msm[1]):-1:1)./length(D_msm[1]), (length(D_other[1]):-1:1)./length(D_other[1])],xaxis=:log,yaxis=:log,label="",
    xlabel = "number of partners over 1 year", ylabel="relative frequency (cumulative)", color=repeat([1 2],inner=(1,100)),linealpha=0.2)
if isdefined(@__MODULE__,:partnersdata)
    freqs=[sum.([(==).(1:49);>(49)],Ref(vcat(fill.(1:length(x),x)...)))./sum(x) for x in eachcol(partnersdata)]
    groupedbar(reduce(hcat,freqs),label=labels, linealpha=0, 
    xlabel = "number of partners over 1 year", ylabel="relative frequency")|>display

    y=[begin
            nonzero=col.>0
            1 .-((cumsum(col[nonzero]))[1:end-1]./sum(col))
        end for col in eachcol(partnersdata)]
    x=[(1:length(col))[col.>0][2:end] for col in eachcol(partnersdata)]
    plot!(loglogplot, x,y,linealpha=0,markershape=:circ,markersize=3,color=(1:2)',markerstrokewidth=1,label=labels)
end
display(loglogplot)


# Likelihood of an outbreak of current size or greater
if printtable
for SAR in [0.05, 0.1, 0.2, 0.5]
    outbreaklikelihood1(edeg3w_msm, edeg3w_other, currentcases, SAR, iter=outbreak_iters)
    print("\n\n")
    outbreaklikelihood2(edeg3w_msm, edeg3w_other, deg3w_msm, deg3w_msmavg, deg3w_otheravg, currentcases, SAR, iter=outbreak_iters)
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
label=string.(seedrange').*[" initial case" fill("",1,4)],legend=:topleft, ylim=(0,1),
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
initdegs = repeat([nothing, deg3w_otheravg],inner=10)
sarrange = repeat(0.1:0.1:1,2)
majobrisk=[simulatetillsize(edeg3w_other,major_threshold,sar;seed=1000, init_deg=initdeg, total=true,iter=1000) for (initdeg, sar) in zip(initdegs, sarrange)]

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
R0values=[R0(deg3w_msm) R0(deg3w_other)].*sarrange
R0valuesplot=plot(sarrange,R0values,yaxis=:log,ylim=(0.01,100),linewidth=1.5,
xlab="SAR per sexually-associated contact", ylab="R₀ over sexual network",
label=["MSM" "non-MSM"],legend=:topleft)
hline!(R0valuesplot,[1],linestyle=:dot,color=:black,label="") |> display

R0valuesplot=plot(sarrange,max.(0,min.(1,1 .- 1 ./R0values)),ylim=(0,1),linewidth=1.5,
xlab="SAR per sexually-associated contact", ylab="relative reduction in R₀ for control",
label=["MSM" "non-MSM"],legend=:topleft)|>display

# Degree distributions with superspreading control
Random.seed!(2022)
αrange=0.04:0.01:100
κrange=0.6
qs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=1),0.99) for α in αrange, κ in κrange]
findq(q,qs)=findmin(x->abs(x-q),qs)
αs=αrange[last.(findq.([25,15],Ref(qs)))]
origdeg=floor.(Int,rand(truncated(deg3w_msm.untruncated,lower=1),100))
orig_1p=round(quantile(truncated(deg3w_msm.untruncated,lower=1),0.99),digits=1)#round(quantile(origdeg,0.99),digits=0)
degs=[floor.(Int,rand(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=1),100)) for α in αs, κ in κrange]
freqs=[sum.([(==).(1:49);>(49)],Ref(x))./1000 for x in [[origdeg];degs]]
plot(bar.(freqs, linealpha=0, 
     ylabel="",legend=:none)...,
    xticks=[:none :none :native],
    ylim=(-0.01,hist_upp),xlim=(-0.4,51),layout=(3,1),size=(400,300),
label=["baseline, 1%tile: $orig_1p" "1%tile: 15" "1%tile: 25"],legend=:topright,
xlabel = ["" "" "number of partners over 21 days"],bottom_margin=[-2mm -2mm 4mm]) |> display

# R0 vs 1%ile
Random.seed!(2022)
αrange=0.04:0.01:100
sarrange=[0.05,0.1,0.2,0.5]
κrange=repeat([κ_msm],length(sarrange))
R0s=repeat(sarrange,inner=length(κrange)÷length(sarrange))'.*max.(1e-15,[R0(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365),lower=inf_period/365),logvalue=false) for α in αrange, κ in κrange])
xs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=1),0.99) for α in αrange, κ in κrange]
R0capplot=plot(xs,R0s,yaxis=:log,ylim=(0.01,100),xlim=(0,50), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="R₀ over MSM sexual network",
legend=:topleft, label="SAR: ".*string.(sarrange'),
palette=:vikO,color=(1:4)'.*50,linewidth=1.5)
hline!(R0capplot,[1],linestyle=:dot,color=:black,linealpha=1,label="")
vline!(R0capplot,[quantile(truncated(deg3w_msm.untruncated,lower=1),0.99)],color=:green,label="",linealpha=0.5,linestyle=:dash) |> display

# Reduction in R0 for control
Random.seed!(2022)
αrange=0.04:0.001:100
sarrange=[0.05,0.1,0.2,0.5]
κrange=repeat([κ_msm],length(sarrange))
R0s=repeat(sarrange,inner=length(κrange)÷length(sarrange))'.*max.(1e-15,[R0(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365),lower=inf_period/365),logvalue=false) for α in αrange, κ in κrange])
xs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=1),0.99) for α in αrange, κ in κrange]
R0capplot=plot(xs,max.(0,min.(1,1 .- 1 ./R0s)),ylim=(0,1),xlim=(0,50), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="relative reduction in R₀ for control",
legend=:topleft, label="SAR: ".*string.(sarrange'),
palette=:vikO,color=(1:4)'.*50,linewidth=1.5)
vline!(R0capplot,[quantile(truncated(deg3w_msm.untruncated,lower=1),0.99)],color=:green,label="",linealpha=0.5,linestyle=:dash)|>display

# Eigenvector analysis
R_msm=R0(deg3w_msm)#meanm1sampling(edeg3w_msm)
R_other=R0(deg3w_other)#meanm1sampling(edeg3w_other)
r_other=mean(deg3w_other)#mean(rand(deg3w_otheravg,1000000))

Rlist=[R_msm,R_other,r_other]
relRrange=0.:0.1:0.5
prange=[0,0.01,0.05,0.1]
@time eigs=[eigenNGM(1,R,p_m2o,Rlist...) for R in relRrange, p_m2o in prange];
evectors=getindex.(eigs,:vector)
plot(relRrange,last.(evectors),
xlabel="reproduction number ratio",ylabel="non-sexually-associated case ratio",
label=["MSM to non-MSM mixing: " "" "" ""].*string.(prange'),legend=:topleft) |> display

## Sensitivity analysis of the degree distribution shape
# MSM
αrange=0.04:0.01:100
κrange=repeat(0.6:0.1:1,2)
R0s=repeat([0.1,0.5],inner=length(κrange)÷2)'.*max.(1e-15,[R0(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365),lower=inf_period/365),logvalue=false) for α in αrange, κ in κrange])
xs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=1),0.99) for α in αrange, κ in κrange]
R0capplot=plot(xs,R0s,yaxis=:log,ylim=(0.01,100),xlim=(0,100), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="R₀ over sexual network",
legend=:bottomright, label="κ: ".*string.(κrange'),
palette=:vikO,color=([1:5;10:-1:6])'.*22,linewidth=1.5)
hline!(R0capplot,[1],linestyle=:dot,linealpha=1,color=:black,label="")
vline!(R0capplot,[quantile(truncated(deg3w_msm.untruncated,lower=1),0.99)],
    color=:green,label="",linealpha=0.5,linestyle=:dash) |> display 

# non-MSM

αrange=BigFloat.(10 .^(-3:0.1:0))
κrange=repeat(1.5:0.1:1.8,2)
@time R0s=repeat([0.1,0.5],inner=length(κrange)÷2)'.*max.(1e-15,[R0(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365),lower=inf_period/365)) for α in αrange, κ in κrange])
xs= [quantile(truncated(Weibull(α,κ2θ(α,κ)*inf_period/365);lower=1),0.99) for α in αrange, κ in κrange]
R0capplot=plot(xs,R0s,yaxis=:log,ylim=(0.001,10),xlim=(0,20), 
xlabel="# of partners over 21 days at 1st percentile", ylabel="R₀ over sexual network",
legend=:bottomright, label="κ: ".*string.(κrange'),
palette=:vikO,color=([1:4;8:-1:5])'.*27,linewidth=1.5)
hline!(R0capplot,[1],linestyle=:dot,linealpha=1,color=:black,label="")
vline!(R0capplot,[TWquantile(truncated(deg3w_other.untruncated,lower=1),0.99)],
    color=:green,label="",linealpha=0.5,linestyle=:dash) |> display
    
end