# Checking consistency of Pareto distribution to empirical MSM partnership datasets

using Pkg
Pkg.activate("../")
using DataFrames
using Distributions
using Plots
using StatsPlots
using Random

studies=DataFrame(
    datasets=["Natsal 2: Schneeberger et al. 2004","Natsal 3: Whittles et al. 2019","GRASP: Whittles et al. 2019"],
    maxpartners=[250,100,1000], samplesize = [138, 188, 691], 
    Pareto_parameter = [1.6, 1.81, 1.60])
display(studies)

function checkmax(datarow)
    pd=Pareto(datarow.Pareto_parameter-1)
    maxs=[maximum(rand(pd,datarow.samplesize)) for i in 1:1000]
    (datamax=datarow.maxpartners, samplemax=maxs)
end
        
Random.seed!(2022)
checkmaxs=checkmax.(eachrow(studies))
DataFrame(quantile.(getindex.(checkmaxs,2), Ref([0,0.025,0.25,0.5,0.75,0.975,1]))|>Base.Fix1(reduce,hcat)|>permutedims,Symbol.(string.([0,2.5,25,50,75,97.5,100]).*"%")) |> display