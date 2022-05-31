using LogExpFunctions
using StatsBase
using Distributions
PoissonCap(λ)=Poisson(min(λ,typemax(Int32)))

function simulatetillsize(edeg::Any, size::Integer, SAR::Number = 0.1; seed = 1, iter = 100000, init_deg=nothing, incidence=false, total=false, samples = 1000)
    n_extinct=0
    incidencerecord = [Int[] for i in 1:samples]
    totalrecord = fill(0,samples)
    for i in 1:iter
        totalcases=seed
        activecases=seed
        if !isnothing(init_deg)
            # first one generation with init_deg instead of edeg
            #sum_deg = max.(0, floor.(Int,rand(init_deg,activecases))) |>sum # sum of excess degrees of seeds. prob = cdf(x+1)-cdf(x) if edeg is a continuous distribution.
            sum_deg = max.(0, rand.(PoissonCap.(rand(init_deg,activecases)))) |>sum # sum of excess degrees of seeds. prob = cdf(x+1)-cdf(x) if edeg is a continuous distribution.
            offsprings=Binomial(sum_deg,SAR)|>rand
            if(offsprings==0) # extinct
                n_extinct+=1
                if total totalrecord[i]=totalcases end
                continue
            end
            totalcases+=offsprings
            activecases=offsprings
            if incidence push!(incidencerecord[i],activecases) end
        end
        while(totalcases<size)
            #sum_edeg = max.(0, floor.(Int,rand(edeg,activecases)).-1) |>sum # sum of excess degrees of seeds. prob = cdf(x+1)-cdf(x) if edeg is a continuous distribution.
            sum_edeg = max.(0, rand.(PoissonCap.(rand(edeg,activecases))).-1) |>sum # sum of excess degrees of seeds. prob = cdf(x+1)-cdf(x) if edeg is a continuous distribution.
            offsprings=Binomial(sum_edeg,SAR)|>rand
            if(offsprings==0) # extinct
                n_extinct+=1
                break
            end
            totalcases+=offsprings
            activecases=offsprings
            if incidence push!(incidencerecord[i],activecases) end
        end
            if total totalrecord[i]=totalcases end
    end
    if incidence||total return((prob=1-n_extinct/iter, incidencerecord=incidencerecord, totalrecord = totalrecord)) end
    1-n_extinct/iter
end

simulatetillsize(edeg::Any, size::Integer, SAR::AbstractArray; seed = 1, iter = 100000, init_deg=nothing, incidence=false, samples = 1000) = simulatetillsize.(Ref(edeg), size,SAR; seed = seed, iter = iter, init_deg=init_deg, incidence=incidence, samples = samples)

struct EDTWeibull{T<:Truncated{<:Weibull},F}
    d::T
    weightfun::F
    EDTWeibull(d, weightfun = identity)=new{typeof(d),typeof(weightfun)}(d, weightfun)
end
Base.rand(tw::Truncated{<:Weibull})=rand_TWeibull(tw)
Base.rand(rng::AbstractRNG,tw::Truncated{<:Weibull})=rand_TWeibull(rng,tw)
Base.rand(tw::Truncated{<:Weibull},n::Int)=rand_TWeibull(tw,n)

Base.rand(edtw::EDTWeibull,n::Integer)=ISR_Weibull(edtw, n)
function rand_TWeibull(d::Truncated{<:Weibull},n::Integer) # inverse method: useful when α is near 1 or greater
    if(d.lower == d.upper) return(d.lower) end
    α, θ = d.untruncated.α, d.untruncated.θ
    logqlower = (-(d.lower/θ)^α)
    logqupper = (-(d.upper/θ)^α)
    us = rand(n)
    log_uoverl = logqupper-logqlower
    logqs = logqlower.+log1mexp.(log.(us).+log1mexp(log_uoverl)) # 1 .- us.*(1-uoverl)
    (-θ^α.*logqs).^(1/α)
end
function rand_TWeibull(d::Truncated{<:Weibull}) # inverse method: useful when α is near 1 or greater
    if(d.lower == d.upper) return(d.lower) end
    α, θ = d.untruncated.α, d.untruncated.θ
    logqlower = (-(d.lower/θ)^α)
    logqupper = (-(d.upper/θ)^α)
    us = rand()
    log_uoverl = logqupper-logqlower
    logqs = logqlower+log1mexp(log(us)+log1mexp(log_uoverl)) # 1 .- us.*(1-uoverl)
    (-θ^α*logqs)^(1/α)
end
function rand_TWeibull(rng::AbstractRNG, d::Truncated{<:Weibull}) # inverse method: useful when α is near 1 or greater
    if(d.lower == d.upper) return(d.lower) end
    α, θ = d.untruncated.α, d.untruncated.θ
    logqlower = (-(d.lower/θ)^α)
    logqupper = (-(d.upper/θ)^α)
    us = rand(rng)
    log_uoverl = logqupper-logqlower
    logqs = logqlower+log1mexp(log(us)+log1mexp(log_uoverl)) # 1 .- us.*(1-uoverl)
    (-θ^α*logqs)^(1/α)
end

function ISR(d::Sampleable, n::Integer, weightfun::Function, samplesize = 1000)# importance sampling resampling
    ss=max(2n,samplesize)
    samples=rand(d, ss)
    weights=Weights(weightfun.(samples))
    sample(samples,weights,n,replace=true)
end
function ISR_Weibull(d::EDTWeibull, n::Integer, weightfun::Function = identity, samplesize = 1000)
    # importance sampling resampling for excess degree with truncated weibull degree dist. Default: ~ x * pdf(x)
    ss=max(2n,samplesize)
    samples=rand_TWeibull(d.d, ss)
    weights= Weights(weightfun.(samples))
    sample(samples,weights,n,replace=true)
end
Weibull2(α,κ)=Weibull(α, (α/κ)^(1/α)) # alternative parameterization where κ is the approximate Pareto shape parameter
function Weibull3(λ,κ) # alternative parameterization where λ = θ^α
    α=λ*κ
    θ=(α/κ)^(1/α)
    Weibull(α, θ)
end

function R0(tweibull::Truncated{<:Weibull}, SAR=1; logvalue=false)
    α,θ,tr_lower, tr_1 = tweibull.untruncated.α, tweibull.untruncated.θ, (tweibull.lower/tweibull.untruncated.θ)^tweibull.untruncated.α, 1/tweibull.untruncated.θ^tweibull.untruncated.α
    if tweibull.lower≥1
        logR0p1=log(θ)+loggamma((α+2)/α,tr_lower)-loggamma((α+1)/α,tr_lower)
        if logvalue return(logexpm1(logR0p1)) end
        return(exp(logR0p1)-1)
    else
        logR0=loggamma((α+1)/α,tr_1) - loggamma((α+1)/α,tr_lower) +
        logexpm1(max(0.0,log(θ) + loggamma((α+2)/α,tr_1)-loggamma((α+1)/α,tr_1) ))
        if logvalue return(logR0) end
        return(exp(logR0))
    end
end
R0(α,κ,lower)=R0(truncated(Weibull2(α,κ),lower=lower))
meanm1sampling(d;iter=100000)=mean(max.(0,rand(d,iter).-1))
κ2θ(α,κ)=(α/κ)^(1/α)

function eigenNGM(SAR,relR,p_m2o,R_msm,R_other,r_other)
    ngm=[SAR*R_msm*(1-p_m2o) 0 0
        SAR*R_msm*p_m2o SAR*R_other SAR*r_other
        relR*R_msm relR*R_msm relR*R_msm]
    ev=eigen(ngm)
    (value=ev.values[end],vector=normalize(ev.vectors[:,end],1))
end