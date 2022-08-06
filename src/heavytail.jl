using LogExpFunctions
using LinearAlgebra
using StatsBase
using SpecialFunctions
using Distributions
using Random
PoissonCap(λ)=Poisson(min(λ,typemax(Int32)))

function simulatetillsize(edeg::Union{<:Sampleable,<:AbstractArray}, size::Integer, SAR::Number = 0.1; seed = 1, iter = 100000, init_deg=nothing, incidence=false, total=false, samples = 1000)
    n_extinct=0
    incidencerecord = [Int[] for i in 1:samples]
    totalrecord = fill(0,samples)
    for i in 1:iter
        totalcases=seed
        activecases=seed
        if !isnothing(init_deg)
            # first one generation with init_deg instead of edeg
            sum_deg = max.(0, rand.(PoissonCap.(rand(init_deg,activecases)))) |>sum # sum of excess degrees of seeds. prob = cdf(x+1)-cdf(x) if edeg is a continuous distribution.
            offsprings=Binomial(sum_deg,SAR)|>rand
            if(offsprings==0) # extinct
                n_extinct+=1
                if total totalrecord[Int((i-1)÷(iter/samples)+1)]=totalcases end
                continue
            end
            totalcases+=offsprings
            activecases=offsprings
            if incidence push!(incidencerecord[i],activecases) end
        end
        while(totalcases<size)
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
            if total totalrecord[Int((i-1)÷(iter/samples)+1)]=totalcases end
    end
    if incidence||total return((prob=1-n_extinct/iter, incidencerecord=incidencerecord, totalrecord = totalrecord)) end
    1-n_extinct/iter
end

function simulatetillsizehetero(edeg1::Union{<:Sampleable,<:AbstractArray},edeg2::Union{<:Sampleable,<:AbstractArray}, size::Integer, SAR::Number = 0.1; seed = 1, iter = 100000, init_deg=nothing, incidence=false, total=false, samples = 1000)
    n_extinct=0
    incidencerecord = [Int[] for i in 1:samples]
    totalrecord = fill(0,samples)
    for i in 1:iter
        totalcases=seed
        activecases=seed
        issex1 = true
        if !isnothing(init_deg)
            # first one generation with init_deg instead of edeg
            sum_deg = max.(0, rand.(PoissonCap.(rand(init_deg,activecases)))) |>sum # sum of excess degrees of seeds. prob = cdf(x+1)-cdf(x) if edeg is a continuous distribution.
            offsprings=Binomial(sum_deg,SAR)|>rand
            if(offsprings==0) # extinct
                n_extinct+=1
                if total totalrecord[Int((i-1)÷(iter/samples)+1)]=totalcases end
                continue
            end
            totalcases+=offsprings
            activecases=offsprings
            if incidence push!(incidencerecord[i],activecases) end
            issex1 = !issex1
        end
        while(totalcases<size)
            edeg=ifelse(issex1,edeg1,edeg2)
            sum_edeg = max.(0, rand.(PoissonCap.(rand(edeg,activecases))).-1) |>sum # sum of excess degrees of seeds. prob = cdf(x+1)-cdf(x) if edeg is a continuous distribution.
            offsprings=Binomial(sum_edeg,SAR)|>rand
            if(offsprings==0) # extinct
                n_extinct+=1
                break
            end
            totalcases+=offsprings
            activecases=offsprings
            if incidence push!(incidencerecord[i],activecases) end
            issex1 = !issex1
        end
            if total totalrecord[Int((i-1)÷(iter/samples)+1)]=totalcases end
    end
    if incidence||total return((prob=1-n_extinct/iter, incidencerecord=incidencerecord, totalrecord = totalrecord)) end
    1-n_extinct/iter
end

function simulatetillsize(edeg1::Union{<:Sampleable,<:AbstractArray},edeg2::Union{<:Sampleable,<:AbstractArray}, size::Integer, SAR::Number = 0.1; seed = 1, iter = 100000, init_deg1=nothing,init_deg2=nothing, incidence=false, total=false, samples = 1000)
    seed1=rand(Binomial(seed,0.5))
    sim1=simulatetillsizehetero(edeg1,edeg2, size, SAR; seed = seed1, iter = iter, init_deg=init_deg1, incidence=incidence, total=true, samples = samples)
    sim2=simulatetillsizehetero(edeg1,edeg2, size, SAR; seed = seed-seed1, iter = iter, init_deg=init_deg2, incidence=incidence, total=true, samples = samples)
    totalrecord=sim1.totalrecord.+sim2.totalrecord
    if incidence||total return((prob=mean(<(size),totalrecord), incidencerecord=[sim1.incidencerecord;sim2.incidencerecord], totalrecord = totalrecord)) end
    mean(≥(size),totalrecord)
end

simulatetillsize(edeg::Union{<:Sampleable,<:AbstractArray}, size::Integer, SAR::AbstractArray; seed = 1, iter = 100000, init_deg=nothing, incidence=false, samples = 1000) = simulatetillsize.(Ref(edeg), size,SAR; seed = seed, iter = iter, init_deg=init_deg, incidence=incidence, samples = samples)

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

TWquantile(tw::Truncated{<:Weibull},x::Real)=(1/tw.untruncated.α)*(tw.untruncated.α*log(tw.lower)+log1pexp(tw.untruncated.α*log(tw.untruncated.θ)+log(-log(1-x))-tw.untruncated.α*log(tw.lower)))|>exp
Distributions.quantile(tw::Truncated{<:Weibull},x::Real)= if tw.lower>0 TWquantile(tw,x) else quantile(tw.untruncated,x) end
Distributions.truncated(d::UnivariateDistribution; lower)=truncated(d, lower, Inf)
function nonzeroPoissonquantile(d::UnivariateDistribution,x::Real, samplesize=50000) # quantile of Poisson(x) where x ~ tw
    zeroprob=mean((exp∘(-)),(rand(d,samplesize)))
    quantile(d,zeroprob+(1-zeroprob)*x)
end

function ISR(d::Sampleable, n::Integer, weightfun::Function, samplesize = 10000)# importance sampling resampling
    ss=max(2n,samplesize)
    samples=rand(d, ss)
    weights=Weights(weightfun.(samples))
    sample(samples,weights,n,replace=true)
end
function ISR_Weibull(d::EDTWeibull, n::Integer, weightfun::Function = identity, samplesize = 10000)
    # importance sampling resampling for excess degree with truncated weibull degree dist. Default: ~ x * pdf(x)
    ss=max(2n,samplesize)
    samples=rand_TWeibull(d.d, ss)
    weights= Weights(weightfun.(samples))
    sample(samples,weights,n,replace=true)
end

function ISR_Weibull_finite(d::EDTWeibull, n::Integer, weightfun::Function = identity, samplesize = 1000000)
    samples=rand_TWeibull(d.d, samplesize)
    weights= Weights(weightfun.(samples))
    sample(samples,weights,n,replace=false)
end

Weibull2(α,κ)=Weibull(α, (α/κ)^(1/α)) # alternative parameterization where κ is the approximate Pareto shape parameter

log_gamma(a, x) = log(gamma_inc(a, x,0)[2]) + loggamma(a) # to avoid error with small alpha

function R0(tweibull::Truncated{<:Weibull}, SAR=1; logvalue=false)
    α,θ,tr_lower, tr_1 = tweibull.untruncated.α, tweibull.untruncated.θ, (tweibull.lower/tweibull.untruncated.θ)^tweibull.untruncated.α, 1/tweibull.untruncated.θ^tweibull.untruncated.α
    if tweibull.lower≥1
        logR0p1=log(θ)+log_gamma((α+2)/α,tr_lower)-log_gamma((α+1)/α,tr_lower)
        if logvalue return(logexpm1(logR0p1)) end
        return(exp(logR0p1)-1)
    else
        logR0=log_gamma((α+1)/α,tr_1) - log_gamma((α+1)/α,tr_lower) +
        logexpm1(max(0.0,log(θ) + log_gamma((α+2)/α,tr_1)-log_gamma((α+1)/α,tr_1) ))
        if logvalue return(logR0) end
        return(exp(logR0))
    end
end

function Distributions.mean(tweibull::Truncated{<:Weibull})
    α,θ,tr_lower = tweibull.untruncated.α, tweibull.untruncated.θ, (tweibull.lower/tweibull.untruncated.θ)^tweibull.untruncated.α
    exp(log(θ)+log_gamma((α+1)/α,tr_lower)-log_gamma(1,tr_lower))
end

meanm1sampling(d;iter=100000)=mean(max.(0,rand(d,iter).-1)) # evaluates E[max(0,x-1)|x ~ d]
κ2θ(α,κ)=(α/κ)^(1/α)

function eigenNGM(SAR,relR,p_m2o,R_msm,R_hem,R_hew,r_hem,r_hew)
    ngm=[SAR*R_msm*(1-p_m2o) 0 0 0
        SAR*R_msm*p_m2o/2 0 SAR*R_hew SAR*r_hew/2
        SAR*R_msm*p_m2o/2 SAR*R_hem 0 SAR*r_hem/2
        relR*R_msm relR*R_msm relR*R_msm relR*R_msm]
    ev=eigen(ngm)
    (value=ev.values[end],vector=normalize(ev.vectors[:,end],1))
end
