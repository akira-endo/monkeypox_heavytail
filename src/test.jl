using Pkg
Pkg.activate("../")
using Test

include("heavytail.jl")

Random.seed!(2022)

# PoissonCap
@test PoissonCap.(0:100) == Poisson.(0:100)
@test PoissonCap(typemax(Int64)) == Poisson(typemax(Int32))

# ISR for sampling from x * pdf(truncatedWeibull, x)
for θ in 0.1:0.1:2
    edtw = EDTWeibull(truncated(Weibull(1, θ); lower = θ))
    edtwsamples=rand(edtw,100000)
    gammad = truncated(Gamma(2, θ); lower = θ)
    gammasamples=rand(gammad,100000)
    @test isapprox.(quantile(edtwsamples, (0.25,0.5,0.75)), quantile(gammasamples, (0.25,0.5,0.75)); rtol=0.01) |> all
    @test isapprox(mean(edtwsamples), mean(gammasamples); rtol=0.01)
end

# Alternative Weibull constuctor
@test (log1mexp.(logcdf.(Weibull2.(0.1:0.1:1,0.1:0.1:1),1)) .≈ - 1.0) |> all
@test Weibull2.(0.1:0.1:1,0.1:0.1:1) == Weibull.(0.1:0.1:1,κ2θ.(0.1:0.1:1,0.1:0.1:1))

# analytic R0 computation for truncated Weibull
# continuity around lower = 1
for θ in 0.1:0.1:2
    w = Weibull(θ, θ)
    @test R0(truncated(w,lower=1)) ≈  R0(truncated(w,lower=1-1e-10))
end

# Compare with sample mean
for θ in 0.1:0.1:2
    edtw = truncated(Weibull(1, θ); lower = θ)
    gammad = truncated(Gamma(2, θ); lower = θ)
    gammasamples=rand(gammad,100000)
    @test isapprox(R0(edtw), mean(max.(0,gammasamples.-1)); rtol=0.01) || 
        (isapprox(R0(edtw), mean(max.(0,gammasamples.-1)); atol=1e-3) && isapprox(R0(edtw), mean(max.(0,gammasamples.-1)); rtol=0.2))
end
  
for θ in 0.1:0.1:1
    tw = truncated(Weibull(1, θ); lower = θ)
    gammad2 = Gamma(2, θ)
    gammad3 = Gamma(3, θ)
    @test R0(tw) ≈ ((2/θ)*(1-cdf(gammad3,max(1,θ)))-(1-cdf(gammad2,max(1,θ))))/(1-cdf(gammad2,θ))
end

