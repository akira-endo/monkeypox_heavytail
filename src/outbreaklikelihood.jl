function outbreaklikelihood1(edeg_msm, edeg_other, currentcases, SAR; major_outbreak = 10000, major_seed = 1000, iter=100000)
    print("Assumed SAR = $SAR:\n")
    print("prob a single MSM case causes an outbreak of the current size ($currentcases cases) or greater: ")
    simulatetillsize(edeg_msm,currentcases,SAR; iter=iter)|>print
    print("\nprob of a major outbreak (≥ $major_outbreak cases) given the current cases: ")
    totals=simulatetillsize(edeg_msm,major_outbreak,SAR;seed=1,total=true, samples =100000,iter=iter).totalrecord
    sum(≥(major_outbreak),totals)/sum(≥(currentcases),totals) |>print # conditional prob
print("\nprob $major_seed non-MSM cases causing a major outbreak: ")
    simulatetillsize(edeg3w_other,major_outbreak,SAR;seed=major_seed,iter=iter)|>print
end
function outbreaklikelihood2(edeg_msm, edeg_other, init_deg_msm, init_deg_genmsm,init_deg_genother, currentcases, SAR ;major_outbreak = 10000, major_seed = 1000,iter=100000)
    print("Assumed SAR = $SAR:\n")
    print("prob a random MSM causes an outbreak of the current size ($currentcases cases) or greater: ")
    simulatetillsize(edeg_msm,currentcases,SAR; init_deg=init_deg_msm,iter=iter)|>print
    print("\nprob a random case from gen. pop. causes an outbreak of the current size or greater in MSM pop.: ")
    simulatetillsize(edeg_msm,currentcases,SAR; init_deg = init_deg_genmsm,iter=iter)|>print
    print("\nprob of $major_seed cases from gen. pop. causing a major outbreak (≥ $major_outbreak cases) in MSM population: ")
    simulatetillsize(edeg_msm,major_outbreak, SAR;seed=major_seed, init_deg=init_deg_genmsm)|>print
    print("\nprob $major_seed random cases causing a major outbreak in non-MSM population: ")
    simulatetillsize(edeg_other,major_outbreak, SAR;seed=major_seed, init_deg=init_deg_genother,iter=iter)|>print;
end