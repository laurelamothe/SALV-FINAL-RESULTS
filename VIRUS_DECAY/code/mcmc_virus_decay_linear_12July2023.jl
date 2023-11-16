######################################################
#####  Laure L. & David D. / Main Exp. SALINITY  #####
########  mcmc to fit a linear regression of  ########  
####### virus decay as a function of salinity ########  
############# 25 July 2023 - Banyuls #################
######################################################

## PACKAGES
using Distributed
addprocs(4)

@everywhere using Dates, Turing, Distributions, DifferentialEquations, Statistics, CSV, DataFrames, Glob


@everywhere begin
    # charge the data (virus decay at different salinity, calculated from cytometry counts overtime)
    data = DataFrame(CSV.File("./VIRUS_DECAY/data/data_input_virus_decay.csv"))

    # definition of the linear regression
    function regression(a, b, Sal)
        return a .* Sal .+ b
    end

end


## function to fit the regression to the data
@everywhere begin

    @model function fit_virus(data_virus, data_salinity)

        ## Priors

        a ~ truncated(LogNormal(log(0.0005), 1), 0, 0.001)
        b ~ truncated(LogNormal(log(0.01), 1), 0, 0.05)


        ## standard deviation of the likelyhood
        σ ~ InverseGamma(10, 0.001)

        ## modelling virus decay
        #predicted_reg = regression(a, b, data_salinity)
        predicted_reg = regression(a, b, data_salinity)

        ## Likelihood
        for i = 1:length(data_virus)
            data_virus[i] ~ Normal(predicted_reg[i], σ)
        end

    end
end

## MCMC run
nburn = 2000 # number of iteration available for the chains to converge
nstep = 2000 # number of iteration after the burn
nchain = 4 # number of independant chains

## Fittin procedure
modelfit = fit_virus(data.obscurity, data.Salinity)
chainfit = reduce(chainscat, pmap(x -> sample(modelfit, NUTS(nburn, 0.65), nstep, save_state=false), 1:nchain))

## EXPORT
chainarray = Array(chainfit)

chain_df = DataFrame(
    a=chainarray[:, 1],
    b=chainarray[:, 2])

stats_chain = ess_rhat(chainfit)

# Save Arrays as CSV
CSV.write("./VIRUS_DECAY/results/SALEXP_virus_decay_dark_" * string(today()) * ".csv", chain_df, writeheader=true)
CSV.write("./VIRUS_DECAY/results/SALEXP_statchains_virus_decay_dark_" * string(today()) * ".csv", DataFrame(stats_chain), writeheader=true)


ll = [loglikelihood(modelfit, j) for j in eachrow(DataFrame(chainfit))]
k = size(chain_df)[2] - 1
n = length(data.light)

AIC = 2 * k - 2 * mean(ll)
AICc = 2 * k - 2 * mean(ll) + (2 * k) * (k + 1) / (n - k - 1)
BIC = k * log(n) - 2 * mean(ll)

