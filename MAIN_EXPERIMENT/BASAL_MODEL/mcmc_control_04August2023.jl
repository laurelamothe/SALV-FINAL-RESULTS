######################################################
########## Laure L. / Main Exp. SALINITY  ############
#################### Modelling #######################
##### System - mcmc to estimate all parameters   #####
########### David/Laure 30 May 2023 - Banyuls ########
######################################################

#%% PACKAGES
using Distributed
addprocs(4)
@everywhere using Dates, Turing, Distributions, DifferentialEquations, Statistics, CSV, DataFrames, Glob

#files to use in fitting procedure
@everywhere begin
    #experimental data
    ostreo_list = [
        "A_S1_host_infected.csv",
        "A_S2_host_infected.csv",
        "A_S3_host_infected.csv",
        "A_S4_host_infected.csv",
        "A_S5_host_infected.csv",
        "A_S6_host_infected.csv",
        "B_S1_host_infected.csv",
        "B_S2_host_infected.csv",
        "B_S3_host_infected.csv",
        "B_S4_host_infected.csv",
        "B_S5_host_infected.csv",
        "B_S6_host_infected.csv",
        "C_S1_host_infected.csv",
        "C_S2_host_infected.csv",
        "C_S3_host_infected.csv",
        "C_S4_host_infected.csv",
        "C_S5_host_infected.csv",
        "C_S6_host_infected.csv",
        "D_S1_host_infected.csv",
        "D_S2_host_infected.csv",
        "D_S3_host_infected.csv",
        "D_S4_host_infected.csv",
        "D_S5_host_infected.csv",
        "D_S6_host_infected.csv"]

    @everywhere virus_list = [
        "A_S1_virus.csv",
        "A_S2_virus.csv",
        "A_S3_virus.csv",
        "A_S4_virus.csv",
        "A_S5_virus.csv",
        "A_S6_virus.csv",
        "B_S1_virus.csv",
        "B_S2_virus.csv",
        "B_S3_virus.csv",
        "B_S4_virus.csv",
        "B_S5_virus.csv",
        "B_S6_virus.csv",
        "C_S1_virus.csv",
        "C_S2_virus.csv",
        "C_S3_virus.csv",
        "C_S4_virus.csv",
        "C_S5_virus.csv",
        "C_S6_virus.csv",
        "D_S1_virus.csv",
        "D_S2_virus.csv",
        "D_S3_virus.csv",
        "D_S4_virus.csv",
        "D_S5_virus.csv",
        "D_S6_virus.csv"]
end



@everywhere begin
    #model (differential equation)
    include("./MAIN_EXPERIMENT/BASAL_MODEL/ode_control_02May2023.jl")
    #initial conditions of the model
    t0 = 0.0
    dt = 0.01
    u0 = 1E6
    tf = 14
    tmod = collect(t0:dt:tf)
    tspan = (t0, tf)
    p = [0.7, 1E8, 0]
    #definition of the initial ODE problem
    prob_system = ODEProblem(ODE_ostreo, u0, tspan, p)
end


@everywhere begin

    @model function fit_system(data_ostreo_rep1, data_ostreo_rep2, time_ostreo, prob_system, tmod, dt)

        ## Priors
        r ~ truncated(LogNormal(log(0.6), 1), 0, 2)
        K ~ truncated(LogNormal(log(1E8), 1), 1E6, 1E10)
        ω = 0

        N0_rep1 ~ truncated(LogNormal(data_ostreo_rep1[1], 0.2), 5E5, 1E7)
        N0_rep2 ~ truncated(LogNormal(data_ostreo_rep2[1], 0.2), 5E5, 1E7)
        u0_rep1 = N0_rep1
        u0_rep2 = N0_rep2

        σ_N ~ InverseGamma(10, 1)

        ## parameters
        p_mcmc = [r, K, ω]

        ## Integration
        # Rep1
        prob_rep1 = remake(prob_system, p=p_mcmc, u0=u0_rep1)
        ymod_rep1 = solve(prob_rep1, Tsit5(), abstol=1E-8, reltol=1E-8, saveat=dt)
        ymod_rep1 = reverse(rotl90(hcat(ymod_rep1.u...)), dims=1)
        id_ostreo = indexin(floor.(time_ostreo, digits=2), tmod)
        predicted_ostreo_rep1 = ymod_rep1[id_ostreo, 1]

        # Rep2
        prob_rep2 = remake(prob_system, p=p_mcmc, u0=u0_rep2)
        ymod_rep2 = solve(prob_rep2, Tsit5(), abstol=1E-8, reltol=1E-8, saveat=dt)
        ymod_rep2 = reverse(rotl90(hcat(ymod_rep2.u...)), dims=1)
        id_ostreo = indexin(floor.(time_ostreo, digits=2), tmod)
        predicted_ostreo_rep2 = ymod_rep2[id_ostreo, 1]


        ## LIKELIHOOD
        # rep1
        #ostreo
        for i = 1:length(predicted_ostreo_rep1)
            temp = predicted_ostreo_rep1[i]
            if temp < 0
                temp = 0
            end
            log_predicted_ostreo_rep1 = log.(temp)
            data_ostreo_rep1[i] ~ Normal(log_predicted_ostreo_rep1, σ_N)
        end

        # rep2
        #ostreo
        for i = 1:length(predicted_ostreo_rep2)
            temp = predicted_ostreo_rep2[i]
            if temp < 0
                temp = 0
            end
            log_predicted_ostreo_rep2 = log.(temp)
            data_ostreo_rep2[i] ~ Normal(log_predicted_ostreo_rep2, σ_N)
        end
    end
end

## Download all data and store it internally
# function to download the ostreo data
function ostreo_df(k)
    #ostreo
    dirdata = "./MAIN_EXPERIMENT/BASAL_MODEL/data/"
    data_ostreo = DataFrame(CSV.File(dirdata * ostreo_list[k]))
    data_ostreo_rep1 = log.(data_ostreo.REPLICATE_1)
    data_ostreo_rep2 = log.(data_ostreo.REPLICATE_2)
    time_ostreo = data_ostreo.TIME
    return (data_ostreo_rep1, data_ostreo_rep2, time_ostreo)
end
# function to download the virus data
function virus_df(k)
    #virus
    dirdata = "./MAIN_EXPERIMENT/BASAL_MODEL/data/"
    data_virus = DataFrame(CSV.File(dirdata * virus_list[k]))
    data_virus_rep1 = log.(data_virus.REPLICATE_1)
    data_virus_rep2 = log.(data_virus.REPLICATE_2)
    time_virus = data_virus.TIME
    return (data_virus_rep1, data_virus_rep2, time_virus)
end
# download all the data
df_ostreo_rep1 = [];
df_virus_rep1 = [];
df_time_ostreo = [];

df_ostreo_rep2 = [];
df_virus_rep2 = [];
df_time_virus = [];

for i in 1:length(ostreo_list)
    push!(df_ostreo_rep1, ostreo_df(i)[1])
    push!(df_ostreo_rep2, ostreo_df(i)[2])
    push!(df_time_ostreo, ostreo_df(i)[3])

    push!(df_virus_rep1, virus_df(i)[1])
    push!(df_virus_rep2, virus_df(i)[2])
    push!(df_time_virus, virus_df(i)[3])
end



## MCMC run
nburn = 2000
nstep = 2000
nchain = 4


for i in 1:length(ostreo_list)
    modelfit = fit_system(df_ostreo_rep1[i], df_ostreo_rep2[i], df_time_ostreo[i], prob_system, tmod, dt)

    chainfit = reduce(chainscat, pmap(x -> sample(modelfit, NUTS(nburn, 0.65), nstep, save_state=false), 1:nchain))

    # EXPORT
    chainarray = Array(chainfit)

    chain_df = DataFrame(
        r=chainarray[:, 1],
        K=chainarray[:, 2],
        N0_rep1=chainarray[:, 3],
        N0_rep2=chainarray[:, 4])
    stats_chain = ess_rhat(chainfit)

    # Save Arrays as CSV
    CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/results/SALEXP_chains_control_" * ostreo_list[i][1:4] * "_" * string(today()) * "_genotoul.csv", chain_df, writeheader=true)#
    CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/results/SALEXP_statchains_control_" * ostreo_list[i][1:4] * "_" * string(today()) * "_genotoul.csv", DataFrame(stats_chain), writeheader=true)

    println(ostreo_list[i][1:4] * " is done!")
end