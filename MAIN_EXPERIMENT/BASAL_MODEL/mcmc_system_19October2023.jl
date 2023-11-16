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

    #estimated parameters of r and K from control cultures
    control_results = [
        "SALEXP_chains_control_A_S1_2023-07-12.csv",
        "SALEXP_chains_control_A_S2_2023-07-12.csv",
        "SALEXP_chains_control_A_S3_2023-07-12.csv",
        "SALEXP_chains_control_A_S4_2023-07-12.csv",
        "SALEXP_chains_control_A_S5_2023-07-12.csv",
        "SALEXP_chains_control_A_S6_2023-07-12.csv",
        "SALEXP_chains_control_B_S1_2023-07-12.csv",
        "SALEXP_chains_control_B_S2_2023-07-12.csv",
        "SALEXP_chains_control_B_S3_2023-07-12.csv",
        "SALEXP_chains_control_B_S4_2023-07-12.csv",
        "SALEXP_chains_control_B_S5_2023-07-12.csv",
        "SALEXP_chains_control_B_S6_2023-07-12.csv",
        "SALEXP_chains_control_C_S1_2023-07-12.csv",
        "SALEXP_chains_control_C_S2_2023-07-12.csv",
        "SALEXP_chains_control_C_S3_2023-07-12.csv",
        "SALEXP_chains_control_C_S4_2023-07-12.csv",
        "SALEXP_chains_control_C_S5_2023-07-12.csv",
        "SALEXP_chains_control_C_S6_2023-07-12.csv",
        "SALEXP_chains_control_D_S1_2023-07-12.csv",
        "SALEXP_chains_control_D_S2_2023-07-12.csv",
        "SALEXP_chains_control_D_S3_2023-07-12.csv",
        "SALEXP_chains_control_D_S4_2023-07-12.csv",
        "SALEXP_chains_control_D_S5_2023-07-12.csv",
        "SALEXP_chains_control_D_S6_2023-07-12.csv"]

    # mesurements of β and η from experimenntal dynamics
    β_η_priors = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/data/beta_and_eta_priors.csv"))
    salinity = [5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40]
end



@everywhere begin
    #model (system of differential equations)
    include("./ode_system_13July2023.jl")
    #function to obtain the estmated parameter on control cultures
    include("./moyenne_param_control_22May2023.jl")
    #mured values of virus decay
    result_virus_decay = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/data/SALEXP_virus_decay_linear_2023-07-12.csv"))
    #initial conditions of the model
    t0 = 0.0
    dt = 0.01
    u0 = [1E6, 0, 0, 1E5]
    tf = 14
    tmod = collect(t0:dt:tf)
    tspan = (t0, tf)
    p = [6.6, 0.01, 1E8, 1E6, 0, 0.5, 1E-9, 1 / 5, 25, 0.0001]
    #definition of the initial ODE problem
    prob_system = ODEProblem(ODE_system, u0, tspan, p)
end


#fitting procedure
@everywhere begin

    @model function fit_system(fixed_r, fixed_K, fixed_δ, β_mean, η_mean, data_ostreo_rep1, data_ostreo_rep2, data_virus_rep1, data_virus_rep2, time_ostreo, time_virus, prob_system, tmod, dt)

        ## Priors
        # fixed params
        r = fixed_r
        Ks = fixed_K
        δ = fixed_δ
        ω = 0
        Rper = 0

        # parameters to optimize
        µ ~ truncated(LogNormal(log(0.001), 1), 0, 0.04)
        Kr ~ truncated(LogNormal(log(1e8), 2), 1e4, 1e9)
        ϵ ~ truncated(LogNormal(log(0.7), 2), 0, 1)
        ϕlog ~ Uniform(log(1E-13), log(8E-8))
        ϕ = exp(ϕlog)
        η ~ truncated(LogNormal(log(η_mean), 0.5), 0.25, 20)
        β ~ truncated(LogNormal(log(β_mean), 0.5), 5, 400)

        N0_rep1 ~ truncated(LogNormal(data_ostreo_rep1[1], 0.2), 5E5, 1E7)
        N0_rep2 ~ truncated(LogNormal(data_ostreo_rep2[1], 0.2), 5E5, 1E7)
        V0_rep1 ~ truncated(LogNormal(data_virus_rep1[1], 0.2), 5E5, 1E7)
        V0_rep2 ~ truncated(LogNormal(data_virus_rep2[1], 0.2), 5E5, 1E7)
        u0_rep1 = [(1 - Rper) * N0_rep1, Rper * N0_rep1, 0, V0_rep1]
        u0_rep2 = [(1 - Rper) * N0_rep2, Rper * N0_rep2, 0, V0_rep2]


        σ_N ~ InverseGamma(10, 1)
        σ_V ~ InverseGamma(10, 1)

        ## parameters
        p_mcmc = [r, μ, Ks, Kr, ω, ϵ, ϕ, η, β, δ]

        ## Integration
        # Rep1
        prob_rep1 = remake(prob_system, p=p_mcmc, u0=u0_rep1)
        ymod_rep1 = solve(prob_rep1, AutoTsit5(Rosenbrock23()), maxiters=1e7, abstol=1E-12, reltol=1E-12, saveat=dt) # On integre ce probleme
        ymod_rep1 = reverse(rotl90(hcat(ymod_rep1.u...)), dims=1)
        id_ostreo = indexin(floor.(time_ostreo, digits=2), tmod)
        id_virus = indexin(floor.(time_virus, digits=2), tmod)
        predicted_ostreo_rep1 = ymod_rep1[id_ostreo, 1] .+ ymod_rep1[id_ostreo, 2] .+ ymod_rep1[id_ostreo, 3]
        predicted_virus_rep1 = ymod_rep1[id_virus, 4]

        # Rep2
        prob_rep2 = remake(prob_system, p=p_mcmc, u0=u0_rep2)
        ymod_rep2 = solve(prob_rep2, AutoTsit5(Rosenbrock23()), maxiters=1e7, abstol=1E-12, reltol=1E-12, saveat=dt)           # On integre ce probleme.
        ymod_rep2 = reverse(rotl90(hcat(ymod_rep2.u...)), dims=1)
        id_ostreo = indexin(floor.(time_ostreo, digits=2), tmod)
        id_virus = indexin(floor.(time_virus, digits=2), tmod)
        predicted_ostreo_rep2 = ymod_rep2[id_ostreo, 1] .+ ymod_rep2[id_ostreo, 2] .+ ymod_rep2[id_ostreo, 3]
        predicted_virus_rep2 = ymod_rep2[id_virus, 4]


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

        #virus
        for i = 1:length(predicted_virus_rep1)
            temp = predicted_virus_rep1[i]
            if temp < 0
                temp = 0
            end
            log_predicted_virus_rep1 = log.(temp)
            data_virus_rep1[i] ~ Normal(log_predicted_virus_rep1, σ_V)
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


        # 
        for i = 1:length(predicted_virus_rep2)
            temp = predicted_virus_rep2[i]
            if temp < 0
                temp = 0
            end
            log_predicted_virus_rep2 = log.(temp)
            data_virus_rep2[i] ~ Normal(log_predicted_virus_rep2, σ_V)
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
    fixed_r = get_mean_control(control_results).r[i]
    fixed_K = get_mean_control(control_results).K[i]
    fixed_δ = mean(result_virus_decay.a) * salinity[i] + mean(result_virus_decay.b)
    modelfit = fit_system(fixed_r, fixed_K, fixed_δ, β_η_priors.beta[i], β_η_priors.eta[i], df_ostreo_rep1[i], df_ostreo_rep2[i], df_virus_rep1[i], df_virus_rep2[i], df_time_ostreo[i], df_time_virus[i], prob_system, tmod, dt)

    chainfit = reduce(chainscat, pmap(x -> sample(modelfit, NUTS(nburn, 0.65), nstep, save_state=false), 1:nchain))

    # EXPORT
    chainarray = Array(chainfit)

    chain_df = DataFrame(
        Kr=chainarray[:, 2],
        μ=chainarray[:, 1],
        ϵ=chainarray[:, 3],
        ϕ=chainarray[:, 4],
        η=chainarray[:, 5],
        β=chainarray[:, 6],
        Rper=0,
        N0_rep1=chainarray[:, 7],
        N0_rep2=chainarray[:, 8],
        V0_rep1=chainarray[:, 9],
        V0_rep2=chainarray[:, 10],
        σ=chainarray[:, 11])

    stats_chain = ess_rhat(chainfit)

    # Save Arrays as CSV
    CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/results/SALEXP_chains_system_final_basal_" * ostreo_list[i][1:4] * "_" * string(today()) * "_genotoul.csv", chain_df, writeheader=true)
    CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/results/SALEXP_statchains_system_final_basal_" * ostreo_list[i][1:4] * "_" * string(today()) * "_genotoul.csv", DataFrame(stats_chain), writeheader=true)

    println(ostreo_list[i][1:4] * " is done!")
end
