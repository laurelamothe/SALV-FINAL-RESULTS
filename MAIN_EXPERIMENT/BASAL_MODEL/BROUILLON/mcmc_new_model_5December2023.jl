##########################################################################
#################### Laure L. / Main Exp. SALINITY  ######################
############################## Modelling #################################
############### System - mcmc to estimate all parameters   ###############
################# David/Laure 5 December 2023 - Roscoff ##################
##########################################################################

#%% PACKAGES
using Distributed
addprocs(4)
@everywhere using Dates, Turing, Distributions, DifferentialEquations, Statistics, CSV, DataFrames, Glob

##########################################################################
##########################################################################
################            PRELIMINARY STEPS             ################
##########################################################################
##########################################################################
@everywhere begin
    # list of experimental data files
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
    virus_list = [
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

    # result of fitting procedure on control cultures
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

    #### IMPORT DATA COMMON TO EVERY CONDITIONS ####
    # Mesured virus decay when available and calculated from linear regression otherwise
    data_virus_decay = DataFrame(CSV.File("./VIRUS_DECAY/data/data_input_virus_decay.csv"))[:, 1:2]

    #Calculated e
    data_e = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/data/e_new_model_06December2023.csv"))

    #### CHARGE EXPERIMENTAL DATA OF ALL CONDITIONS ####
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
    df_ostreo_rep1 = []
    df_virus_rep1 = []
    df_time_ostreo = []
    df_ostreo_rep2 = []
    df_virus_rep2 = []
    df_time_virus = []

    # download all the data
    for i in 1:length(ostreo_list)
        push!(df_ostreo_rep1, ostreo_df(i)[1])
        push!(df_ostreo_rep2, ostreo_df(i)[2])
        push!(df_time_ostreo, ostreo_df(i)[3])

        push!(df_virus_rep1, virus_df(i)[1])
        push!(df_virus_rep2, virus_df(i)[2])
        push!(df_time_virus, virus_df(i)[3])
    end

    #### IMPORT INPUT SCRIPTS ####
    #model (system of differential equations)
    include("./ode_new_model_25November2023.jl")
    #function to obtain the estmated parameter on control cultures
    include("./moyenne_param_control_22May2023.jl")


    #### CREATE THE ODE PROBLEM ####
    #initial conditions of the model
    t0 = 0.0
    dt = 0.01
    u0 = [1E6, 0, 0, 1E5]
    tf = 14
    tmod = collect(t0:dt:tf)
    tspan = (t0, tf)
    p = [6.6, 0.01, 1E8, 0.5, 1e-8, 0.5, 50, 0.0001]
    #definition of the initial ODE problem
    prob_system = ODEProblem(ODE_system, u0, tspan, p)
end

##########################################################################
##########################################################################
################   DEFINITION OF CALIBRATION PROCEDURE    ################
##########################################################################
##########################################################################

@everywhere begin
    @model function fit_system(fixed_r, fixed_K, fixed_δ, fixed_e, data_ostreo_rep1, data_ostreo_rep2, data_virus_rep1, data_virus_rep2, time_ostreo, time_virus, prob_system, tmod, dt)
        #### SET PARAMETERS ####
        # fixed parameters
        r = fixed_r
        K = fixed_K
        e = fixed_e
        δ = fixed_δ

        # parameters to calibrate
        µ ~ truncated(LogNormal(log(0.08), 1), 0, 1)
        α ~ truncated(LogNormal(log(0.1), 1), 0, 1 / r) # upper limit set to stay under 1
        ϵ ~ truncated(LogNormal(log(0.7), 1), 0, 1)
        η ~ truncated(LogNormal(log(3), 1), 0.01, 20)
        β ~ truncated(LogNormal(log(50), 0.5), 5, 100)

        N0_rep1 ~ truncated(LogNormal(data_ostreo_rep1[1], 0.2), 5E5, 1E7)
        N0_rep2 ~ truncated(LogNormal(data_ostreo_rep2[1], 0.2), 5E5, 1E7)
        V0_rep1 ~ truncated(LogNormal(data_virus_rep1[1], 0.2), 5E5, 1E7)
        V0_rep2 ~ truncated(LogNormal(data_virus_rep2[1], 0.2), 5E5, 1E7)
        u0_rep1 = [N0_rep1, 0, 0, V0_rep1]
        u0_rep2 = [N0_rep2, 0, 0, V0_rep2]

        σ_N ~ InverseGamma(10, 1)
        σ_V ~ InverseGamma(10, 1)

        # parameters
        p_mcmc = [r, μ, K, e, α, ϵ, η, β, δ]

        ## INTERGRATION
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

        #### LIKELIHOOD ####
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

        #virus
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


##########################################################################
##########################################################################
################     CALIBRATION IN EVERY CONDITIONS      ################
##########################################################################
##########################################################################

#### SET CHAIN PARAMETERS ####
@everywhere begin
    nburn = 2000
    nstep = 2000
    nchain = 4
end


for i in 2:6 #length(ostreo_list)
    #### EXTRACT FIXED PARAMETERS AT CURRENT SALINITY ####
    fixed_r = get_mean_control(control_results).r[i]
    fixed_K = get_mean_control(control_results).K[i]
    fixed_δ = data_virus_decay.light[i]
    fixed_e = data_e.e[i]

    #### MCMC ####
    #update ODE problem
    modelfit = fit_system(fixed_r, fixed_K, fixed_δ, fixed_e, df_ostreo_rep1[i], df_ostreo_rep2[i], df_virus_rep1[i], df_virus_rep2[i], df_time_ostreo[i], df_time_virus[i], prob_system, tmod, dt)

    #run chains
    chainfit = reduce(chainscat, pmap(x -> sample(modelfit, NUTS(nburn, 0.65), nstep, save_state=false), 1:nchain))

    #export
    chainarray = Array(chainfit)
    chain_df = DataFrame(
        μ=chainarray[:, 1],
        α=chainarray[:, 2],
        ϵ=chainarray[:, 3],
        η=chainarray[:, 4],
        β=chainarray[:, 5],
        N0_rep1=chainarray[:, 6],
        N0_rep2=chainarray[:, 7],
        V0_rep1=chainarray[:, 8],
        V0_rep2=chainarray[:, 9],
        σ=chainarray[:, 10])
    stats_chain = ess_rhat(chainfit)
    CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/SALEXP_chains_new_model_" * ostreo_list[i][1:4] * "_" * string(today()) * "_genotoul.csv", chain_df, writeheader=true)
    CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/SALEXP_statchains_new_model_" * ostreo_list[i][1:4] * "_" * string(today()) * "_genotoul.csv", DataFrame(stats_chain), writeheader=true)

    println(ostreo_list[i][1:4] * " is done!")
end
