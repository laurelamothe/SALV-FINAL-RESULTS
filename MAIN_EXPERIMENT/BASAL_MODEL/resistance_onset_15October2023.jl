########################################################
########## Laure L. / Main Exp. SALINITY  ##############
#################### Modelling #########################
######## analyse the emergence of resistance ###########
####### David/Laure 12 October 2023 - Roscoff ##########
########################################################

using CSV, DataFrames, Measures, StatsBase, Dates, Distributions, DifferentialEquations

data_list_host = [
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
    "D_S6_host_infected.csv"
]
result_list_control = [
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
    "SALEXP_chains_control_D_S6_2023-07-12.csv"
]

result_list_system = [
    "SALEXP_chains_system_final_basal_A_S1_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S2_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S3_2023-07-29_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S4_2023-07-29_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S5_2023-07-29_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S6_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S1_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S2_2023-07-29_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S3_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S4_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S5_2023-07-29_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S6_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S1_2023-07-29_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S2_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S3_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S4_2023-07-29_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S5_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S6_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S1_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S2_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S3_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S4_2023-07-29_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S5_2023-07-30_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S6_2023-07-30_genotoul.csv"
]

color_palette = ["#72DB23", "#12B4E0", "#DF321B", "#EB8203", "#000000"]
salinity = [5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40]

#PLOT DYNAMICS: initialisation
ω = 0
δ = 0.028
t0 = 0.0
dt = 0.01
tf = 14
tmod = collect(t0:dt:tf)
tspan = (t0, tf)
nrep = 2;
nstep = 8000; #For the true plot you can have 4000 steps
populations = 5
sol_theo = []
mresol = NaN * ones(nstep, length(tmod), populations);
mquantile = NaN * ones(populations, length(tmod), 3, nrep);
mquantile_mean = NaN * ones(5, length(tmod));
onset_R = DataFrame(emerg=ones(nstep), pop=ones(nstep))
include("./ode_system_13July2023.jl")
include("./moyenne_param_control_22May2023.jl")
result_virus_decay = DataFrame(CSV.File("./VIRUS_DECAY/results/SALEXP_virus_decay_light_2023-07-12.csv"))



for i in 1:length(result_list_system)

    #import data
    data_system = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/" * result_list_system[i]))
    data_system.ϕ = exp.(data_system.ϕ)
    data_system = hcat(DataFrame(
            Kr=data_system.Kr,
            µ=data_system.µ
        ), data_system[:, 3:11])
    data_control = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/control/" * result_list_control[i]))

    ##DYNAMICS
    r = get_mean_control(result_list_control).r[i]
    Ks = get_mean_control(result_list_control).K[i]
    δ = mean(result_virus_decay.a) * salinity[i] + mean(result_virus_decay.b)
    dt_sols = []

    for k in 1:nstep
        id = rand((1:size(data_system)[1]))
        #defined parameter spcific to the condition
        µ = data_system.µ[id]
        Kr = data_system.Kr[id]
        ϵ = data_system.ϵ[id]
        ϕ = data_system.ϕ[id]
        η = data_system.η[id]
        β = data_system.β[id]
        Rper = data_system.Rper[id]

        N0 = (data_system.N0_rep1[id] + data_system.N0_rep2[id]) / 2
        V0 = (data_system.V0_rep1[id] + data_system.V0_rep2[id]) / 2

        u0 = [(1 - Rper) * N0, Rper * N0, 0, V0]
        p = [r, μ, Ks, Kr, ω, ϵ, ϕ, η, β, δ]

        #solve system
        prob_system = ODEProblem(ODE_system, u0, tspan, p)
        sol_theo = solve(prob_system, AutoTsit5(Tsit5()), abstol=1E-8, reltol=1E-8, saveat=dt)
        sol_theo = DataFrame(S=sol_theo[1, :], R=sol_theo[2, :], I=sol_theo[3, :], V=sol_theo[4, :], N=sol_theo[1, :] .+ sol_theo[2, :] .+ sol_theo[3, :])
        for l in 1:populations
            mresol[k, :, l] = sol_theo[:, l]
        end
        #find the time when R reach 1% of N
        for o in 1:length(tmod)
            if mresol[k, o, 2] < last(mresol[1, :, 5]) / 100
                onset_R.emerg[k] = tmod[o]
            end
            if mresol[k, o, 2] > last(mresol[1, :, 5]) / 2
                onset_R.pop[k] = tmod[o] - onset_R.emerg[k]
                break
            end
        end
    end

    println(data_list_host[i][1:4] * " IS DONE !")
    CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/results/SALEXP_R_onset_" * data_list_host[i][1:4] * "_" * string(today()) * ".csv", onset_R, writeheader=true)

end