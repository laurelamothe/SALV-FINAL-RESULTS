########################################################
########## Laure L. / Main Exp. SALINITY  ##############
#################### Modelling #########################
########## make all plot and save in a file ############
########## David/Laure 03 July 2023 - Paris ############
########################################################

using StatsPlots, CSV, DataFrames, Plots, Measures, StatsBase, Distributions, Dates, DifferentialEquations
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
Rper_list = collect(0:0.15:0.7)
#PLOT DYNAMICS: initialisation
ω = 0
δ = 0.028
t0 = 0.0
dt = 0.01
tf = 14
tmod = collect(t0:dt:tf)
tspan = (t0, tf)
nrep = 2;
nstep = 10; #For the true plot you can have 4000 steps
populations = 5
sol_theo = []
mresol = NaN * ones(nstep, length(tmod), populations);
mquantile = NaN * ones(populations, length(tmod), 3, nrep);
mquantile_mean = NaN * ones(5, length(tmod), 3);
include("../BASAL_MODEL/ode_system_13July2023.jl")
include("../BASAL_MODEL/moyenne_param_control_22May2023.jl")
result_virus_decay = DataFrame(CSV.File("./VIRUS_DECAY/results/SALEXP_virus_decay_light_2023-11-05.csv"))
for i in 7:9#length(result_list_system)
    pl = plot(layout=(5, 1), plot_title="from " * result_list_control[i][23:26], size=(500, 1500), margins=7mm, dpi=600, yscale=:log10, ylims=(1e1, 1e10))

    data_ostreo = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/data/" * ostreo_list[i]))
    data_virus = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/data/" * virus_list[i]))

    scatter!(pl, subplot=5, data_ostreo.TIME, (data_ostreo.REPLICATE_1 .+ data_ostreo.REPLICATE_2) ./ 2, color=:gray, alpa=0.3, markerstrokewidth=0, markersize=3)
    scatter!(pl, subplot=4, data_virus.TIME, (data_virus.REPLICATE_1 .+ data_virus.REPLICATE_2) ./ 2, color=:gray, alpa=0.3, markerstrokewidth=0, markersize=3)

    for R in 1:length(Rper_list)
        #import data
        data_system = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/" * result_list_system[i]))
        data_system.ϕ = exp.(data_system.ϕ)
        data_system = hcat(DataFrame(
                Kr=data_system.Kr,
                µ=data_system.µ
            ), data_system[:, 3:11])
        data_control = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/control/" * result_list_control[i]))

        ### 4: DYNAMICS ###
        r = get_mean_control(result_list_control).r[i]
        Ks = get_mean_control(result_list_control).K[i]
        δ = mean(result_virus_decay.a) * salinity[i] + mean(result_virus_decay.b)
        dt_sols = []

        for j = 1:1#nrep
            for k in 1:nstep
                id = rand((1:size(data_system)[1]))
                #defined parameter spcific to the condition
                µ = data_system.µ[id]
                Kr = data_system.Kr[id]
                ϵ = data_system.ϵ[id]
                ϕ = data_system.ϕ[id]
                η = data_system.η[id]
                β = data_system.β[id]
                Rper = Rper_list[R]
                if j == 1
                    N0 = data_system.N0_rep1[id]
                    V0 = data_system.V0_rep1[id]
                else
                    N0 = data_system.N0_rep2[id]
                    V0 = data_system.V0_rep2[id]
                end

                u0 = [(1 - Rper) * N0, Rper * N0, 0, V0]
                p = [r, μ, Ks, Kr, ω, ϵ, ϕ, η, β, δ]

                #solve system
                prob_system = ODEProblem(ODE_system, u0, tspan, p)
                sol_theo = solve(prob_system, AutoTsit5(Tsit5()), abstol=1E-8, reltol=1E-8, saveat=dt)
                sol_theo = DataFrame(S=sol_theo[1, :], R=sol_theo[2, :], I=sol_theo[3, :], V=sol_theo[4, :], N=sol_theo[1, :] .+ sol_theo[2, :] .+ sol_theo[3, :])
                for l in 1:populations
                    mresol[k, :, l] = sol_theo[:, l]
                end
            end


            for l in 1:populations
                #calculate 5% and 95% quantiles
                for k in 1:length(tmod)
                    mquantile[l, k, :, j] = abs.(quantile(Array(mresol[:, k, l]), [0.025, 0.5, 0.975]))
                    mquantile[l, k, 2, j] = mean(mresol[:, k, l])
                end
                plot!(pl, subplot=l, tmod, mquantile[l, :, 2, j], color=color_palette[R], label="Rper = " * string(Rper_list[R] * 100) * "%", legend=false)

            end
        end
        println(ostreo_list[i][1:4] * " IS DONE !")
    end
    plot!(pl, subplot=1, legend=true)
    savefig(pl, "./MAIN_EXPERIMENT/THEORICAL ANALYSIS RPER/figures/plot_resistance_théo_" * ostreo_list[i][1:4] * "_" * string(today()) * ".png")

end

pl