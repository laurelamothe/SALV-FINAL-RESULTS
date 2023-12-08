########################################################
########## Laure L. / Main Exp. SALINITY  ##############
#################### Modelling #########################
########## make all plot and save in a file ############
########## David/Laure 03 July 2023 - Paris ############
########################################################
aie = 1
#list of plot : 
#   chains
#   rhats
#   corr
#   dynamics
#   prior posterior
using StatsPlots, CSV, DataFrames, Plots, Measures, StatsBase, Distributions, Dates, DifferentialEquations

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

data_list_virus = [
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
    "D_S6_virus.csv"
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
    "SALEXP_chains_system_final_basal_D_S6_2023-07-30_genotoul.csv"]

stat_result_list_system = [
    "SALEXP_statchains_system_final_basal_A_S1_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_A_S2_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_A_S3_2023-07-29_genotoul.csv",
    "SALEXP_statchains_system_final_basal_A_S4_2023-07-29_genotoul.csv",
    "SALEXP_statchains_system_final_basal_A_S5_2023-07-29_genotoul.csv",
    "SALEXP_statchains_system_final_basal_A_S6_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_B_S1_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_B_S2_2023-07-29_genotoul.csv",
    "SALEXP_statchains_system_final_basal_B_S3_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_B_S4_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_B_S5_2023-07-29_genotoul.csv",
    "SALEXP_statchains_system_final_basal_B_S6_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_C_S1_2023-07-29_genotoul.csv",
    "SALEXP_statchains_system_final_basal_C_S2_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_C_S3_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_C_S4_2023-07-29_genotoul.csv",
    "SALEXP_statchains_system_final_basal_C_S5_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_C_S6_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_D_S1_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_D_S2_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_D_S3_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_D_S4_2023-07-29_genotoul.csv",
    "SALEXP_statchains_system_final_basal_D_S5_2023-07-30_genotoul.csv",
    "SALEXP_statchains_system_final_basal_D_S6_2023-07-30_genotoul.csv"]

result_list_resistance = [
    "SALEXP_R_onset_A_S1_2023-10-15.csv",
    "SALEXP_R_onset_A_S2_2023-10-15.csv",
    "SALEXP_R_onset_A_S3_2023-10-15.csv",
    "SALEXP_R_onset_A_S4_2023-10-15.csv",
    "SALEXP_R_onset_A_S5_2023-10-15.csv",
    "SALEXP_R_onset_A_S6_2023-10-15.csv",
    "SALEXP_R_onset_B_S1_2023-10-15.csv",
    "SALEXP_R_onset_B_S2_2023-10-15.csv",
    "SALEXP_R_onset_B_S3_2023-10-15.csv",
    "SALEXP_R_onset_B_S4_2023-10-15.csv",
    "SALEXP_R_onset_B_S5_2023-10-15.csv",
    "SALEXP_R_onset_B_S6_2023-10-15.csv",
    "SALEXP_R_onset_C_S1_2023-10-15.csv",
    "SALEXP_R_onset_C_S2_2023-10-15.csv",
    "SALEXP_R_onset_C_S3_2023-10-15.csv",
    "SALEXP_R_onset_C_S4_2023-10-15.csv",
    "SALEXP_R_onset_C_S5_2023-10-15.csv",
    "SALEXP_R_onset_C_S6_2023-10-15.csv",
    "SALEXP_R_onset_D_S1_2023-10-15.csv",
    "SALEXP_R_onset_D_S2_2023-10-15.csv",
    "SALEXP_R_onset_D_S3_2023-10-15.csv",
    "SALEXP_R_onset_D_S4_2023-10-15.csv",
    "SALEXP_R_onset_D_S5_2023-10-15.csv",
    "SALEXP_R_onset_D_S6_2023-10-15.csv",
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
nstep = 1000; #For the true plot you can have 4000 steps
populations = 5
sol_theo = []
mresol = NaN * ones(nstep, length(tmod), populations);
mquantile = NaN * ones(populations, length(tmod), 3, nrep);
mquantile_mean = NaN * ones(5, length(tmod), 3);
include("./ode_system_13July2023.jl")
include("./moyenne_param_control_22May2023.jl")
result_virus_decay = DataFrame(CSV.File("./VIRUS_DECAY/results/SALEXP_virus_decay_light_2023-11-05.csv"))

pchain = plot(layout=(2, 4), size=(800, 400), margin=0mm, dpi=600)
prhat = plot(size=(200, 200), margin=2mm, dpi=600)
pcor = plot(layout=(6, 6), size=(500, 500), margin=0mm, left_margin=1mm, xscale=:log, yscale=:log, dpi=600)
pdyn = plot(layout=(2, 4), size=(1100, 500), margin=0mm, yscale=:log10, dpi=600)
p_all = plot(layout=(4, 6), size=(1200, 800), margins=0mm, dpi=600, yscale=:log10, ylims=(1e3, 1e10), xlims=(0, 14), tickfontcolor=:white, legend=false, grid=true)



for i in 1:length(result_list_system)
    #create new directory
    dir_path = "./MAIN_EXPERIMENT/BASAL_MODEL/figures/" * string(today()) * ""
    mkpath(dir_path)
    mkpath(dir_path * "/chains")
    mkpath(dir_path * "/rhats")
    mkpath(dir_path * "/covariances")
    mkpath(dir_path * "/dynamics")

    #import data
    data_system = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/" * result_list_system[i]))
    data_system.ϕ = exp.(data_system.ϕ)
    data_system = hcat(DataFrame(
            Kr=data_system.Kr,
            µ=data_system.µ
        ), data_system[:, 3:11])
    data_control = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/control/" * result_list_control[i]))
    stat_system = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/" * stat_result_list_system[i]))
    exp_points_host = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/data/" * data_list_host[i]))
    exp_points_host = exp_points_host[completecases(exp_points_host), :]
    exp_points_virus = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/data/" * data_list_virus[i]))
    exp_points_virus = exp_points_virus[completecases(exp_points_virus), :]
    #initialize Plots

    pchain = plot(layout=(2, 4), size=(800, 400), margin=0mm, dpi=600)
    prhat = plot(size=(200, 200), margin=2mm, dpi=600)
    pcor = plot(layout=(6, 6), size=(500, 500), margin=0mm, left_margin=1mm, xscale=:log, yscale=:log, dpi=600)
    pdyn = plot(layout=(2, 5), size=(1100, 500), margin=0mm, yscale=:log10, dpi=600)

    if aie == 1
        ### 1: CHAINS ###
        for j in 1:8
            plot!(pchain, subplot=j, data_system[:, j], title=names(data_system)[j])
        end
        plot!(pchain, legend=false)
        savefig(pchain, dir_path * "/chains/" * data_list_host[i][1:4] * ".png")

        ### 2: RHATS ###
        scatter!(prhat, stat_system.parameters[1:6], stat_system.rhat[1:6], legend=false, ylims=(0.95, 1.05))
        savefig(prhat, dir_path * "/rhats/" * data_list_host[i][1:4] * ".png")
        ### 3: CORR ###
        index = rand(1:2000, 500)

        for i in 1:6
            for j in 1:6
                scatter!(pcor, subplot=j + 6 * (i - 1),
                    data_system[index, j], data_system[index, i],
                    markeralpha=0.2,
                    markersize=2.5,
                    markercolor=:black)
                plot!(ticks=false, legend=:none)
                if i == 1
                    plot!(pcor, subplot=j, title=names(data_system)[j])
                end
                if j == 1
                    plot!(pcor, subplot=1 + 6 * (i - 1), ylabel=names(data_system)[i])
                end
            end
        end
        savefig(pcor, dir_path * "/covariances/" * data_list_host[i][1:4] * ".png")

        ### 4: DYNAMICS ###
        r = get_mean_control(result_list_control).r[i]
        Ks = get_mean_control(result_list_control).K[i]
        δ = mean(result_virus_decay.a) * salinity[i] + mean(result_virus_decay.b)
        dt_sols = []

        for j = 1:nrep
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
                #Plots
                plot!(pdyn, subplot=l + populations * (j - 1), tmod, mquantile[l, :, 1, j], fillrange=mquantile[l, :, 3, j], fillalpha=0.4,
                    color=color_palette[l], linealpha=0, legend=false, title=names(sol_theo)[l])

                plot!(pdyn, subplot=l + populations * (j - 1), tmod, mquantile[l, :, 2, j], lw=3, color=color_palette[l], legend=false)
            end
            scatter!(pdyn, subplot=4, exp_points_virus.TIME, exp_points_virus.REPLICATE_1, color=color_palette[4], markerstrokewidth=0)
            scatter!(pdyn, subplot=9, exp_points_virus.TIME, exp_points_virus.REPLICATE_2, color=color_palette[4], markerstrokewidth=0)

            scatter!(pdyn, subplot=5, exp_points_host.TIME, exp_points_host.REPLICATE_1, color=color_palette[5], markerstrokewidth=0)
            scatter!(pdyn, subplot=10, exp_points_host.TIME, exp_points_host.REPLICATE_2, color=color_palette[5], markerstrokewidth=0)

            plot!(pdyn, ylims=(1e3, 1e10))
        end
        savefig(pdyn, dir_path * "/dynamics/" * data_list_host[i][1:4] * ".png")
    end
    ## ALL DYNAMICS
    """
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
    end

    #calculate 5% and 95% quantiles
    for k in 1:length(tmod)
        mquantile_mean[1, k, :] = abs.(quantile(Array(mresol[:, k, 5]), [0.025, 0.5, 0.975]))
        mquantile_mean[1, k, 2] = mean(mresol[:, k, 5])

        mquantile_mean[2, k, :] = abs.(quantile(Array(mresol[:, k, 4]), [0.025, 0.5, 0.975]))
        mquantile_mean[2, k, 2] = mean(mresol[:, k, 4])

        mquantile_mean[3, k, :] = abs.(quantile(Array(mresol[:, k, 1]), [0.025, 0.5, 0.975]))
        mquantile_mean[3, k, 2] = mean(mresol[:, k, 1])

        mquantile_mean[4, k, :] = abs.(quantile(Array(mresol[:, k, 2]), [0.025, 0.5, 0.975]))
        mquantile_mean[4, k, 2] = mean(mresol[:, k, 2])

        mquantile_mean[5, k, :] = abs.(quantile(Array(mresol[:, k, 3]), [0.025, 0.5, 0.975]))
        mquantile_mean[5, k, 2] = mean(mresol[:, k, 3])

    end
    #Plots

    temp_S = DataFrame(Time=tmod, Pop=mquantile_mean[3, :, 2])
    temp_S = filter(row -> row.Time in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], temp_S)
    temp_R = DataFrame(Time=tmod, Pop=mquantile_mean[4, :, 2])
    temp_R = filter(row -> row.Time in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], temp_R)
    temp_I = DataFrame(Time=tmod, Pop=mquantile_mean[5, :, 2])
    temp_I = filter(row -> row.Time in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], temp_I)


    bar!(p_all, subplot=i, temp_R.Time, (temp_R.Pop .+ temp_I.Pop .+ temp_S.Pop) ./ (temp_R.Pop .+ temp_I.Pop .+ temp_S.Pop), fillrange=temp_R.Pop .+ temp_I.Pop .+ temp_S.Pop .+ 1, color=color_palette[1], linecolor=nothing)
    bar!(p_all, subplot=i, temp_R.Time, (temp_R.Pop .+ temp_I.Pop) ./ (temp_R.Pop .+ temp_I.Pop), fillrange=temp_R.Pop .+ temp_I.Pop .+ 1, color=color_palette[3], linecolor=nothing)
    bar!(p_all, subplot=i, temp_R.Time, (temp_R.Pop) ./ (temp_R.Pop), fillrange=temp_R.Pop .+ 1, color=color_palette[2], linecolor=nothing)


    #plot!(p_all, subplot=i, temp_R.Time, (temp_R.Pop .+ temp_I.Pop .+ temp_S.Pop) ./ (temp_R.Pop .+ temp_I.Pop .+ temp_S.Pop), fillrange=temp_R.Pop .+ temp_I.Pop .+ temp_S.Pop, color="#6CF0A1", linecolor=nothing)
    #plot!(p_all, subplot=i, temp_R.Time, (temp_R.Pop .+ temp_I.Pop) ./ (temp_R.Pop .+ temp_I.Pop), fillrange=temp_R.Pop .+ temp_I.Pop, color="#95B5F5", linecolor=nothing)
    #plot!(p_all, subplot=i, temp_R.Time, (temp_I.Pop) ./ (temp_I.Pop), fillrange=temp_I.Pop, color="#F6A78E", linecolor=nothing)



    plot!(p_all, subplot=i, tmod, mquantile_mean[1, :, 1], fillrange=mquantile_mean[1, :, 3], fillalpha=0.4,
        color=color_palette[5], linealpha=0, legend=false)
    plot!(p_all, subplot=i, tmod, mquantile_mean[1, :, 2], lw=2, color=color_palette[5], legend=false)

    plot!(p_all, subplot=i, tmod, mquantile_mean[2, :, 1], fillrange=mquantile_mean[2, :, 3], fillalpha=0.4,
        color=color_palette[4], linealpha=0, legend=false)

    plot!(p_all, subplot=i, tmod, mquantile_mean[2, :, 2], lw=2, color=color_palette[4], legend=false)

    scatter!(p_all, subplot=i, exp_points_virus.TIME, (exp_points_virus.REPLICATE_1 + exp_points_virus.REPLICATE_2) / 2,
        yerr=abs.(exp_points_virus.REPLICATE_1 .- exp_points_virus.REPLICATE_2), color=color_palette[4], markerstrokecolor=color_palette[4], markeralpha=0, markerstrokewidth=1.5)
    scatter!(p_all, subplot=i, exp_points_virus.TIME, (exp_points_virus.REPLICATE_1 + exp_points_virus.REPLICATE_2) / 2, color=color_palette[4], markerstrokewidth=0)

    scatter!(p_all, subplot=i, exp_points_host.TIME, (exp_points_host.REPLICATE_1 + exp_points_host.REPLICATE_2) / 2,
        yerr=abs.(exp_points_host.REPLICATE_1 .- exp_points_host.REPLICATE_2), color=color_palette[3], markerstrokecolor=color_palette[5], markeralpha=0, markerstrokewidth=1.5)
    scatter!(p_all, subplot=i, exp_points_host.TIME, (exp_points_host.REPLICATE_1 + exp_points_host.REPLICATE_2) / 2, color=color_palette[5], markerstrokewidth=0)

    if i < 7
        plot!(p_all, subplot=i, title=string(salinity[i]) * " PSU")
    end
    if i in [1, 7, 13, 19]
        plot!(p_all, subplot=i, ytickfontcolor=:black, ylabel="particle.mL⁻¹")
    end
    if 18 < i
        plot!(p_all, subplot=i, xtickfontcolor=:black, xlabel="day")
    end
    """
    println(data_list_host[i][1:4] * " IS DONE !")
end
plot!(p_all, subplot=1, left_margin=5mm)
plot!(p_all, subplot=19, bottom_margin=5mm)
plot!(axis=true, xticks=true, tickfontsize=9, axis_color=:black, legend=false, grid=true)

#savefig(p_all, "./MAIN_EXPERIMENT/BASAL_MODEL/figures/all_conditions/dynamics_all_fit_system_" * string(today()) * ".png")

pdyn
