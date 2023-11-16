########################################################
########## Laure L. / Main Exp. SALINITY  ##############
#################### Modelling #########################
# Control - plotting mcmc result according to salinity #
########### David/Laure 03 May 2023 - Banyuls ##########
########################################################

using StatsPlots, CSV, DataFrames, Measures, StatsBase, Plots, Dates


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
    "SALEXP_chains_system_final_basal_A_S1_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S2_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S3_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S4_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S5_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_A_S6_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S1_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S2_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S3_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S4_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S5_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_B_S6_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S1_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S2_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S3_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S4_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S5_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_C_S6_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S1_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S2_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S3_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S4_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S5_2023-10-21_genotoul.csv",
    "SALEXP_chains_system_final_basal_D_S6_2023-10-21_genotoul.csv"
]


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

result_list_hyper = [
    ["SALEXP_chains_hyperfunction_polynomial_r_A_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_polynomial_r_B_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_polynomial_r_C_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_polynomial_r_D_2023-11-03.csv"],
    ["SALEXP_chains_hyperfunction_constant_µ_A_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_µ_B_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_µ_C_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_µ_D_2023-11-03.csv"],
    ["SALEXP_chains_hyperfunction_normal_K_A_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_normal_K_B_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_normal_K_C_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_normal_K_D_2023-11-03.csv"],
    ["SALEXP_chains_hyperfunction_polynomial_Kr_A_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_polynomial_Kr_B_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_polynomial_Kr_C_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_polynomial_Kr_D_2023-11-03.csv"],
    ["SALEXP_chains_hyperfunction_constant_ϵ_A_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_constant_ϵ_B_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_constant_ϵ_C_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_constant_ϵ_D_2023-11-04.csv"],
    ["SALEXP_chains_hyperfunction_linear_ϕ_A_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_linear_ϕ_B_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_linear_ϕ_C_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_linear_ϕ_D_2023-11-04.csv"],
    ["SALEXP_chains_hyperfunction_normal_η_A_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_normal_η_B_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_normal_η_C_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_normal_η_D_2023-11-03.csv"],
    ["SALEXP_chains_hyperfunction_constant_β_A_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_β_B_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_β_C_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_β_D_2023-11-03.csv"],
    ["SALEXP_chains_hyperfunction_constant_ρ_A_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_ρ_B_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_ρ_C_2023-11-03.csv"
        "SALEXP_chains_hyperfunction_constant_ρ_D_2023-11-03.csv"],
    ["SALEXP_chains_hyperfunction_constant_Δ_A_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_constant_Δ_B_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_constant_Δ_C_2023-11-04.csv"
        "SALEXP_chains_hyperfunction_constant_Δ_D_2023-11-04.csv"]
]

name_salinity = ["S1", "S2", "S3", "S4", "S5", "S6"];
name_strain = ["RCC4221", "RCC1116", "RCC1123", "RCC1108"];
salinities = [5, 10, 15, 25, 35, 40]
p_list = ["r", "µ", "K", "Kr", "ϵ", "ϕ", "η", "β", "ρ", "Δ"]
color_palette = ["#F79135", "#70AD47", "#4DA2E1", "#A06DDD"]
include("./selected_hyperfunction.jl")

#initialze df in which we will store all the data
df = DataFrame(Salinity=[], Strain=[], r=[], μ=[], K=[], Kr=[], ϵ=[], ϕ=[], η=[], β=[], ρ=[], Δ=[])
means = DataFrame(Salinity=[], Strain=[], r=[], μ=[], K=[], Kr=[], ϵ=[], ϕ=[], η=[], β=[], ρ=[], Δ=[])
bars = DataFrame(Salinity=[], Strain=[], r=[], μ=[], K=[], Kr=[], ϵ=[], ϕ=[], η=[], β=[], ρ=[], Δ=[])
pl = Plots.plot(layout=(2, 5), margins=0mm, xtickfontcolor=:white)
sal = 1;
strain = 1;
for i in 1:length(result_list_system)
    if result_list_system[i][22:25] == "X_S4"
        sal += 1
    else

        println(result_list_system[i][36:37])
        #import data
        data_system = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/" * result_list_system[i]))
        data_system.ϕ = exp.(data_system.ϕ)

        data_control = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/control/" * result_list_control[i]))
        data_resistance = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/resistance/" * result_list_resistance[i]))

        df = vcat(df, hcat(
            DataFrame(
                Salinity=repeat([salinities[sal]], size(data_control)[1]),
                Strain=repeat([result_list_system[i][34]], size(data_control)[1]),
                r=data_control[:, 1],
                µ=data_system[:, 2],
                K=data_control[:, 2],
                Kr=data_system[:, 1]),
            hcat(data_system[:, 3:6],
                DataFrame(ρ=data_resistance[!, 1], Δ=data_resistance[!, 2]))))

        means = vcat(means, DataFrame(
            Salinity=[salinities[sal]],
            Strain=[result_list_system[i][34]],
            r=[mean(df.r)],
            μ=[mean(df.µ)],
            K=[mean(df.K)],
            Kr=[mean(df.Kr)],
            ϵ=[mean(df.ϵ)],
            ϕ=[mean(df.ϕ)],
            η=[mean(df.η)],
            β=[mean(df.β)],
            ρ=[mean(df.ρ)],
            Δ=[mean(df.Δ)]))

        bars = vcat(bars, DataFrame(
            Salinity=[salinities[sal]],
            Strain=[result_list_system[i][34]],
            r=[abs.(quantile(df.r, 0.75) .- quantile(df.r, 0.25)) / 2],
            μ=[abs.(quantile(df.µ, 0.75) .- quantile(df.µ, 0.25)) / 2],
            K=[abs.(quantile(df.K, 0.75) .- quantile(df.K, 0.25)) / 2],
            Kr=[abs.(quantile(df.Kr, 0.75) .- quantile(df.Kr, 0.25)) / 2],
            ϵ=[abs.(quantile(df.ϵ, 0.75) .- quantile(df.ϵ, 0.25)) / 2],
            ϕ=[abs.(quantile(df.ϕ, 0.75) .- quantile(df.ϕ, 0.25)) / 2],
            η=[abs.(quantile(df.η, 0.75) .- quantile(df.η, 0.25)) / 2],
            β=[abs.(quantile(df.β, 0.75) .- quantile(df.β, 0.25)) / 2],
            ρ=[abs.(quantile(df.ρ, 0.75) .- quantile(df.ρ, 0.25)) / 2],
            Δ=[abs.(quantile(df.Δ, 0.75) .- quantile(df.Δ, 0.25)) / 2]))

        df = DataFrame(Salinity=[], Strain=[], r=[], μ=[], K=[], Kr=[], ϵ=[], ϕ=[], η=[], β=[], ρ=[], Δ=[])

        #count the current strain and salinity we are at
        sal += 1

        if result_list_system[i][36:37] == "S6"

            for p in 1:length(p_list)

                if p == 2 || p == 5
                    scatter!(subplot=p,
                        means.Salinity, means[:, p+2] .* 100,
                        yerr=bars[:, p+2] .* 100,
                        label=false,
                        markerstrokecolor=color_palette[strain], markercolor=color_palette[strain], markeralpha=0.3
                    )
                    scatter!(subplot=p,
                        means.Salinity, means[:, p+2] .* 100,
                        label=name_strain[strain],
                        legend=false,
                        color=color_palette[strain],
                        markerstrokewidth=0,
                        markersize=7,
                        title=p_list[p],
                        titlefontvalign=:center
                    )
                else
                    scatter!(subplot=p,
                        means.Salinity, means[:, p+2],
                        yerr=bars[:, p+2],
                        label=false,
                        markerstrokecolor=color_palette[strain], markercolor=color_palette[strain], markeralpha=0.3
                    )
                    scatter!(subplot=p,
                        means.Salinity, means[:, p+2],
                        label=name_strain[strain],
                        legend=false,
                        color=color_palette[strain],
                        markerstrokewidth=0,
                        markersize=7,
                        title=p_list[p],
                        titlefontvalign=:center
                    )

                end


                ## PLOT HYPERFUNCTIONS
                if p_list[p] == "r"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    a = data_hyper.a
                    b = data_hyper.b
                    c = data_hyper.c
                    pol_mean = zeros(2, length(samp))
                    pol_025 = zeros(2, length(samp))
                    pol_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        pol = poly(a, b, c, [samp[t]])
                        pol_mean[:, t] = [samp[t], mean(pol)]
                        pol_025[:, t] = [samp[t], quantile(pol, 0.025)]
                        pol_975[:, t] = [samp[t], quantile(pol, 0.975)]
                    end
                    plot!(pl, subplot=p, pol_mean[1, :], pol_mean[2, :], color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, pol_025[1, :], pol_025[2, :], fillrange=pol_975[2, :], fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="cell.mL⁻¹.day⁻¹")

                end
                if p_list[p] == "η"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    a = data_hyper.a
                    σ = data_hyper.σm
                    µ = data_hyper[:, 3]
                    pol_mean = zeros(2, length(samp))
                    pol_025 = zeros(2, length(samp))
                    pol_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        pol = normal(a, σ, µ, [samp[t]])
                        pol_mean[:, t] = [samp[t], mean(pol)]
                        pol_025[:, t] = [samp[t], quantile(pol, 0.025)]
                        pol_975[:, t] = [samp[t], quantile(pol, 0.975)]
                    end
                    plot!(pl, subplot=p, pol_mean[1, :], pol_mean[2, :], color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, pol_025[1, :], pol_025[2, :], fillrange=pol_975[2, :], fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="day⁻¹")

                end

                if p_list[p] == "ϕ"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    a = data_hyper.a
                    b = data_hyper.b
                    logi_mean = zeros(2, length(samp))
                    logi_025 = zeros(2, length(samp))
                    logi_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        logi = linear(a, b, [samp[t]])
                        logi_mean[:, t] = [samp[t], mean(logi)]
                        logi_025[:, t] = [samp[t], quantile(logi, 0.025)]
                        logi_975[:, t] = [samp[t], quantile(logi, 0.975)]
                    end
                    plot!(pl, subplot=p, logi_mean[1, :], logi_mean[2, :], color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, logi_025[1, :], logi_025[2, :], fillrange=logi_975[2, :], fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="mL.day⁻¹")

                end

                if p_list[p] == "K"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    a = data_hyper.a
                    σ = data_hyper.σm
                    µ = data_hyper[:, 3]
                    logi_mean = zeros(2, length(samp))
                    logi_025 = zeros(2, length(samp))
                    logi_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        logi = normal(a, σ, µ, [samp[t]]) .* 1e7
                        logi_mean[:, t] = [samp[t], mean(logi)]
                        logi_025[:, t] = [samp[t], quantile(logi, 0.025)]
                        logi_975[:, t] = [samp[t], quantile(logi, 0.975)]
                    end
                    plot!(pl, subplot=p, logi_mean[1, :], logi_mean[2, :], color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, logi_025[1, :], logi_025[2, :], fillrange=logi_975[2, :], fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="cell.mL⁻¹")
                end

                if p_list[p] == "Kr"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    a = data_hyper.a
                    b = data_hyper.b
                    c = data_hyper.c
                    lin_mean = zeros(2, length(samp))
                    lin_025 = zeros(2, length(samp))
                    lin_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        lin = poly(a, b, c, [samp[t]]) .* 1e7
                        lin_mean[:, t] = [samp[t], mean(lin)]
                        lin_025[:, t] = [samp[t], quantile(lin, 0.025)]
                        lin_975[:, t] = [samp[t], quantile(lin, 0.975)]
                    end
                    plot!(pl, subplot=p, lin_mean[1, :], lin_mean[2, :], color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, lin_025[1, :], lin_025[2, :], fillrange=lin_975[2, :], fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="cell.mL⁻¹")

                end

                if p_list[p] == "µ"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    c = data_hyper.c
                    cons_mean = zeros(2, length(samp))
                    cons_025 = zeros(2, length(samp))
                    cons_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        cons = constant(c, [samp[t]])
                        cons_mean[:, t] = [samp[t], mean(cons)]
                        cons_025[:, t] = [samp[t], quantile(cons, 0.025)]
                        cons_975[:, t] = [samp[t], quantile(cons, 0.975)]
                    end
                    plot!(pl, subplot=p, cons_mean[1, :], cons_mean[2, :] .* 100, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, cons_025[1, :], cons_025[2, :] .* 100, fillrange=cons_975[2, :] .* 100, fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="%")

                end

                if p_list[p] == "ϵ"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    c = data_hyper.c
                    #b = data_hyper.b
                    lin_mean = zeros(2, length(samp))
                    lin_025 = zeros(2, length(samp))
                    lin_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        lin = constant(c, [samp[t]])
                        lin_mean[:, t] = [samp[t], mean(lin)]
                        lin_025[:, t] = [samp[t], quantile(lin, 0.025)]
                        lin_975[:, t] = [samp[t], quantile(lin, 0.975)]
                    end
                    plot!(pl, subplot=p, lin_mean[1, :], lin_mean[2, :] .* 100, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, lin_025[1, :], lin_025[2, :] .* 100, fillrange=lin_975[2, :] .* 100, fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="%")

                end
                if p_list[p] == "β"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    c = data_hyper.c
                    cons_mean = zeros(2, length(samp))
                    cons_025 = zeros(2, length(samp))
                    cons_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        cons = constant(c, [samp[t]])
                        cons_mean[:, t] = [samp[t], mean(cons)]
                        cons_025[:, t] = [samp[t], quantile(cons, 0.025)]
                        cons_975[:, t] = [samp[t], quantile(cons, 0.975)]
                    end
                    plot!(pl, subplot=p, cons_mean[1, :], cons_mean[2, :], color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, cons_025[1, :], cons_025[2, :], fillrange=cons_975[2, :], fillalpha=0.2, lw=0, color=color_palette[strain], label=false)

                end

                if p_list[p] == "ρ"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    c = data_hyper.c
                    cons_mean = zeros(2, length(samp))
                    cons_025 = zeros(2, length(samp))
                    cons_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        cons = constant(c, [samp[t]])
                        cons_mean[:, t] = [samp[t], mean(cons)]
                        cons_025[:, t] = [samp[t], quantile(cons, 0.025)]
                        cons_975[:, t] = [samp[t], quantile(cons, 0.975)]
                    end
                    plot!(pl, subplot=p, cons_mean[1, :], cons_mean[2, :], color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, cons_025[1, :], cons_025[2, :], fillrange=cons_975[2, :], fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="day")

                end

                if p_list[p] == "Δ"
                    samp = collect(5:0.2:40)
                    data_hyper = DataFrame(CSV.File("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/" * result_list_hyper[p][strain]))
                    c = data_hyper.c
                    cons_mean = zeros(2, length(samp))
                    cons_025 = zeros(2, length(samp))
                    cons_975 = zeros(2, length(samp))
                    for t in 1:length(samp)
                        cons = constant(c, [samp[t]])
                        cons_mean[:, t] = [samp[t], mean(cons)]
                        cons_025[:, t] = [samp[t], quantile(cons, 0.025)]
                        cons_975[:, t] = [samp[t], quantile(cons, 0.975)]
                    end
                    plot!(pl, subplot=p, cons_mean[1, :], cons_mean[2, :], color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, cons_025[1, :], cons_025[2, :], fillrange=cons_975[2, :], fillalpha=0.2, lw=0, color=color_palette[strain], label=false)
                    plot!(pl, subplot=p, ylabel="day")

                end
                if 5 < p
                    plot!(pl, subplot=p, xlabel="Salinity (PSU)", xtickfontcolor=:black)

                end
            end



            println(result_list_system[i][34] * " is done !")
            df = DataFrame(Salinity=[], Strain=[], r=[], μ=[], K=[], Kr=[], ϵ=[], ϕ=[], η=[], β=[], ρ=[], Δ=[])
            means = DataFrame(Salinity=[], Strain=[], r=[], μ=[], K=[], Kr=[], ϵ=[], ϕ=[], η=[], β=[], ρ=[], Δ=[])
            bars = DataFrame(Salinity=[], Strain=[], r=[], μ=[], K=[], Kr=[], ϵ=[], ϕ=[], η=[], β=[], ρ=[], Δ=[])

            sal = 1 #initialize salinity
            strain += 1 #move on to the next strain


            if strain > length(name_strain)
                strain = 1
            end

        end
    end
end


plot!(pl, size=(1300, 600), dpi=600, leftmargins=9mm)
plot!(subplot=1, top_margin=5mm)
plot!(subplot=10, bottom_margin=5mm, legend=true)


## PLOT HYPERFUNCTIONS
##IMPORT DATA


savefig(pl, "./MAIN_EXPERIMENT/HYPERFUNCTIONS/figures/" * "param_according_to_salinity_" * string(today()) * ".png")
pl
