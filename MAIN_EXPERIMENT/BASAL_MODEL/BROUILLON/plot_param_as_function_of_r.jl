########################################################
########## Laure L. / Main Exp. SALINITY  ##############
#################### Modelling #########################
########## make all plot and save in a file ############
########## David/Laure 03 July 2023 - Paris ############
########################################################

using StatsPlots, CSV, DataFrames, Plots, Measures, StatsBase, Distributions, Dates, DifferentialEquations

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

color_palette = ["#72DB23", "#12B4E0", "#DF321B", "#EB8203", "#000000"]
strain = ["A", "B", "C", "D"]
param = ["r", "K", "Kr", "μ", "ϵ", "ϕ", "η", "β"]
include("./moyenne_param_control_22May2023.jl")
include("./get_stat_param_system_04September2023.jl")
data = get_stat_system(result_list_control, result_list_system)
pl = plot(layout=(2, 4), margins=5mm, size=(1200, 600), legend=false, dpi=600)

for st in 1:length(strain)
    data_mean = filter(row -> row.condition[1:2] == strain[st] * "_", data[1])
    sort!(data_mean, :r)
    for p in 1:length(param)
        plot!(pl, subplot=p, title=param[p])
        plot!(pl, subplot=p, data_mean.r, data_mean[:, p+2], lw=2, color=color_palette[st], label=false)
        scatter!(pl, subplot=p, data_mean.r, data_mean[:, p+2], markersize=4, markerstrokewidth=0, color=color_palette[st], label=strain[st])
        if p > 4
            plot!(pl, subplot=p, xlabel="r (cell.ml⁻¹.day⁻¹)")
        end
    end
    println(data_mean.condition[1][1:1] * " IS DONE !")
end
plot!(pl, subplot=1, legend=true)
pl