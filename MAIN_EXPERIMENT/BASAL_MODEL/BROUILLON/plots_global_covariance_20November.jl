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
salinity = [5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40, 5, 10, 15, 25, 35, 40]
param = ["r", "K", "Kr", "μ", "ϵ", "ϕ", "η", "β"]
include("./moyenne_param_control_22May2023.jl")
result_virus_decay = DataFrame(CSV.File("./VIRUS_DECAY/results/SALEXP_virus_decay_light_2023-11-05.csv"))
pl = Plots.plot(layout=(size(data)[2], size(data)[2]), size=(1000, 850), margin=0mm, left_margin=5mm, xscale=:log, yscale=:log)
global_data = DataFrame()
correlation_table = DataFrame(r=zeros(8), K=zeros(8), Kr=zeros(8), μ=zeros(8), ϵ=zeros(8), ϕ=zeros(8), η=zeros(8), β=zeros(8))
cospe = DataFrame(r=zeros(8), K=zeros(8), Kr=zeros(8), μ=zeros(8), ϵ=zeros(8), ϕ=zeros(8), η=zeros(8), β=zeros(8))


for file in 1:length(result_list_system)

    #import data
    data_system = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/" * result_list_system[file]))
    data_system.ϕ = exp.(data_system.ϕ)
    data_system = hcat(DataFrame(
            Kr=data_system.Kr,
            µ=data_system.µ
        ), data_system[:, 3:11])
    data_control = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/control/" * result_list_control[file]))
    data = hcat(data_control[:, 1:2], data_system[:, 1:6])
    #global_data = vcat(global_data, data)


    for i in 1:size(data)[2]
        for j in 1:size(data)[2]
            cospe[i, j] = corspearman(data[:, j], data[:, i])
        end
    end
    correlation_table = correlation_table .+ abs.(cospe)

    if salinity[file] == 40
        CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/BROUILLON/results/spearman_corr_" * result_list_system[file][34:34] * "_" * string(today()) * ".csv", correlation_table, writeheader=true)
        correlation_table = DataFrame(r=zeros(8), K=zeros(8), Kr=zeros(8), μ=zeros(8), ϵ=zeros(8), ϕ=zeros(8), η=zeros(8), β=zeros(8))
    end
    println(result_list_system[file][34:37] * " IS DONE !")
end


"""
index = rand(1:size(global_data)[1], 3000)

for i in 1:size(global_data)[2]
    for j in 1:size(global_data)[2]
        scatter!(pl, subplot=j + size(global_data)[2] * (i - 1),
            global_data[index, j], global_data[index, i],
            markeralpha=0.05,
            markercolor=:black)
        plot!(pl, ticks=false, legend=:none)

        correlation_table[i, j+1] = corspearman(global_data[:, j], global_data[:, i])

        if i == 1
            plot!(pl, subplot=j, title=names(global_data)[j])
        end
        if j == 1
            plot!(pl, subplot=1 + size(global_data)[2] * (i - 1), ylabel=names(global_data)[i])
        end
    end
end
savefig(pl, "./MAIN_EXPERIMENT/BASAL_MODEL/BROUILLON/figures/plot_covariance" * string(today()) * ".png")
CSV.write("./MAIN_EXPERIMENT/BASAL_MODEL/BROUILLON/results/spearman_corr_" * string(today()) * ".csv", correlation_table, writeheader=true)
"""