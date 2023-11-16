######################################################
########## Laure L. / Main Exp. SALINITY  ############
################## Hyperfunction #####################
############## test the significance of ##############
##### parameter variation according to salinity ######
######### Laure 08 October 2023 - Roscoff ############
###################################################### 

using Plots, Measures, StatsPlots, Dates, Statistics, CSV, DataFrames, Distributions, HypothesisTests

color_palette = ["#F79135", "#70AD47", "#4DA2E1", "#A06DDD"]
letter_strain = ["A", "B", "C", "D"]
name_strain = ["RCC4221", "RCC1116", "RCC1123", "RCC1108"]
salinities = [5, 10, 15, 25, 35, 40]

list_parameters = ["r", "K", "Kr", "μ", "ϵ", "ϕ", "η", "β"]
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

#1) charge the data of one parameter for one strain and every salinity

strain = 'B'
param = "Kr"
df = DataFrame(S5=ones(8000), S10=ones(8000), S15=ones(8000), S25=ones(8000), S35=ones(8000), S40=ones(8000))
x = filter(row -> row[34] == strain, result_list_system)
y = filter(row -> row[23] == strain, result_list_control)
tmp = []
pl = plot();
#for strain in letter_strain
#   for param in list_parameters
for sal in 1:length(x)
    data = hcat(DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/control/" * y[sal]))[:, 1:2], DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/" * x[sal])))
    index = rand(1:8000, 100)
    df[:, sal] = data[:, param]
    append!(tmp, [data[index, param]])
    density!(pl, df[index, sal])
end
#CSV.write("./MAIN_EXPERIMENT/basal_model_result_" * param * "_" * strain * "_" * string(today()) * ".csv", df, writeheader=true)
#  end
#end

pl











"""
for sal in 1:length(tmp)
    println("LEVENES TEST FOR STRAIN " * strain * " PARAMETER " * param * " AND SALINITY " * string(salinities[sal]) * " : ")
    println(ExactOneSampleKSTest(tmp[sal], Normal(mean(tmp[sal]), sqrt(var(tmp[sal])))))
end
tmp = append!(tmp[1:3], tmp[5:6])
LeveneTest(tmp[3:6]...)


result = OneWayANOVATest(tmp...)
"""

