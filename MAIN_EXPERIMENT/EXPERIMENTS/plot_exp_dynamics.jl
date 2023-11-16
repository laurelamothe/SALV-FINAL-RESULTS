########################################################
########## Laure L. / Main Exp. SALINITY  ##############
#################### Modelling #########################
########## make all plot and save in a file ############
########## David/Laure 03 July 2023 - Paris ############
########################################################

#list of plot : 
#   chains
#   rhats
#   corr
#   dynamics
#   prior posterior
using StatsPlots, CSV, DataFrames, Plots, Measures, StatsBase, Distributions, Dates, DifferentialEquations

color_palette=["#77DF49", "#75C7E1", "#EB7753", "#FABC22", "#000000", "#CA8DDB"]
grad_host_control = palette(["#A7F381", "#59E115", "#2A6A0A"], 6)
grad_host_infected = palette(["#FCA68E", "#F73A07", "#821F04"], 6)
grad_virus = palette(["#FCDD8E", "#D99C05", "#6A4C02"], 6)

salinity = [5, 10, 15, 25, 35, 40]
strain = ["RCC4221", "RCC1116", "RCC1123", "RCC1108"]
#IMPORT DATA
#define column types
type_list = [String127, Float64, String127, Int128, Int128, String127, String127, String127, String127, Float64, Float64, Float64, Float64]
#read CSV
data = DataFrame(CSV.File("./MAIN_EXPERIMENT/EXPERIMENTS/raw_data_main_experiment.csv", types=type_list))


pl = plot(layout=grid(5, 3, heights=[0.24, 0.24, 0.24, 0.24, 0.04]), size=(900, 900), dpi=600, xlims=(0, 14), ylims=(1e3, 1e10), yscale=:log10, tickfontcolor=:white, legend=false, grid=true)
plot!(grid=:y)
for i in 1:4 #each STRAIN
    plot!(subplot=1 + 3 * (i - 1), ylabel="particle.mL⁻¹", ytickfontcolor=:black)
    for j in 1:6 #each salinity
        ##HOST CONTROLE
        data_temp1 = filter(row -> row.STRAIN == strain[i] && row.SALINITY == salinity[j] && row.POPULATION == "Host" && row.CONTROL == "Yes" && row.REPLICATE == 1, data)
        data_temp2 = filter(row -> row.STRAIN == strain[i] && row.SALINITY == salinity[j] && row.POPULATION == "Host" && row.CONTROL == "Yes" && row.REPLICATE == 2, data)
        means = (data_temp1.CONCENTRATION .+ data_temp1.CONCENTRATION) ./ 2
        bars = abs.(data_temp1.CONCENTRATION .- data_temp2.CONCENTRATION)
        plot!(subplot=1 + 3 * (i - 1), data_temp1.TIME, means, color=grad_host_control[j], lw=1.5, label=:none)
        scatter!(subplot=1 + 3 * (i - 1), data_temp1.TIME, means, color=grad_host_control[j], markerstrokecolor=grad_host_control[j], yerr=bars, markeralpha=0, label=:none, markerstrokewidth=1.5)
        scatter!(subplot=1 + 3 * (i - 1), data_temp1.TIME, means, color=grad_host_control[j], markerstrokewidth=0, label=string(salinity[j]) * "PSU", markersize=4)


        ##HOST INFECTED
        data_temp1 = filter(row -> row.STRAIN == strain[i] && row.SALINITY == salinity[j] && row.POPULATION == "Host" && row.CONTROL == "No" && row.REPLICATE == 1, data)
        data_temp2 = filter(row -> row.STRAIN == strain[i] && row.SALINITY == salinity[j] && row.POPULATION == "Host" && row.CONTROL == "No" && row.REPLICATE == 2, data)
        means = (data_temp1.CONCENTRATION .+ data_temp1.CONCENTRATION) ./ 2
        bars = abs.(data_temp1.CONCENTRATION .- data_temp2.CONCENTRATION)
        plot!(subplot=2 + 3 * (i - 1), data_temp1.TIME, means, color=grad_host_infected[j], lw=1.5, label=:none)
        scatter!(subplot=2 + 3 * (i - 1), data_temp1.TIME, means, color=grad_host_infected[j], markerstrokecolor=grad_host_infected[j], yerr=bars, markeralpha=0, markerstrokewidth=1.5, label=:none)
        scatter!(subplot=2 + 3 * (i - 1), data_temp1.TIME, means, color=grad_host_infected[j], markerstrokewidth=0, label=string(salinity[j]) * "PSU", markersize=4)


        ##VIRUS
        data_temp1 = filter(row -> row.STRAIN == strain[i] && row.SALINITY == salinity[j] && row.POPULATION == "Virus" && row.REPLICATE == 1, data)
        data_temp2 = filter(row -> row.STRAIN == strain[i] && row.SALINITY == salinity[j] && row.POPULATION == "Virus" && row.REPLICATE == 2, data)
        means = (data_temp1.CONCENTRATION .+ data_temp1.CONCENTRATION) ./ 2
        bars = abs.(data_temp1.CONCENTRATION .- data_temp2.CONCENTRATION)
        plot!(subplot=3 + 3 * (i - 1), data_temp1.TIME, means, color=grad_virus[j], lw=1.5, label=:none)
        scatter!(subplot=3 + 3 * (i - 1), data_temp1.TIME, means, color=grad_virus[j], markerstrokecolor=grad_virus[j], yerr=bars, markeralpha=0, markerstrokewidth=1.5, label=:none)
        scatter!(subplot=3 + 3 * (i - 1), data_temp1.TIME, means, color=grad_virus[j], markerstrokewidth=0, label=string(salinity[j]) * "PSU", markersize=4)

    end
end


plot!(subplot=1, title="Host (control cultures)")
plot!(subplot=2, title="Host (infected cultures)")
plot!(subplot=3, title="Virus (infected cultures)")

#plot!(subplot=10, legend=:outerbottom, legend_columns=4)
#plot!(subplot=11, legend=:outerbottom, legend_columns=4)
#plot!(subplot=12, legend=:outerbottom, legend_columns=4)


plot!(subplot=10, xlabel="day", xtickfontcolor=:black)
plot!(subplot=11, xlabel="day", xtickfontcolor=:black)
plot!(subplot=12, xlabel="day", xtickfontcolor=:black)


for j in 1:6
    scatter!(subplot=13, [1:2], [0, 0],
        markerstrokewidth=0,
        label=string(salinity[j]) * "PSU",
        color=grad_host_control[j], legend=:top,
        legend_columns=3,
        axis=false)

    scatter!(subplot=14, [1:2], [0, 0],
        markerstrokewidth=0,
        label=string(salinity[j]) * "PSU",
        color=grad_host_infected[j], legend=:top,
        legend_columns=3,
        axis=false)

    scatter!(subplot=15, [1:2], [0, 0],
        markerstrokewidth=0,
        label=string(salinity[j]) * "PSU",
        color=grad_virus[j], legend=:top,
        legend_columns=3,
        axis=false)
end
pl

savefig(pl, "./MAIN_EXPERIMENT/EXPERIMENTS/figures/experimental_dynamics_" * string(today()) * ".png")