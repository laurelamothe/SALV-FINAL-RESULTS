########################################################
########## Laure L. / Main Exp. SALINITY  ##############
#################### Modelling #########################
########## make all plot and save in a file ############
########## David/Laure 03 July 2023 - Paris ############
########################################################

#list of plot : 
#   chains
#   regression

## PACKAGES
using StatsPlots, CSV, DataFrames, Plots, Measures, StatsBase, Distributions, Dates, DifferentialEquations, Printf

color_palette = ["#77DF49", "#75C7E1", "#EB7753", "#FABC22", "#000000", "#CA8DDB"]
conditions = ["light", "obscurity"]
sal = collect(5:0.1:40)
lin_mean = zeros(2, length(sal))
lin_025 = zeros(2, length(sal))
lin_975 = zeros(2, length(sal))

## CHARGE DATA #
vd_results_lin = [
    DataFrame(CSV.File("./VIRUS_DECAY/results/SALEXP_virus_decay_light_2023-11-05.csv")),
    DataFrame(CSV.File("./VIRUS_DECAY/results/SALEXP_virus_decay_dark_2023-11-05.csv"))]

vd_data = DataFrame(CSV.File("./VIRUS_DECAY/data/data_input_virus_decay.csv"))
preg = scatter(layout=(1, 2), size=(800, 600), dpi=600, xlabel="Salinity (PSU)", legend=:none, ylims=(0.01, 0.05), ytickfontcolor=:white)
scatter!(preg, subplot=1, ylabel="Virus Decay \n(virus.mL⁻¹.day⁻¹)", ytickfontcolor=:black)


for i in 1:length(vd_results_lin)

    ## PLOT CHAINS ##
    pchain_lin = plot(layout=(1, 2), size=(600, 300), dpi=600, legend=false)
    plot!(pchain_lin, subplot=1, vd_results_lin[i].a, title="a", color=color_palette[2])
    plot!(pchain_lin, subplot=2, vd_results_lin[i].b, title="b", color=color_palette[2])
    savefig(pchain_lin, "./VIRUS_DECAY/figures/" * "chains_mcmc_virus_decay_" * conditions[i] * "_" * string(today()) * ".png")

    ## PLOT DATA AND REGRESSIONS ##
    scatter!(preg, subplot=i, vd_data.Salinity, vd_data[:, i+1], label="Experimental points", title=conditions[i], color=color_palette[4], markerstrokewidth=0, markersize=6)


    for t in 1:length(sal)
        lin = vd_results_lin[i].a .* sal[t] .+ vd_results_lin[i].b
        lin_mean[:, t] = [sal[t], mean(lin)]
        lin_025[:, t] = [sal[t], quantile(lin, 0.025)]
        lin_975[:, t] = [sal[t], quantile(lin, 0.975)]
    end

    plot!(preg, subplot=i, lin_mean[1, :], lin_mean[2, :], label="Regression : " * @sprintf("%.1E", mean(vd_results_lin[i].a)) * " * Salinity + " * @sprintf("%.1E", mean(vd_results_lin[i].b)), color=color_palette[4], lw=2)
    plot!(preg, subplot=i, lin_025[1, :], lin_025[2, :], fillrange=lin_975[2, :], fillalpha=0.2, lw=0, color=color_palette[4], label="CI 95%")
end
fake_plot = plot(size=(100, 400), grid=false, axis=false, rand(10), color=color_palette[4], lw=2, ylims=(10, 20), label="Regression")
scatter!(fake_plot, rand(10), color=color_palette[4], markerstrokewidth=0, markersize=6, label="Experimental points")
plot!(fake_plot, rand(10), fillrange=rand(10), fillalpha=0.2, lw=0, color=color_palette[4], label="CI 95%", legend=:top)

pl = plot(preg, fake_plot, layout=(2, 1), margins=0mm)
plot!(pl, subplot=1, leftmargin=5mm)
savefig(pl, "./VIRUS_DECAY/figures/" * "regressions_" * string(today()) * ".png")

pl