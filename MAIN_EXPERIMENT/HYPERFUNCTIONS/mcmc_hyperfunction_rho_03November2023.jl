######################################################
########## Laure L. / Main Exp. SALINITY  ############
#################### Modelling #######################
##### System - mcmc to estimate all parameters #######
####### David/Laure 24 August 2023 - Banyuls #########
###################################################### 
"""
## Package and processors
using Distributed
addprocs(4)
@everywhere using Plots, Measures, StatsPlots, Dates, Turing, Distributions, DifferentialEquations, Statistics, CSV, DataFrames, Glob
"""
## charge data and variables

@everywhere begin
    param = "ρ"
    color_palette = ["#F79135", "#70AD47", "#4DA2E1", "#A06DDD"]
    letter_strain = ["A", "B", "C", "D"]
    name_strain = ["RCC4221", "RCC1116", "RCC1123", "RCC1108"]

    list_parameters = ["r", "K", "Kr", "μ", "ϵ", "ϕ", "η", "β"]

    result_list_resi = [
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
        "SALEXP_R_onset_D_S6_2023-10-15.csv"
    ]

    include("get_stat_param_system_04September2023.jl")
    data = get_stat_resi(result_list_resi)
    data1 = data[1]#filter(row -> row.condition != "C_S2", data[1])

end



## define hyperfunctions
@everywhere begin

    @everywhere begin
        function constant(c, Sal)
            return c .* (Sal .* 0.0 .+ 1)
        end
    end

    @everywhere begin
        function linear(a, b, Sal)
            return a .* Sal .+ b
        end
    end

    @everywhere begin
        function poly(a, b, c, Sal)
            return a .* Sal .^ 2 .+ b .* Sal .+ c
        end
    end


    @everywhere begin
        function logistic(a, k, s0, Sal)
            return k ./ (1 .+ exp.(.-a .* (Sal .- s0)))
        end
    end


    @everywhere begin
        function normal(a, σ, µ, Sal)
            return a .* (1 ./ (σ .* sqrt(2 * pi)) .* exp.(-1 / 2 .* ((Sal .- µ) ./ σ) .^ 2))
        end
    end


    @everywhere begin
        function power(a, b, c, Sal)
            return a .* Sal .^ b .+ c
        end
    end


    @model function fit_hyperfunction_constant(data, salinity, sigma_mean, sigma_std)
        ## Priors
        c ~ truncated(LogNormal(log(6), 0.3), 0.0001, 15)

        σ ~ truncated(Normal(sigma_mean, sigma_std), 0, 10)

        ## Integration
        predicted_reg = constant(c, salinity)

        ## LIKELIHOOD
        for i = 1:length(data)          # Pour chaque point de données regarder ...
            data[i] ~ Normal(predicted_reg[i], σ)    # Likelihood --> compare log-likelihood avec les log(donnée).
        end
    end

    @model function fit_hyperfunction_linear(data, salinity, sigma_mean, sigma_std)
        ## Priors
        a ~ truncated(Cauchy(0, 1), -10, 10)
        b ~ truncated(LogNormal(log(6), 0.3), 0.0001, 15)

        σ ~ truncated(Normal(sigma_mean, sigma_std), 0, 10000)

        ## Integration
        predicted_reg = linear(a, b, salinity)

        ## LIKELIHOOD
        for i = 1:length(data)          # Pour chaque point de données regarder ...
            data[i] ~ Normal(predicted_reg[i], σ)    # Likelihood --> compare log-likelihood avec les log(donnée).
        end
    end


    @model function fit_hyperfunction_polynomial(data, salinity, sigma_mean, sigma_std)
        ## Priors
        a ~ truncated(Cauchy(0, 1), -10, 10)
        b ~ truncated(Cauchy(0, 1), -10, 10)
        c ~ truncated(LogNormal(log(6), 0.3), 0.0001, 15)

        σ ~ truncated(Normal(sigma_mean, sigma_std), 0, 10)

        ## Integration
        predicted_reg = poly(a, b, c, salinity)

        ## LIKELIHOOD
        for i = 1:length(data)          # Pour chaque point de données regarder ...
            data[i] ~ Normal(predicted_reg[i], σ)    # Likelihood --> compare log-likelihood avec les log(donnée).
        end
    end


    @model function fit_hyperfunction_logistic(data, salinity, sigma_mean, sigma_std)
        ## Priors
        a ~ truncated(LogNormal(log(1), 0.8), 1e-13, 5)
        k ~ truncated(LogNormal(log(4), 0.8), 0.01, 30)
        s0 ~ truncated(Cauchy(0, 0.3), -20, 40)

        σ ~ truncated(Normal(sigma_mean, sigma_std), 0, 10)

        ## Integration
        predicted_reg = logistic(a, k, s0, salinity)

        ## LIKELIHOOD
        for i = 1:length(data)          # Pour chaque point de données regarder ...
            data[i] ~ Normal(predicted_reg[i], σ)    # Likelihood --> compare log-likelihood avec les log(donnée).
        end
    end

    @model function fit_hyperfunction_normal(data, salinity, sigma_mean, sigma_std)
        ## Priors
        a ~ truncated(LogNormal(log(1000), 0.8), 10, 1000000)
        σm ~ truncated(LogNormal(log(10), 0.3), 0.001, 100)
        µ ~ truncated(Cauchy(10, 0.3), 1, 40)

        σ ~ truncated(Normal(sigma_mean, sigma_std), 0, 10)

        ## Integration
        predicted_reg = normal(a, σm, µ, salinity)

        ## LIKELIHOOD
        for i = 1:length(data)          # Pour chaque point de données regarder ...
            data[i] ~ Normal(predicted_reg[i], σ)    # Likelihood --> compare log-likelihood avec les log(donnée).
        end
    end

    @model function fit_hyperfunction_power(data, salinity, sigma_mean, sigma_std)
        ## Priors
        a ~ truncated(LogNormal(log(10), 0.8), 0.01, 100)
        b ~ truncated(LogNormal(log(0.8), 0.8), 1e-13, 1)
        c ~ truncated(LogNormal(log(1), 0.8), 1e-3, 30)

        σ ~ truncated(Normal(sigma_mean, sigma_std), 0, 10)

        ## Integration
        predicted_reg = power(a, b, c, salinity)

        ## LIKELIHOOD
        for i = 1:length(data)          # Pour chaque point de données regarder ...
            data[i] ~ Normal(predicted_reg[i], σ)    # Likelihood --> compare log-likelihood avec les log(donnée).
        end
    end
    list_model = [fit_hyperfunction_constant, fit_hyperfunction_linear, fit_hyperfunction_polynomial, fit_hyperfunction_logistic, fit_hyperfunction_normal, fit_hyperfunction_power]
    list_model_names = ["constant", "linear", "polynomial", "logistic", "normal", "power"]
    list_model_hyperparam = [["c", "σ"], ["a", "b", "σ"], ["a", "b", "c", "σ"], ["a", "k", "s0", "σ"], ["a", "σm", "µ", "σ"], ["a", "b", "c", "σ"]]
end

## run MCMC
@everywhere begin
    nburn = 2000
    nstep = 2000
    nchain = 4
end

function runMCMC(fitting_process, name_model, strain, param, data, nburn, nstep, nchain, list_hyper_param)
    data_mean = filter(row -> row.condition[1:2] == strain * "_", data1)[!, param]
    salinities = filter(row -> row.condition[1:2] == strain * "_", data1).salinity
    sigma_mean = data[3][param]
    sigma_std = data[4][param]

    modelfit = fitting_process(data_mean, salinities, sigma_mean, sigma_std)
    chainfit = reduce(chainscat, pmap(x -> sample(modelfit, NUTS(nburn, 0.65), nstep, save_state=false), 1:nchain))
    chain_df = DataFrame(Array(chainfit), :auto)
    rename!(chain_df, list_hyper_param)
    stat_chain = ess_rhat(chainfit)

    CSV.write("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/SALEXP_chains_hyperfunction_" * name_model * "_" * param * "_" * strain * "_" * string(today()) * ".csv", chain_df, writeheader=true)
    CSV.write("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/SALEXP_statchains_hyperfunction_" * name_model * "_" * param * "_" * strain * "_" * string(today()) * ".csv", DataFrame(stat_chain), writeheader=true)

    ll = [loglikelihood(modelfit, j) for j in eachrow(DataFrame(chainfit))]
    k = size(chain_df)[2] - 1
    n = length(data_mean)

    AIC = 2 * k - 2 * mean(ll)
    AICc = 2 * k - 2 * mean(ll) + (2 * k) * (k + 1) / (n - k - 1)
    BIC = k * log(n) - 2 * mean(ll)

    return (chain_df, stat_chain, AIC, AICc, BIC, data_mean, salinities)
end

begin
    AIC_list = DataFrame(strain=[], model=[], AIC=[])
    AICc_list = DataFrame(strain=[], model=[], AICc=[])
    BIC_list = DataFrame(strain=[], model=[], BIC=[])


    p_reg = plot(layout=(2, 3), size=(1000, 600), margins=5mm, dpi=600, ylabel="ρ", xlabel="Salinity (PSU)")
    p_chain_p = plot(layout=(1, 3), size=(900, 300), dpi=600)
    p_chain_lo = plot(layout=(1, 3), size=(900, 300), dpi=600)
    p_chain_l = plot(layout=(1, 2), size=(600, 300), dpi=600)
    p_chain_n = plot(layout=(1, 4), size=(1200, 300), dpi=600)
    p_chain_c = plot(layout=(1, 1), size=(300, 300), dpi=600)

    for j in 1:length(list_model)
        for i in 1:length(letter_strain)
            calib = runMCMC(list_model[j], list_model_names[j], letter_strain[i], param, data, nburn, nstep, nchain, list_model_hyperparam[j])
            AIC_list = vcat(AIC_list, DataFrame(strain=[name_strain[i]], model=list_model_names[j], AIC=[calib[3]]))
            AICc_list = vcat(AICc_list, DataFrame(strain=[name_strain[i]], model=list_model_names[j], AICc=[calib[4]]))
            BIC_list = vcat(BIC_list, DataFrame(strain=[name_strain[i]], model=list_model_names[j], BIC=[calib[5]]))

            if list_model_names[j] == "constant"
                plot!(p_chain_c, calib[1].c[:], title="c", label=name_strain[i])

                scatter!(p_reg, subplot=1, calib[7], calib[6], markerstrokewidth=0, color=color_palette[i], label=name_strain[i], title="Constant regression")
                samp = collect(5:0.2:40)
                c = calib[1].c
                cons_mean = zeros(2, length(samp))
                cons_025 = zeros(2, length(samp))
                cons_975 = zeros(2, length(samp))
                for t in 1:length(samp)
                    cons = constant(c, [samp[t]])
                    cons_mean[:, t] = [samp[t], mean(cons)]
                    cons_025[:, t] = [samp[t], quantile(cons, 0.025)]
                    cons_975[:, t] = [samp[t], quantile(cons, 0.975)]
                end
                plot!(p_reg, subplot=1, cons_mean[1, :], cons_mean[2, :], color=color_palette[i], label=false)
                plot!(p_reg, subplot=1, cons_025[1, :], cons_025[2, :], fillrange=cons_975[2, :], fillalpha=0.2, lw=0, color=color_palette[i], label=false)
            elseif list_model_names[j] == "linear"
                plot!(p_chain_l, subplot=1, calib[1].a[:], title="a", label=name_strain[i])
                plot!(p_chain_l, subplot=2, calib[1].b[:], title="b", label=name_strain[i])

                scatter!(p_reg, subplot=2, calib[7], calib[6], markerstrokewidth=0, color=color_palette[i], label=name_strain[i], title="Linear regression")
                samp = collect(5:0.2:40)
                a = calib[1].a
                b = calib[1].b
                lin_mean = zeros(2, length(samp))
                lin_025 = zeros(2, length(samp))
                lin_975 = zeros(2, length(samp))
                for t in 1:length(samp)
                    lin = linear(a, b, [samp[t]])
                    lin_mean[:, t] = [samp[t], mean(lin)]
                    lin_025[:, t] = [samp[t], quantile(lin, 0.025)]
                    lin_975[:, t] = [samp[t], quantile(lin, 0.975)]
                end
                plot!(p_reg, subplot=2, lin_mean[1, :], lin_mean[2, :], color=color_palette[i], label=false)
                plot!(p_reg, subplot=2, lin_025[1, :], lin_025[2, :], fillrange=lin_975[2, :], fillalpha=0.2, lw=0, color=color_palette[i], label=false)

            elseif list_model_names[j] == "polynomial"
                plot!(p_chain_p, subplot=1, calib[1].a[:], title="a", label=name_strain[i])
                plot!(p_chain_p, subplot=2, calib[1].b[:], title="b", label=name_strain[i])
                plot!(p_chain_p, subplot=3, calib[1].c[:], title="c", label=name_strain[i])

                scatter!(p_reg, subplot=3, calib[7], calib[6], markerstrokewidth=0, color=color_palette[i], label=name_strain[i], title="Polynomial regression")
                samp = collect(5:0.2:40)
                a = calib[1].a
                b = calib[1].b
                c = calib[1].c
                p_mean = zeros(2, length(samp))
                p_025 = zeros(2, length(samp))
                p_975 = zeros(2, length(samp))
                for t in 1:length(samp)
                    pow = poly(a, b, c, [samp[t]])
                    p_mean[:, t] = [samp[t], mean(pow)]
                    p_025[:, t] = [samp[t], quantile(pow, 0.025)]
                    p_975[:, t] = [samp[t], quantile(pow, 0.975)]
                end
                plot!(p_reg, subplot=3, p_mean[1, :], p_mean[2, :], color=color_palette[i], label=false)
                plot!(p_reg, subplot=3, p_025[1, :], p_025[2, :], fillrange=p_975[2, :], fillalpha=0.2, lw=0, color=color_palette[i], label=false)

            elseif list_model_names[j] == "power"
                plot!(p_chain_p, subplot=1, calib[1].a[:], title="a", label=name_strain[i])
                plot!(p_chain_p, subplot=2, calib[1].b[:], title="b", label=name_strain[i])

                scatter!(p_reg, subplot=6, calib[7], calib[6], markerstrokewidth=0, color=color_palette[i], label=name_strain[i], title="Power regression")
                samp = collect(5:0.2:40)
                a = calib[1].a
                b = calib[1].b
                c = calib[1].c
                p_mean = zeros(2, length(samp))
                p_025 = zeros(2, length(samp))
                p_975 = zeros(2, length(samp))
                for t in 1:length(samp)
                    pow = power(a, b, c, [samp[t]])
                    p_mean[:, t] = [samp[t], mean(pow)]
                    p_025[:, t] = [samp[t], quantile(pow, 0.025)]
                    p_975[:, t] = [samp[t], quantile(pow, 0.975)]
                end
                plot!(p_reg, subplot=6, p_mean[1, :], p_mean[2, :], legend=false, color=color_palette[i], label=false)
                plot!(p_reg, subplot=6, p_025[1, :], p_025[2, :], fillrange=p_975[2, :], fillalpha=0.2, lw=0, color=color_palette[i], label=false)

            elseif list_model_names[j] == "logistic"
                plot!(p_chain_lo, subplot=1, calib[1].a[:], title="a", label=name_strain[i])
                plot!(p_chain_lo, subplot=2, calib[1].k[:], title="k", label=name_strain[i])
                plot!(p_chain_lo, subplot=3, calib[1].k[:], title="s0", label=name_strain[i])

                scatter!(p_reg, subplot=4, calib[7], calib[6], markerstrokewidth=0, color=color_palette[i], label=name_strain[i], title="Logistic regression")
                samp = collect(5:0.2:40)
                a = calib[1].a
                k = calib[1].k
                s0 = calib[1].s0
                logi_mean = zeros(2, length(samp))
                logi_025 = zeros(2, length(samp))
                logi_975 = zeros(2, length(samp))
                for t in 1:length(samp)
                    logi = logistic(a, k, s0, [samp[t]])
                    logi_mean[:, t] = [samp[t], mean(logi)]
                    logi_025[:, t] = [samp[t], quantile(logi, 0.025)]
                    logi_975[:, t] = [samp[t], quantile(logi, 0.975)]
                end
                plot!(p_reg, subplot=4, logi_mean[1, :], logi_mean[2, :], color=color_palette[i], label=false)
                plot!(p_reg, subplot=4, logi_025[1, :], logi_025[2, :], fillrange=logi_975[2, :], fillalpha=0.2, lw=0, color=color_palette[i], label=false)
            else
                plot!(p_chain_n, subplot=1, calib[1].a[:], title="a", label=name_strain[i])
                # plot!(p_chain_n, subplot=2, calib[1].b[:], title="b", label=name_strain[i])
                plot!(p_chain_n, subplot=3, calib[1].σm[:], title="σ", label=name_strain[i])
                plot!(p_chain_n, subplot=4, calib[1][:, 3], title="µ", label=name_strain[i])

                scatter!(p_reg, subplot=5, calib[7], calib[6], markerstrokewidth=0, color=color_palette[i], label=name_strain[i], title="Normal regression")
                samp = collect(5:0.2:40)
                a = calib[1].a
                # b = calib[1].b
                σm = calib[1].σm
                µ = calib[1][:, 3]
                norm_mean = zeros(2, length(samp))
                norm_025 = zeros(2, length(samp))
                norm_975 = zeros(2, length(samp))
                for t in 1:length(samp)
                    norm = normal(a, σm, µ, [samp[t]])
                    norm_mean[:, t] = [samp[t], mean(norm)]
                    norm_025[:, t] = [samp[t], quantile(norm, 0.025)]
                    norm_975[:, t] = [samp[t], quantile(norm, 0.975)]
                end
                plot!(p_reg, subplot=5, norm_mean[1, :], norm_mean[2, :], color=color_palette[i], label=false)
                plot!(p_reg, subplot=5, norm_025[1, :], norm_025[2, :], fillrange=norm_975[2, :], fillalpha=0.2, lw=0, color=color_palette[i], label=false)

            end
        end

    end
    plot_crit = plot(layout=(1, 3), size=(900, 300), dpi=600)
    scatter!(plot_crit, subplot=1,
        filter(row -> row.model == "constant", AIC_list).strain,
        filter(row -> row.model == "constant", AIC_list).AIC,
        color=color_palette[1],
        label="constant", title="AIC", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=1,
        filter(row -> row.model == "linear", AIC_list).strain,
        filter(row -> row.model == "linear", AIC_list).AIC,
        color=color_palette[2],
        label="linear", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=1,
        filter(row -> row.model == "polynomial", AIC_list).strain,
        filter(row -> row.model == "polynomial", AIC_list).AIC,
        color=color_palette[3],
        label="polynomial", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=1,
        filter(row -> row.model == "logistic", AIC_list).strain,
        filter(row -> row.model == "logistic", AIC_list).AIC,
        color=color_palette[4],
        label="logistic", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=1,
        filter(row -> row.model == "normal", AIC_list).strain,
        filter(row -> row.model == "normal", AIC_list).AIC,
        color=:black,
        label="normal", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=1,
        filter(row -> row.model == "power", AIC_list).strain,
        filter(row -> row.model == "power", AIC_list).AIC,
        color=:pink,
        label="power", markersize=6, markerstrokewidth=0)

    scatter!(plot_crit, subplot=2,
        filter(row -> row.model == "constant", AICc_list).strain,
        filter(row -> row.model == "constant", AICc_list).AICc,
        color=color_palette[1],
        label="constant", title="AICc", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=2,
        filter(row -> row.model == "linear", AICc_list).strain,
        filter(row -> row.model == "linear", AICc_list).AICc,
        color=color_palette[2],
        label="linear", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=2,
        filter(row -> row.model == "polynomial", AICc_list).strain,
        filter(row -> row.model == "polynomial", AICc_list).AICc,
        color=color_palette[3],
        label="polynomial", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=2,
        filter(row -> row.model == "logistic", AICc_list).strain,
        filter(row -> row.model == "logistic", AICc_list).AICc,
        color=color_palette[4],
        label="logistic", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=2,
        filter(row -> row.model == "normal", AICc_list).strain,
        filter(row -> row.model == "normal", AICc_list).AICc,
        color=:black,
        label="normal", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=2,
        filter(row -> row.model == "power", AICc_list).strain,
        filter(row -> row.model == "power", AICc_list).AICc,
        color=:pink,
        label="power", markersize=6, markerstrokewidth=0)

    scatter!(plot_crit, subplot=3,
        filter(row -> row.model == "constant", BIC_list).strain,
        filter(row -> row.model == "constant", BIC_list).BIC,
        color=color_palette[1],
        label="constant", title="BIC", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=3,
        filter(row -> row.model == "linear", BIC_list).strain,
        filter(row -> row.model == "linear", BIC_list).BIC,
        color=color_palette[2],
        label="linear", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=3,
        filter(row -> row.model == "polynomial", BIC_list).strain,
        filter(row -> row.model == "polynomial", BIC_list).BIC,
        color=color_palette[3],
        label="polynomial", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=3,
        filter(row -> row.model == "logistic", BIC_list).strain,
        filter(row -> row.model == "logistic", BIC_list).BIC,
        color=color_palette[4],
        label="logistic", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=3,
        filter(row -> row.model == "normal", BIC_list).strain,
        filter(row -> row.model == "normal", BIC_list).BIC,
        color=:black,
        label="normal", markersize=6, markerstrokewidth=0)
    scatter!(plot_crit, subplot=3,
        filter(row -> row.model == "power", BIC_list).strain,
        filter(row -> row.model == "power", BIC_list).BIC,
        color=:pink,
        label="power", markersize=6, markerstrokewidth=0)

    mean_crit = DataFrame(
        model=list_model_names,
        AIC=[mean(filter(row -> row.model == "constant", AIC_list).AIC),
            mean(filter(row -> row.model == "linear", AIC_list).AIC),
            mean(filter(row -> row.model == "polynomial", AIC_list).AIC),
            mean(filter(row -> row.model == "logistic", AIC_list).AIC),
            mean(filter(row -> row.model == "normal", AIC_list).AIC),
            mean(filter(row -> row.model == "power", AIC_list).AIC)],
        AICc=[mean(filter(row -> row.model == "constant", AICc_list).AICc),
            mean(filter(row -> row.model == "linear", AICc_list).AICc),
            mean(filter(row -> row.model == "polynomial", AICc_list).AICc),
            mean(filter(row -> row.model == "logistic", AICc_list).AICc),
            mean(filter(row -> row.model == "normal", AICc_list).AICc),
            mean(filter(row -> row.model == "power", AICc_list).AICc)], BIC=[mean(filter(row -> row.model == "constant", BIC_list).BIC),
            mean(filter(row -> row.model == "linear", BIC_list).BIC),
            mean(filter(row -> row.model == "polynomial", BIC_list).BIC),
            mean(filter(row -> row.model == "logistic", BIC_list).BIC),
            mean(filter(row -> row.model == "normal", BIC_list).BIC),
            mean(filter(row -> row.model == "power", BIC_list).BIC)])


    CSV.write("./MAIN_EXPERIMENT/HYPERFUNCTIONS/results/SALEXP_criterion_models_" * param * "_" * string(today()) * ".csv", mean_crit, writeheader=true)

    savefig(p_reg, "./MAIN_EXPERIMENT/HYPERFUNCTIONS/figures/" * "hyperfunction_regression_" * param * "_" * string(today()) * ".png")
    #savefig(p_chain_c, "./MAIN_EXPERIMENT/HYPERFUNCTIONS/figures/" * "hyperfunction_constant_chains_" * param * "_" * string(today()) * ".png")
    #savefig(p_chain_l, "./MAIN_EXPERIMENT/HYPERFUNCTIONS/figures/" * "hyperfunction_linear_chains_" * param * "_" * string(today()) * ".png")
    #savefig(p_chain_p, "./MAIN_EXPERIMENT/HYPERFUNCTIONS/figures/" * "hyperfunction_polynomial_chains_" * param * "_" * string(today()) * ".png")
    savefig(plot_crit, "./MAIN_EXPERIMENT/HYPERFUNCTIONS/figures/" * "hyperfunction_criteria_" * param * "_" * string(today()) * ".png")


end

