function get_stat_system(result_list_c, result_list_s)
    salinities = [5, 10, 15, 25, 35, 40]

    fitted_param_mean = DataFrame(
        condition=(""^length(result_list_c)),
        salinity=(zeros(length(result_list_c))),
        r=(zeros(length(result_list_c))),
        K=(zeros(length(result_list_c))),
        Kr=(zeros(length(result_list_s))),
        μ=(zeros(length(result_list_s))),
        ϵ=(zeros(length(result_list_s))),
        ϕ=(zeros(length(result_list_s))),
        η=(zeros(length(result_list_s))),
        β=(zeros(length(result_list_s))),
        Rper=(zeros(length(result_list_s))),
        N0_rep1=(zeros(length(result_list_s))),
        N0_rep2=(zeros(length(result_list_s))),
        V0_rep1=(zeros(length(result_list_s))),
        V0_rep2=(zeros(length(result_list_s))),
        #σ=(zeros(length(result_list)))
    )


    fitted_param_var = DataFrame(
        condition=(""^length(result_list_c)),
        salinity=(zeros(length(result_list_c))),
        r=(zeros(length(result_list_c))),
        K=(zeros(length(result_list_c))),
        Kr=(zeros(length(result_list_s))),
        μ=(zeros(length(result_list_s))),
        ϵ=(zeros(length(result_list_s))),
        ϕ=(zeros(length(result_list_s))),
        η=(zeros(length(result_list_s))),
        β=(zeros(length(result_list_s))),
        Rper=(zeros(length(result_list_s))),
        N0_rep1=(zeros(length(result_list_s))),
        N0_rep2=(zeros(length(result_list_s))),
        V0_rep1=(zeros(length(result_list_s))),
        V0_rep2=(zeros(length(result_list_s))),
        #σ=(zeros(length(result_list)))
    )

    for i in 1:length(result_list_c)
        df_c = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/control/" * result_list_c[i]))
        df_c.K = df_c.K ./ 1e7
        df_s = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/system/" * result_list_s[i]))
        df_s.Kr = df_s.Kr ./ 1e7
        fitted_param_mean.condition[i] = result_list_c[i][23:26]
        fitted_param_mean.salinity[i] = salinities[parse(Int128, result_list_c[i][26])]
        fitted_param_mean.r[i] = mean(df_c.r)
        fitted_param_mean.K[i] = mean(df_c.K)
        fitted_param_mean.Kr[i] = mean(df_s.Kr)
        fitted_param_mean.μ[i] = mean(df_s.μ)
        fitted_param_mean.ϵ[i] = mean(df_s.ϵ)
        fitted_param_mean.ϕ[i] = mean(df_s.ϕ)
        fitted_param_mean.η[i] = mean(df_s.η)
        fitted_param_mean.β[i] = mean(df_s.β)
        fitted_param_mean.Rper[i] = mean(df_s.Rper)
        fitted_param_mean.N0_rep1[i] = mean(df_s.N0_rep1)
        fitted_param_mean.N0_rep2[i] = mean(df_s.N0_rep2)
        fitted_param_mean.V0_rep1[i] = mean(df_s.V0_rep1)
        fitted_param_mean.V0_rep2[i] = mean(df_s.V0_rep2)
        #fitted_param_mean.σ[i] = mean(df_s.σ)

        fitted_param_var.condition[i] = result_list_c[i][23:26]
        fitted_param_mean.salinity[i] = salinities[parse(Int128, result_list_c[i][26])]
        fitted_param_var.r[i] = var(df_c.r)
        fitted_param_var.K[i] = var(df_c.K)
        fitted_param_var.Kr[i] = var(df_s.Kr)
        fitted_param_var.μ[i] = var(df_s.μ)
        fitted_param_var.ϵ[i] = var(df_s.ϵ)
        fitted_param_var.ϕ[i] = var(df_s.ϕ)
        fitted_param_var.η[i] = var(df_s.η)
        fitted_param_var.β[i] = var(df_s.β)
        fitted_param_var.Rper[i] = var(df_s.Rper)
        fitted_param_var.N0_rep1[i] = var(df_s.N0_rep1)
        fitted_param_var.N0_rep2[i] = var(df_s.N0_rep2)
        fitted_param_var.V0_rep1[i] = var(df_s.V0_rep1)
        fitted_param_var.V0_rep2[i] = var(df_s.V0_rep2)
        #fitted_param_var.σ[i] = var(df.σ)
    end
    mean_of_var = Dict()
    std_of_var = Dict()
    for i in ["A", "B", "C", "D"]
        temp = []
        temp = eachcol(filter(row -> row.condition[1:2] == i * "_", fitted_param_var))
        for j in 3:length(temp)
            mean_of_var[names(fitted_param_var)[j]] = mean(temp[j])
            std_of_var[names(fitted_param_var)[j]] = std(temp[j])
        end
    end

    return (fitted_param_mean, fitted_param_var, mean_of_var, std_of_var)
end




























function get_stat_resi(result_list_r)
    salinities = [5, 10, 15, 25, 35, 40]

    fitted_param_mean = DataFrame(
        condition=(""^length(result_list_r)),
        salinity=(zeros(length(result_list_r))),
        ρ=(zeros(length(result_list_r))),
        Δ=(zeros(length(result_list_r)))
    )


    fitted_param_var = DataFrame(
        condition=(""^length(result_list_r)),
        salinity=(zeros(length(result_list_r))),
        ρ=(zeros(length(result_list_r))),
        Δ=(zeros(length(result_list_r)))
    )

    for i in 1:length(result_list_r)
        df_r = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/resistance/" * result_list_r[i]))
        fitted_param_mean.condition[i] = result_list_r[i][16:19]
        fitted_param_mean.salinity[i] = salinities[parse(Int128, result_list_r[i][19])]
        fitted_param_mean.ρ[i] = mean(df_r.emerg)
        fitted_param_mean.Δ[i] = mean(df_r.pop)

        fitted_param_var.condition[i] = result_list_r[i][16:19]
        fitted_param_var.salinity[i] = salinities[parse(Int128, result_list_r[i][19])]
        fitted_param_var.ρ[i] = var(df_r.emerg)
        fitted_param_var.Δ[i] = var(df_r.pop)
    end
    mean_of_var = Dict()
    std_of_var = Dict()
    for i in ["A", "B", "C", "D"]
        temp = []
        temp = eachcol(filter(row -> row.condition[1:2] == i * "_", fitted_param_var))
        for j in 2:length(temp)
            mean_of_var[names(fitted_param_var)[j]] = mean(temp[j])
            std_of_var[names(fitted_param_var)[j]] = std(temp[j])
        end
    end

    return (fitted_param_mean, fitted_param_var, mean_of_var, std_of_var)
end
