function get_mean_control(result_list)

    df = []
    fitted_param = DataFrame(condition=(""^length(result_list)), r=(zeros(length(result_list))), K=(zeros(length(result_list))))

    for i in 1:length(result_list)
        df = DataFrame(CSV.File("./MAIN_EXPERIMENT/BASAL_MODEL/results/control/" * result_list[i]))
        #df = DataFrame(CSV.File("./results/control/" * result_list[i]))
        fitted_param.condition[i] = result_list[i][23:26]
        fitted_param.r[i] = mean(df.r)
        fitted_param.K[i] = mean(df.K)
    end
    return (fitted_param)
end
