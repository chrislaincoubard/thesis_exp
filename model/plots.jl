using PlotlyJS
using CSV
using Glob
using DataFrames
using Statistics
include("utils.jl")


function computelight!(arrLight, I0, ke, dz)
    
    # ind = findfirst(x -> x == 0, pop)
    for i in eachindex(arrLight)
        arrLight[i] = I0 * exp(-ke*(dz*i))
    end 
    arrLight[1] = I0
end

path = raw"C:\Users\Chrislain\Documents\Results\model_data_test"
# path = raw"C:\Users\Chrislain\Documents\Results\tmp\model_data_test_2"

# save_path_O2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_O2_plots")
# save_path_CO2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_CO2_plots_tests")
# save_path_N = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_N_plots")
# save_path_P = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_P_plots")
# # save_path_pH1 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_pH_plots\ph1")
# # save_path_pH2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_pH_plots\ph2")
# save_path_mu = mkpath(raw"c:\\Users\Chrislain\Documents\Results\tmp\plots\model_mu_test")

save_path_O2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_3\model_O2_plots")
save_path_CO2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_3\model_CO2_plots_tests")
save_path_CO2_50 = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_3\model_CO2_plots_tests_2")
save_path_N = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_3\model_N_plots")
save_path_P = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_3\model_P_plots")
# save_path_pH1 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_pH_plots\ph1")
# save_path_pH2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_pH_plots\ph2")
save_path_mu = mkpath(raw"c:\\Users\Chrislain\Documents\Results\plots_3\model_mu_test_2")

light = zeros(168)

computelight!(light, 100, 11650, 1e-6)

dfs = dataframesfromdir(path)
println(keys(dfs))

# df_O2 = dfs["model_O2_200.csv"]
# df_O2_100 = dfs["model_O2_100.csv"]
df_O2_50 = dfs["model_O2_50.csv"]
# df_CO2 = dfs["model_CO2_200.csv"]
df_CO2_50 = dfs["model_CO2_50.csv"]
# df_CO2_100 = dfs["model_CO2_100.csv"]
# df_N = dfs["model_N_200.csv"]
df_N_50 = dfs["model_N_50.csv"]
# df_P = dfs["model_P_200.csv"]
df_P_50 = dfs["model_P_50.csv"]
# df_mu = dfs["model_mu_200.csv"]
# df_pop_50 = dfs["model_pop_50.csv"]
# df_pop_100 = dfs["model_pop_100.csv"]
# df_pop_200 = dfs["model_pop_200.csv"]
# df_pop_300 = dfs["model_pop_300.csv"]
# df_pop_400 = dfs["model_pop_400.csv"]
# df_pop_800 = dfs["model_pop_800.csv"]

df_mu_50 = dfs["model_mu_50.csv"]
# df_mu_100 = dfs["model_mu_100.csv"]
# df_mu_200 = dfs["model_mu_200.csv"]
# df_mu_300 = dfs["model_mu_300.csv"]
# df_mu_400 = dfs["model_mu_400.csv"]
# df_mu_800 = dfs["model_mu_800.csv"]

# df_N_100 = dfs["model_N_100.csv"]
# df_P_100 = dfs["model_P_100.csv"]

# df_N_300 = dfs["model_N_300.csv"]
# df_P_300 = dfs["model_P_300.csv"]

# df_mu_gross_100 = dfs["model_mu_gross_100.csv"]
# df_R_100 = dfs["model_R_100.csv"]

# dfs_mu = [dfs["model_mu_50.csv"], dfs["model_mu_100.csv"], dfs["model_mu_200.csv"], dfs["model_mu_300.csv"], dfs["model_mu_400.csv"], dfs["model_mu_800.csv"]]
# df_H = dfs["model_H_200.csv"]
# df_H1 = df_H[:,r"^H1"]
# df_H2 = df_H[:,r"^H2"]
# df_pH = dfs["model_pH_200.csv"]
# df_pH1 = df_pH[:, r"^A_"]
# df_pH2 = df_pH[:, r"^B_"]

# height_50 = []

function plot_O2(df_O2)
    for time in names(df_O2)
        file_name_O2 = "O2_profile_$(time)_h.png"
        cleanO2 = removezeros(df_O2[!,time])
        p = plot([scatter(x = 0:1000, y = cleanO2, mode = "line", name = "O2")], 
            Layout(title = "02 concentration profile for 200 µmol<sub>photons</sup>/m<sup>2</sup>/s after $(time) hours",
            xaxis_title = "Depth (µm)",
            yaxis_title = "O2 concentration mol/m<sup>3</sup>",
            xaxis_range = [0,325],
            yaxis_range = [0.25,0.45]))
        # display(p)
        savefig(p, joinpath(save_path_O2, file_name_O2))
    end
end
# for time in names(df_CO2)
#     filename_CO2 = "CO2_profile_$(time)_h.png"
#     cleanCO2 = removezeros(df_CO2[!,time])
#     p = plot(scatter(x = 0:1000, y= cleanCO2, mode = "line"),
#         Layout(title = "CO2 concentration profile $(time) h",
#         xaxis_title = "Depth (µm)",
#         yaxis_title = "CO2 concentration mol/m3",
#         xaxis_range = [0,325],
#         yaxis_range = [-0.1,0.16]))
#     # display(p)
#     savefig(p, joinpath(save_path_CO2, filename_CO2))
# end

function plot_CO2(df_CO2)
    for time in names(df_CO2)
        filename_CO2 = "CO2_profile_$(time)_h.png"
        cleanCO2 = removezeros(df_CO2[!,time])
        p = plot(scatter(x = 0:1000, y= cleanCO2, mode = "line"),
            Layout(title = "CO2 concentration profile $(time) h",
            xaxis_title = "Depth (µm)",
            yaxis_title = "CO2 concentration mol/m3",
            xaxis_range = [0,325],
            yaxis_range = [-0.1,0.16]))
        # display(p)
        savefig(p, joinpath(save_path_CO2_50, filename_CO2))
    end
end

function plot_N(df_N)
    for time in names(df_N)
        filename_N = "NO3_profile_$(time)_h.png"
        cleanN = removezeros(df_N[!,time])
        p = plot(scatter(x = 0:1000, y= cleanN, mode = "line"),
            Layout(title = "NO3 concentration profile $(time) h",
            xaxis_title = "Depth (µm)",
            yaxis_title = "NO3 concentration mol/m3"))
        savefig(p, joinpath(save_path_N, filename_N))
    end
end

function plot_P(df_P)
    for time in names(df_P)
        filename_P = "P_profile_$(time)_h.png"
        cleanN = removezeros(df_P[!,time])
        p = plot(scatter(x = 0:1000, y= cleanN, mode = "line"),
            Layout(title = "P concentration profile $(time) h",
            xaxis_title = "Depth (µm)",
            yaxis_title = "P concentration mol/m3"))
        savefig(p, joinpath(save_path_P, filename_P))
    end
end

plot_O2(df_O2_50)
# plot_CO2(df_CO2_50)
plot_N(df_N_50)
plot_P(df_P_50)

# for time in names(df_pH1)
#     filename_pH = "pH1_profile_$(time)_h.png"
#     cleanN = removezeros(df_pH[!,"$time"])
#     p = plot(scatter(x = 0:1000, y= cleanN, mode = "line"),
#         Layout(title = "Ph1 profile $(time) h",
#         xaxis_title = "Depth (µm)",
#         yaxis_title = "pH"))
#     savefig(p, joinpath(save_path_pH1, filename_pH))
# end

# for time in names(df_mu)
#     filename_mu = "mu_profile_$(time)_h.png"
#     cleanmu_50 = removezeros(df_mu_50[!,"$time"])
#     cleanmu_100 = removezeros(df_mu_100[!,"$time"])
#     cleanmu_200 = removezeros(df_mu_200[!,"$time"])
#     cleanmu_300 = removezeros(df_mu_300[!,"$time"])
#     cleanmu_400 = removezeros(df_mu_400[!,"$time"])
#     cleanmu_800 = removezeros(df_mu_800[!,"$time"])
#     p = plot(
#         [scatter(x = 0:1000, y = cleanmu_50.*86400, mode = "line",name = "50",line = attr(width = 3)),
#         scatter(x = 0:1000, y = cleanmu_100.*86400, mode = "line",name = "100",line = attr(width = 3)),
#         scatter(x = 0:1000, y = cleanmu_200.*86400, mode = "line",name = "200",line = attr(width = 3)),
#         scatter(x = 0:1000, y = cleanmu_300.*86400, mode = "line",name = "300", line = attr(color = "black", width = 3)),
#         scatter(x = 0:1000, y = cleanmu_400.*86400, mode = "line",name = "400",line = attr(width = 3)),
#         scatter(x = 0:1000, y = cleanmu_800.*86400, mode = "line",name = "800",line = attr(width = 3)),],
#     Layout(title = "net mu profile $(time) h",
#     xaxis_title = "Depth (µm)",
#     yaxis_title = "net mu (d-1)",
#     legend_title_text = "Light intensities<br>(µmol<sub>photons</sup>/m<sup>2</sup>/s)"))
#     savefig(p, joinpath(save_path_mu, filename_mu))
#     if time == "168"
#         display(p)
#     end
# end

# for time in names(df_mu)
#     filename = "mu_details_$(time)_h.png"
#     clean_net_mu = removezeros(df_mu_100[!,"$time"])
#     clean_gross_mu = removezeros(df_mu_gross_100[!,"$time"])
#     clean_R = removezeros(df_R_100[!,"$time"])
#     p = plot([
#         scatter( x = 0:length(clean_net_mu), y = clean_net_mu.*86400, mode = "line", line = attr(width = 3), name = "Net µ"),
#         scatter( x = 0:length(clean_gross_mu), y = clean_gross_mu.*86400, mode = "line", line = attr(width = 3), name = "Gross µ"),
#         scatter( x = 0:length(clean_R), y = clean_R.*86400, mode = "line", line = attr(width = 3), name = "Respiration"),
#         scatter( x = 0:length(clean_net_mu), y = light, mode = "line", line = attr(dash = "dash", color = "black"), name = "light intensity",yaxis = "y2")],
#         Layout(title = "Growth rate details : 7 days, light 100 µmol<sub>photons</sup>/m<sup>2</sup>/s",
#         yaxis_title = "Growth and respiration rates (d<sup>-1</sup>)",
#         xaxis_title = "Depths (µm)",
#         yaxis2=attr(
#             title="Light intensity<br>(µmol<sub>photons</sup>/m<sup>2</sup>/s)",
#             overlaying="y", side="right")))
#     if time == "168"
#         display(p)
#     end
# end

# for time in names(df_N)
#     filename_N = "ions_profile_$(time)_h.png"
#     cleanN = removezeros(df_N_100[!,time])
#     cleanP = removezeros(df_P_100[!,time])
#     p = plot(
#         [scatter(x = 0:1000, y= cleanN, mode = "line", name = "N"),
#         scatter(x = 0:1000, y = cleanP, mode = "line", name = "P")],
#         Layout(title = "NO3 concentration profile $(time) h",
#         xaxis_title = "Depth (µm)",
#         yaxis_title = "NO3 concentration mol/m3"))
#     if time == "168"
#         display(p)
#     end
#     # savefig(p, joinpath(save_path_N, filename_N))
# end

# for time in names(df_N)
#     if time == "168"
#     filename = "Gases_profiles_$(time)_h.png"
#     cleanO2 = removezeros(df_O2_100[!,time])
#     cleanCO2 = removezeros(df_CO2_100[!,time])
#     cleanmu = removezeros(df_mu_100[!,time])
#     p = plot([
#         scatter(x = 0:1000, y = cleanO2, mode = "line", name = "O2", line=attr(width = 3)),
#         scatter(x = 0:1000, y = cleanCO2, mode = "line", name = "CO2", line=attr(width = 3)),
#         scatter(x = 0:1000, y = cleanmu.*86400, mode = "line", name = "mu", line=attr(width = 3), yaxis = "y2")
#     ],Layout(title = "Gases profiles $(time) h", xaxis_title = "Depths (µm)", yaxis_title = "Gases concentration (mM)",
#     yaxis2 = attr(title="Growth rate (d<sup>-1</sup>)",overlaying="y", side="right"), font = attr(size = 16)))
#     display(p)
#     end
# end

# for time in names(df_pH2)
#     filename_pH = "pH2_profile_$(time)_h.png"
#     cleanN = removezeros(df_pH[!,time])
#     p = plot(scatter(x = 0:1000, y= cleanN, mode = "line"),
#         Layout(title = "Ph2 profile $(time) h",
#         xaxis_title = "Depth (µm)",
#         yaxis_title = "pH"))
#     savefig(p, joinpath(save_path_pH2, filename_pH))
# end

println("End of script")