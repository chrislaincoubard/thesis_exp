using PlotlyJS
using CSV
using Glob
using DataFrames
using Statistics
include("utils.jl")


path = raw"C:\Users\Chrislain\Documents\Results\model_data_test"
# path = raw"C:\Users\Chrislain\Documents\Results\tmp\model_data_test_2"

# save_path_O2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_O2_plots")
# save_path_CO2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_CO2_plots_tests")
# save_path_N = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_N_plots")
# save_path_P = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_P_plots")
# # save_path_pH1 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_pH_plots\ph1")
# # save_path_pH2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_pH_plots\ph2")
# save_path_mu = mkpath(raw"c:\\Users\Chrislain\Documents\Results\tmp\plots\model_mu_test")

save_path_O2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_2\model_O2_plots")
save_path_CO2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_2\model_CO2_plots_tests")
save_path_N = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_2\model_N_plots")
save_path_P = mkpath(raw"C:\Users\Chrislain\Documents\Results\plots_2\model_P_plots")
# save_path_pH1 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_pH_plots\ph1")
# save_path_pH2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\tmp\plots\model_pH_plots\ph2")
save_path_mu = mkpath(raw"c:\\Users\Chrislain\Documents\Results\plots_2\model_mu_test")



dfs = dataframesfromdir(path)

df_O2 = dfs["model_O2_200.csv"]
df_CO2 = dfs["model_CO2_200.csv"]
df_N = dfs["model_N_200.csv"]
df_P = dfs["model_P_200.csv"]
# df_H = dfs["model_H_200.csv"]
# df_H1 = df_H[:,r"^H1"]
# df_H2 = df_H[:,r"^H2"]
# df_pH = dfs["model_pH_200.csv"]
# df_pH1 = df_pH[:, r"^A_"]
# df_pH2 = df_pH[:, r"^B_"]
df_mu = dfs["model_mu_200.csv"]



for time in names(df_O2)
    file_name_O2 = "O2_profile_$(time)_h.png"
    # cleanO2 = removezeros(df_O2[!,time])
    p = plot([scatter(x = 0:1000, y = df_O2[!,time], mode = "line", name = "O2")], 
        Layout(title = "02 concentration profile for 200 µmol<sub>photons</sup>/m<sup>2</sup>/s after $(time) hours",
        xaxis_title = "Depth (µm)",
        yaxis_title = "O2 concentration mol/m<sup>3</sup>",
        xaxis_range = [0,325],
        yaxis_range = [0.25,0.45]))
    # display(p)
    savefig(p, joinpath(save_path_O2, file_name_O2))
end

for time in names(df_CO2)
    filename_CO2 = "CO2_profile_$(time)_h.png"
    # cleanCO2 = removezeros(df_CO2[!,time])
    p = plot(scatter(x = 0:1000, y= df_CO2[!,time], mode = "line"),
        Layout(title = "CO2 concentration profile $(time) h",
        xaxis_title = "Depth (µm)",
        yaxis_title = "CO2 concentration mol/m3",
        xaxis_range = [0,325]))
    # display(p)
    savefig(p, joinpath(save_path_CO2, filename_CO2))
end

# for time in names(df_N)
#     filename_N = "NO3_profile_$(time)_h.png"
#     cleanN = removezeros(df_N[!,time])
#     p = plot(scatter(x = 0:1000, y= cleanN, mode = "line"),
#         Layout(title = "NO3 concentration profile $(time) h",
#         xaxis_title = "Depth (µm)",
#         yaxis_title = "NO3 concentration mol/m3"))
#     savefig(p, joinpath(save_path_N, filename_N))
# end

# for time in names(df_P)
#     filename_P = "P_profile_$(time)_h.png"
#     cleanN = removezeros(df_N[!,time])
#     p = plot(scatter(x = 0:1000, y= cleanN, mode = "line"),
#         Layout(title = "P concentration profile $(time) h",
#         xaxis_title = "Depth (µm)",
#         yaxis_title = "P concentration mol/m3"))
#     savefig(p, joinpath(save_path_P, filename_P))
# end

# for time in names(df_pH1)
#     filename_pH = "pH1_profile_$(time)_h.png"
#     cleanN = removezeros(df_pH[!,"$time"])
#     p = plot(scatter(x = 0:1000, y= cleanN, mode = "line"),
#         Layout(title = "Ph1 profile $(time) h",
#         xaxis_title = "Depth (µm)",
#         yaxis_title = "pH"))
#     savefig(p, joinpath(save_path_pH1, filename_pH))
# end

for time in names(df_mu)
    filename_mu = "mu_profile_$(time)_h.png"
    # cleanmu = removezeros(df_mu[!,"$time"])
    p = plot(scatter(x = 0:1000, y = df_mu[!,"$time"].*86400, mode = "line"),
    Layout(title = "Mu profile $(time) h",
    xaxis_title = "Depth (µm)",
    yaxis_title = "mu (d-1)"))
    savefig(p, joinpath(save_path_mu, filename_mu))
end

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