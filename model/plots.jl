using PlotlyJS
using CSV
using Glob
using DataFrames
using Statistics
include("utils.jl")


path = raw"C:\Users\Chrislain\Documents\Results\model_data"
save_path_O2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\model_O2_plots")
save_path_CO2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\model_CO2_plots")
save_path_N = mkpath(raw"C:\Users\Chrislain\Documents\Results\model_N_plots")



dfs = dataframesfromdir(path)

df_O2 = dfs["model_O2_200.csv"]
df_CO2 = dfs["model_CO2_200.csv"]
df_N = dfs["model_N_200.csv"]

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
    savefig(p, joinpath(save_path, file_name_O2))
end

for time in names(df_CO2)
    filename_CO2 = "CO2_profile_$(time)_h.png"
    cleanCO2 = removezeros(df_CO2[!,time])
    p = plot(scatter(x = zplt*10^6, y= cleanCO2, mode = "line"),
        Layout(title = "CO2 concentration profile $(time) h",
        xaxis_title = "Depth (µm)",
        yaxis_title = "CO2 concentration mol/m3"))
    display(p)
    savefig(p, joinpath(save_path, filename_CO2))
end

for time in names(df_N)
    filename_N = "NO3_profile_$(time)_h.png"
    cleanN = removezeros(df_N[!,time])
    p = plot(scatter(x = zplt*10^6, y= cleanN, mode = "line"),
        Layout(title = "NO3 concentration profile $(time) h",
        xaxis_title = "Depth (µm)",
        yaxis_title = "NO3 concentration mol/m3"))
    savefig(p, joinpath(save_path, filename_N))
end