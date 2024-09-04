using PlotlyJS
using CSV
using Glob
using DataFrames
using Statistics
include("utils.jl")

path = raw"C:\Users\LGPM Bamako\Documents\Results\result_model"
dfs = dataframesfromdir(path)
df = dfs["data_model_300.csv"]
df_height = dfs["height_model_300.csv"]
df_integral = dfs["data_full_integral_300.csv"]
df_mu_1s = dfs["data_model_300_1s.csv"]
df_mu_150s = dfs["data_model_300_150s.csv"]
df_mu_1440 = dfs["data_model_300_1444s.csv"]
filter!(row -> !(row.mu_net == 0), df)

zplt = 1:nrow(df)

pltR = plot([
    scatter(x = zplt, y = df[!,"R"].*86400, mode = "line", name = "Respiration"),
    scatter(x = zplt, y = df[!,"gross_mu"].*86400, mode = "line", name = "µ brut"),
    scatter(x = zplt, y = df[!,"mu_net"].*86400, mode = "line", name = "µ net")],
    Layout(title = "Growth rate", xaxis_title = "Depth of biofilm (µm)", yaxis_title = "s-1"))


pltLight = plot(scatter(x = zplt, y = df[!,"light"], mode = "markers"), 
Layout(title = "Light", xaxis_title = "Depth (µm)", yaxis_title = "PPFD (µmol/m²/s)"))


pltheight = plot([
    scatter(x = df_height[!,"time"], y = df_height[!,"Height"].*10^6, mode = "line"),
    scatter(x = df_integral[!,"time"], y = df_integral[!,"height"].*10^6, mode = "line")], 
    Layout(title = "Biofilm Growth", xaxis_title = "Time (h)", yaxis_title = "Height (µm)"))
display(pltheight)
display(pltLight)

plt_mean = plot([
    scatter(x = df_height[!,"time"], y = df_height[!,"mean_mu"].*86400, mode = "markers"),
    scatter(x = df_integral[!,"time"], y = df_integral[!,"mu"].*86400, mode = "markers")],
Layout(title = "mean mu", xaxis_title = "Time (h)", yaxis_title = "Mean mu (d-1)"))
# pltOxygen = plot(scatter(x = zplt*10^6, y = clean_oxygen, mode = "markers"), 
# Layout(title = "Oxygen", xaxis_title = "Depth (µm)", yaxis_title = "O2 concentration in mM"))

display(pltR)
display(plt_mean)


plotcomp = plot([
    scatter(x = zplt, y = df[!,"mu_net"], mode = "markers", name = "50s"),
    scatter(x = zplt, y = df_mu_1s[!,"mu_net"], mode = "markers", name = "1s"),
    scatter(x = zplt, y = df_mu_150s[!,"mu_net"], mode = "markers", name = "150s"),
    scatter(x = zplt, y = df_mu_1440[!,"mu_net"], mode = "markers", name = "1440s")],
    Layout(title = "test different DeltaT", xaxis_title="dephts (µm)", yaxis_title="net mu"))

display(plotcomp)

# for i in eachindex(df[!,"mu_net"])
#     if df_mu_1s[!,"mu_net"][i] != df_mu_150s[!,"mu_net"]
#         bit = df_mu_1s[!,"mu_net"][i] - df_mu_150s[!,"mu_net"][i]
#         println(bit)
#         # println("1s = $(df_mu_1s[!,"mu_net"][i]) vs 150s = $(df_mu_150s[!,"mu_net"][i])")
#     end
# end
println(isequal(df_mu_1s[!,"mu_net"], df_mu_150s[!,"mu_net"]))

println("finished")