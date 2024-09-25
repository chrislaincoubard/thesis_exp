using PlotlyJS
using CSV
using Glob
using DataFrames
using Statistics
include("utils.jl")

path = raw"C:\Users\chris\OneDrive\Documents\thèse\results"
dfs = dataframesfromdir(path)
df = dfs["data_model_300.csv"]
dff = dfs["data_model_100.csv"]
df_height_300 = dfs["height_model_300_7.csv"]
df_height_100 = dfs["height_model_100_7.csv"]
df_height_500 = dfs["height_model_500_7.csv"]
df_height_800 = dfs["height_model_800_7.csv"]
df_integral_300 = dfs["data_full_integral_300.csv"]
df_integral_100 = dfs["data_full_integral_100.csv"]
df_integral_500 = dfs["data_full_integral_500.csv"]
filter!(row -> !(row.mu_net == 0), df)

zplt = 1:nrow(df)
display(df)
pltR = plot([
    scatter(x = zplt, y = df[!,"R"].*86400, mode = "line", name = "Respiration"),
    scatter(x = zplt, y = df[!,"gross_mu"].*86400, mode = "line", name = "µ brut"),
    scatter(x = zplt, y = df[!,"mu_net"].*86400, mode = "line", name = "µ net"),
    scatter(x = zplt, y = df[!,"light"], mode = "line", line =attr(dash="dash", color = "red"), name = "light", yaxis = "y2")],
    Layout(title = "Growth rate", xaxis_title = "Depth of biofilm (µm)", yaxis_title = "s-1",
    yaxis2 = attr(title = "Light (µmol/m2/s)",overlaying = "y", side ="right"),
    legend=attr(x=0.7, y=1,),
    modebar_add = ["drawline", "drawopenpath", "eraseshape"]))

pltRR = plot([
scatter(x = zplt, y = dff[!,"R"].*86400, mode = "line", name = "Respiration"),
scatter(x = zplt, y = dff[!,"gross_mu"].*86400, mode = "line", name = "µ brut"),
scatter(x = zplt, y = dff[!,"mu_net"].*86400, mode = "line", name = "µ net"),
scatter(x = zplt, y = dff[!,"light"], mode = "line", line =attr(dash="dash", color = "red"), name = "light", yaxis = "y2")],
Layout(title = "Growth rate", xaxis_title = "Depth of biofilm (µm)", yaxis_title = "s-1",
yaxis2 = attr(title = "Light (µmol/m2/s)",overlaying = "y", side ="right"),
legend=attr(x=0.855, y=1,),
modebar_add = ["drawline", "drawopenpath", "eraseshape"]))

plt_all_height = plot([
    scatter(x = df_height_300[!,"time"], y = df_height_300[!,"Height"].*10^6, mode = "line", name = "300"),
    scatter(x = df_height_100[!,"time"], y = df_height_100[!,"Height"].*10^6, mode = "line", name = "100"),
    scatter(x = df_height_500[!,"time"], y = df_height_500[!,"Height"].*10^6, mode = "line", name = "500"),
    scatter(x = df_height_800[!,"time"], y = df_height_800[!,"Height"].*10^6, mode = "line", name = "800", line=attr(color = "black"))],
    Layout(title = "Growth for different light", xaxis_title ="Time (h)", yaxis_title = "Height of biofilm (µm)",
    modebar_add = ["drawline", "drawopenpath", "eraseshape"])
)
display(plt_all_height)


pltheight = plot([
    scatter(x = df_height_300[!,"time"], y = df_height_300[!,"Height"].*10^6, mode = "line", name ="FVM"),
    scatter(x = df_integral_300[!,"time"], y = df_integral_300[!,"height"].*10^6, mode = "line", name = "integral")], 
    Layout(title = "Biofilm Growth, light = 300", xaxis_title = "Time (h)", yaxis_title = "Height (µm)",
    modebar_add = ["drawline", "drawopenpath", "eraseshape"]))
display(pltheight)
pltheight2 = plot([
    scatter(x = df_height_100[!,"time"], y = df_height_100[!,"Height"].*10^6, mode = "line", name ="FVM"),
    scatter(x = df_integral_100[!,"time"], y = df_integral_100[!,"height"].*10^6, mode = "line", name = "integral")], 
    Layout(title = "Biofilm Growth, light = 100", xaxis_title = "Time (h)", yaxis_title = "Height (µm)",
    modebar_add = ["drawline", "drawopenpath", "eraseshape"]))
display(pltheight2)
pltheight3 = plot([
    scatter(x = df_height_500[!,"time"], y = df_height_500[!,"Height"].*10^6, mode = "line", name ="FVM"),
    scatter(x = df_integral_500[!,"time"], y = df_integral_500[!,"height"].*10^6, mode = "line", name = "integral")], 
    Layout(title = "Biofilm Growth, light = 500", xaxis_title = "Time (h)", yaxis_title = "Height (µm)",
    modebar_add = ["drawline", "drawopenpath", "eraseshape"]))
display(pltheight3)
# display(pltLight)

# plt_mean = plot([
#     scatter(x = df_height_300[!,"time"], y = df_height_300[!,"mean_mu"].*86400, mode = "markers"),
#     scatter(x = df_integral_300[!,"time"], y = df_integral_300[!,"mu"].*86400, mode = "markers")],
# Layout(title = "mean mu", xaxis_title = "Time (h)", yaxis_title = "Mean mu (d-1)"))
# pltOxygen = plot(scatter(x = zplt*10^6, y = clean_oxygen, mode = "markers"), 
# Layout(title = "Oxygen", xaxis_title = "Depth (µm)", yaxis_title = "O2 concentration in mM"))

display(pltR)
display(pltRR)
# display(plt_mean)
savefig(pltR, raw"C:\Users\chris\OneDrive\Documents\thèse\plots\growthrate300.png")
savefig(plt_all_height, raw"C:\Users\chris\OneDrive\Documents\thèse\plots\different_light_growth.png")
savefig(pltheight, raw"C:\Users\chris\OneDrive\Documents\thèse\plots\FVM_vs_Integration_300.png")
savefig(pltheight2, raw"C:\Users\chris\OneDrive\Documents\thèse\plots\FVM_vs_Integration_100.png")
savefig(pltheight3, raw"C:\Users\chris\OneDrive\Documents\thèse\plots\FVM_vs_Integration_500.png")


println("finished")