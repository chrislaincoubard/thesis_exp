using PlotlyJS
using CSV
using DataFrames
using Statistics

data = CSV.File(raw"C:\Users\Chrislain\Documents\exp_data\rate_photo_Bminus.csv") |> DataFrame
data2 = CSV.File(raw"C:\Users\Chrislain\Documents\exp_data\rate_photoB+.csv") |> DataFrame

data_control = CSV.File(raw"C:\Users\Chrislain\Documents\exp_data\control_Oxygen.csv",delim = ";") |> DataFrame
data_025B = CSV.File(raw"C:\Users\Chrislain\Documents\exp_data\060924_0,5MBicarbonate.csv", header =75, silencewarnings = true, footerskip = 2) |> DataFrame
data_05B = CSV.File(raw"C:\Users\Chrislain\Documents\exp_data\060924_0,25MBicarbonate.csv", header = 75, silencewarnings = true, footerskip = 2) |> DataFrame


println(first(data_025B,5))
rates_control = Float64[]
rates_b025 = Float64[]
rates_b05 = Float64[]

println(data_025B."Oxygen 1"[12])

for i in eachindex(data_05B."Oxygen 1")[Not(1:30)]
    rate = (data_05B."Oxygen 1"[i] - data_05B."Oxygen 1"[i-30])*2
    push!(rates_b05, rate)
end

for i in eachindex(data_025B."Oxygen 1")[Not(1:30)]
    rate = (data_025B."Oxygen 1"[i] - data_025B."Oxygen 1"[i-30])*2
    push!(rates_b025, rate)
end

for i in eachindex(data_control.Oxygen)[Not(1:30)]
    rate = (data_control.Oxygen[i] - data_control.Oxygen[i-30])*2
    push!(rates_control, rate) 
end

println(rates_b05)

p = plot([
    scatter(x = eachindex(rates_control), y = rates_control, mode = "line", name = "Control"),
    scatter(x = eachindex(rates_b025), y = rates_b025, mode = "line", name = "Bicarb 0.25"),
    scatter(x = eachindex(rates_b05), y = rates_b05, mode = "line", name = "Bicarb 0.5")
    ],
Layout(title = "rate O2 control", xaxis_title = "time", yaxis_title = "O2"))
display(p)

df25 = DataFrame()
df50 = DataFrame()

# df25[!,"0"] = rates_b025[1:301]
# df25[!,"20"] = rates_b025[302:601]
# df25[!,"50"] = rates_b025[602:901]
# df25[!,"100"] = rates_b025[902:1201]
# df25[!,"180"] = rates_b025[1202:1501]
# df25[!,"250"] = rates_b025[1502:1801]
# df25[!,"500"] = rates_b025[1802:2101]
# df25[!,"1000"] = rates_b025[2102:2401]

# df50[!,"0"] = rates_b05[1:301]
# df50[!,"20"] = rates_b05[302:601]
# df50[!,"50"] = rates_b05[602:901]
# df50[!,""]

means = mean.(eachcol(data))
light = parse.(Int64, names(data))
means2 = mean.(eachcol(data2))
println(data2)
mean0 = mean(skipmissing(data[!,"0"]))
means[1] = mean0
mean02 = mean(skipmissing(data2[!,"0"]))
means2[1] = mean02
display(means)
means .= means ./ 1000 .*60 ./ 0.00159
means2 .= means2 ./ 1000 .*60 ./ 0.00198
plot(
    [scatter(x = light,y = means, name = "bicarbonate -"),
    scatter(x = light, y = means2, name = "bicarbonate +")],
    Layout(title = "PI-curve",xaxis_title = "light µmol/m2/s", yaxis_title = "Oxygen rate (µmolO2/mgChl/h)",
    modebar_add = ["drawopenpath", "eraseshape"]))

    