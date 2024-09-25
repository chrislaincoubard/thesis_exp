using PlotlyJS
using CSV
using DataFrames
using Statistics

OD190 = CSV.File(raw"C:\Users\chris\OneDrive\Documents\thesis_code\src\OD190.csv", header = 28, footerskip =2) |> DataFrame
OD415 = CSV.File(raw"C:\Users\chris\OneDrive\Documents\thesis_code\src\OD415.csv", header = 28, footerskip =2) |> DataFrame
OD643 = CSV.File(raw"C:\Users\chris\OneDrive\Documents\thesis_code\src\OD643.csv", header = 28, footerskip =2) |> DataFrame
OD862 = CSV.File(raw"C:\Users\chris\OneDrive\Documents\thesis_code\src\OD862.csv", header = 28, footerskip =2) |> DataFrame
OD923 = CSV.File(raw"C:\Users\chris\OneDrive\Documents\thesis_code\src\OD923.csv", header = 28, footerskip =2) |> DataFrame


pl = plot([
    scatter(x = OD190.Time./60, y = OD190."Oxygen 1", mode = "line", name = "OD 190"),
    scatter(x = OD415.Time./60, y = OD415."Oxygen 1", mode = "line", name = "OD 415"),
    scatter(x = OD643.Time./60, y = OD643."Oxygen 1", mode = "line", name = "OD 643"),
    scatter(x = OD862.Time./60, y = OD862."Oxygen 1", mode = "line", name = "OD 862"),
    scatter(x = OD923.Time./60, y = OD923."Oxygen 1", mode = "line", name = "OD 923") 
], Layout(title = "O2 production at light = 100", xaxis_title = "time (min)", yaxis_title = "O2 (nmol/mL)"))

rate_190 = Float64[]
println(length(rate_190))
OD190."Rate 1"[OD190."Rate 1".=="INVALID"] .= "0"
OD190."Rate 1" = parse.(Float64, OD190."Rate 1")
println(OD190."Rate 1")
for i in eachindex(OD190.Time)[Not(1:30)]
    inst_rate = (OD190."Oxygen 1"[i] - OD190."Oxygen 1"[i-29])*2
    push!(rate_190, inst_rate)
end

pp = plot([
    scatter(x = OD190.Time./60, y = OD190."Rate 1", mode = "line", name = "auto rate"),
    scatter(x = OD190.Time./60, y = rate_190, mode = "line", name = "manual rate")
])



