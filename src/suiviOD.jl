using PlotlyJS
using CSV
using DataFrames

data = CSV.File(raw"C:\Users\chris\OneDrive\Documents\thèse\results\suivicroissance.csv", delim =";" ) |> DataFrame

data50 = filter(:Intensite => i -> i == 50, data)
data100 = filter(:Intensite => i -> i == 100, data)
println(data50)

function plotOD(data, lambda, intensite)
    p = plot(scatter( x = data.Temps, y= data[!,lambda], mode = "line"),
    Layout(title = "I$intensite OD$lambda", xaxis_title = "time (h)", yaxis_title="OD"))
    display(p)
end

function plotCells(data, intensite)
    p = plot(scatter(x = data.Temps, y = data.Cells, mode = "line"),
    Layout(title ="I$intensite Cell Number", xaxis_title = "time (h)", yaxis_title = "Cells/mL"))
    display(p)
end

plotOD(data50, "678", 50)
plotOD(data50, "750",50)
plotOD(data50,"800",50)
plotOD(data50, "540",50)
plotOD(data50, "600",50)
plotOD(data50, "480",50)
plotOD(data50, "530",50)

plotOD(data100, "678",100)
plotOD(data100, "750",100)
plotOD(data100, "800",100)
plotOD(data100, "540",100)
plotOD(data100, "600",100)
plotOD(data100, "480",100)
plotOD(data100, "530",100)
plotCells(data50, 50)
plotCells(data100, 100)





# p = plot(scatter( x = data50.Temps, y= data50."678", mode = "line"),
# Layout(title = "I50 OD678", xaxis_title = "time (h)", yaxis_title="OD"))
# display(p)
# p1 = plot(scatter( x = data50.Temps, y= data50."750", mode = "line"),
# Layout(title = "I50 OD750", xaxis_title = "time (h)", yaxis_title="OD"))
# display(p1)
# p2 = plot(scatter( x = data50.Temps, y= data50."750", mode = "line"),
# Layout(title = "I50 OD750", xaxis_title = "time (h)", yaxis_title="OD"))
# display(p1)
µ1 = (log(0.42834) - log(0.1257))/(72/24-4/24)
µ2 = (log(0.33032) - log(0.10822))/(139/24-72/24)
µ3 = (log(0.3027) - log( 0.13094))/(185/24-139/24)
µ4 = (log(0.556) - log( 0.205))/(72/24-4/24)
µ5 = (log(0.397) - log( 0.125))/(139/24-72/24)
µ6 = (log(0.404) - log( 0.158))/(185/24-139/24)
µ7 = (log(0.613) - log(0.404))/(210/24-185/24)

µµ = (log(0.2445) - log(0.0738))/(72/24-4/24)

µ = (log())

println("taux de croissance pour I50 :\n Batch 1 : $µ1\n Batch 2 : $µ2\n Batch 3 : $µ3")

println("taux de croissance pour I100 :\n Batch 1 : $µ4\n Batch 2 : $µ5\n Batch 3 : $µ6\n")

println(µ7)