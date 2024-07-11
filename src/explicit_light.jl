#Biofilm development with light integration
#Author : Chrislain Coubard
#Date 06/24

using PlotlyJS
using LinearAlgebra
using ColorSchemes
using DataFrames
using ColorTypes
using CSV
using Statistics
using DelimitedFiles
include("model_parameters.jl")
include("functions_model.jl")
include("functions_plots.jl")



tp = TimeParams()
hmp = HanModelParams()

##Space discretization
z = 1e-3
nz = 1000
dz = z/nz
xplt = 0:dz:z
light_intensities = [100,200,300,500,1000]

# pop_save = zeros(tp.n_save, nz)
# height_save = zeros(tp.n_save)
time_save = zeros(tp.n_save)
df_mu = DataFrame()
df_height = DataFrame(Height = Float64[], Intensity=String[], Time = Int64[])
df = DataFrame()


for I0 in light_intensities
    LI = zeros(nz)
    µ = zeros(nz)
    pop = zeros(nz)
    O2 = zeros(nz)
    height = zeros(tp.n_save)
    X0 = hmp.rho * dz
    pop[1:30] .= X0
    println("Start for $I0")
for time_step in 1:tp.n_save
    for i_inner in 1:tp.n_inner
        computelight!(LI, I0, hmp.ke, dz, nz, pop)
        updatemu!(µ, hmp.RD, hmp.RL, LI, hmp.n, hmp.Ik,hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, nz, pop)
        updateheight!(pop, µ, tp.dt)
        smootharray!(pop, nz, X0)

    end
    time = tp.dt*time_step*tp.n_inner/3600
    currheight = sum(pop) / hmp.rho
    time_save[time_step] = time
    height[time_step] = currheight
    push!(df_height, [currheight, "$(I0)", time])
    end
    colname = "$I0"
    df_mu[!,colname] = µ
    df[!,colname] = height
end


########## Plots ##############

pl = plot([
    scatter(x = time_save, y = df."100", mode = "markers", name = "100"), 
    scatter(x = time_save, y = df."200", mode = "markers", name = "200"),
    scatter(x = time_save, y = df."300", mode = "markers", name = "300"),
    scatter(x = time_save, y = df."500", mode = "markers", name = "500"),
    scatter(x = time_save, y = df."1000", mode = "markers", name = "1000"),
])
display(pl)

