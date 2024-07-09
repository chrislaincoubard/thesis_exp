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
mp = ModelParams() 
light_intensities = [100,200,300,500,1000]
z = 1e-3
nz = 1000
dz = z/nz
zplt = 0:dz:z

µ_save = zeros(tp.n_save, nz)
time_save = zeros(tp.n_save)
pop_save = zeros(tp.n_save, nz)
height_save = zeros(tp.n_save)

df_mu = DataFrame()
df_height = DataFrame(Height = Float64[], Intensity=String[], Time = Int64[])
df = DataFrame()

for I0 in light_intensities
    LI = zeros(nz)
    µ = zeros(nz)
    pop = zeros(nz)
    height = zeros(tp.n_save)
    X0 = mp.rho * dz
    pop[1:30] .= X0
    println("Start for $I0")
    for time_step in 1:tp.n_save
        for i_inner in 1:tp.n_inner
            computelight!(LI, I0, mp.ke, dz, nz, pop)
            updatemu!(µ, mp.RD, mp.RL, LI, mp.n, mp.Ik,mp.k, mp.sigma, mp.tau, mp.kd, mp.kr, nz, pop)
            pop .= solvematrix(µ, tp.dt, nz, pop)
            smootharray!(pop, nz, X0)
        end
        time = tp.dt*time_step*tp.n_inner/3600
        currheight = sum(pop) / mp.rho
        time_save[time_step] = time
        height[time_step] = currheight
        push!(df_height, [currheight, "$I0", time])
    end
    colname = "$I0"
    df_mu[!,colname] = µ
    df[!,colname] = height
end

######## Plots #############
# println(first(df_mu,5))
mu_traces = buildmutraces(df_mu, zplt)
plotmutraces(mu_traces, "C:/Users/LGPM Bamako/Documents/Results", "mu_plot.png")
# writedlm( "height.csv",  df_height, ',')
println(first(df,5))

pl = plot([
    scatter(x = time_save, y = df."100", mode = "markers", name = "100"), 
    scatter(x = time_save, y = df."200", mode = "markers", name = "200"),
    scatter(x = time_save, y = df."300", mode = "markers", name = "300"),
    scatter(x = time_save, y = df."500", mode = "markers", name = "500"),
    scatter(x = time_save, y = df."1000", mode = "markers", name = "1000"),
])
display(pl)

