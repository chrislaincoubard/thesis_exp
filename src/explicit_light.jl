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

function updateO2!(O, dt, dz, D,pop,O2atm)
    ind = findfirst(x -> x == 0, pop)-1
    O[1] = ((O[2]-O2atm)*D*2/dz)*dt/dz + O[1]
    for i in 2:ind-1
        O[i] = ((O[i+1]-O[i])*D/dz+(O[i]-O[i-1])*D/dz)*dt/dz + O[i]
    end
    O[ind] = ((O[ind]-O[ind-1])*D/(dz/2))*dt/dz + O[ind]
end

tp = TimeParams()
hmp = HanModelParams()
gmp = GasesParams()

##Space discretization
z = 1e-3
nz = 1000
dz = z/nz
xplt = 0:dz:z
# light_intensities = [100,200,300,500,1000]
light_intensities = [300]

# pop_save = zeros(tp.n_save, nz)
# height_save = zeros(tp.n_save)
time_save = zeros(tp.n_save)
df_mu = DataFrame()
df_height = DataFrame(Height = Float64[], Intensity=String[], Time = Float64[])
df = DataFrame()
df_O2 = DataFrame()

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
    println("we're at $time_step")
    for i_inner in 1:tp.n_inner
        computelight!(LI, I0, hmp.ke, dz, pop)
        updatemu!(µ, hmp.RD, hmp.RL, LI, hmp.n, hmp.Ik,hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, pop)
        updateheight!(pop, µ, tp.dt)
        smootharray!(pop, nz, X0)
        updateO2!(O2, tp.dt, dz, gmp.D_oxygen, pop, gmp.O2atm)
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
    df_O2[!,colname] = O2
end
println("Computation Done")


########## Plots ##############

pl = plot([
    scatter(x = time_save, y = df."300", mode = "markers", name = "300"),
])
display(pl)

function buildO2traces(df, x)
    traces = Vector{GenericTrace}(undef, length(names(df)))
    for (index, n) in enumerate(names(df))
    dfcol_clean = filter(x -> x!= 0, df[!,n])
    sct = scatter(y = dfcol_clean, x = x, mode = "markers",
     marker = attr(size = 3), name = n)
    traces[index] = sct
    end
    return traces
end

O2_traces = buildO2traces(df_O2, zplt)

plot(O2_traces, Layout(title = "test"))
# print(first(df_O2, 15))
# pl1 = plot([
#     scatter(x = time_save, y = df_O2."100", mode = "markers", name = "100"), 
#     scatter(x = time_save, y = df_O2."200", mode = "markers", name = "200"),
#     scatter(x = time_save, y = df_O2."300", mode = "markers", name = "300"),
#     scatter(x = time_save, y = df_O2."500", mode = "markers", name = "500"),
#     scatter(x = time_save, y = df_O2."1000", mode = "markers", name = "1000"),
# ])
# display(pl1)
