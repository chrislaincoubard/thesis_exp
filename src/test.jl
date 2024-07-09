using PlotlyJS
using Statistics
using LinearAlgebra
using ColorSchemes
using DataFrames
include("model_parameters.jl")

mp = ModelParams() 
light_intensities = 0:5:1000
length(light_intensities)
z = 1e-4
nz = 10000
dz = z/nz
zplt = 0:dz:z
df_mu = DataFrame()


RR = Float64[]
µµ = Float64[]
µgg = Float64[]
for I in light_intensities
    R = mp.RD + (mp.RL - mp.RD) * (I^mp.n/(I^mp.n + mp.Ik^mp.n))
    a = mp.k*mp.sigma*I
    b = mp.tau*mp.sigma*I
    c = mp.kd*mp.tau*mp.sigma^2*I^2/mp.kr
    µg = a / (1 + b +c)
    µ = (µg - R)
    push!(RR, R)
    push!(µµ, µ)
    push!(µgg, µg)
end
pl1 = plot(scatter(x = light_intensities, y = µgg*86400, mode = "line"), Layout(title = "µgg"))

pl2 = plot(scatter(x = light_intensities, y = RR, mode = "line"), Layout(title = "RR"))
pl3 = plot(scatter(x = light_intensities, y = µµ*86400, mode = "line"),  Layout(title = "µµ"))
# pl4 = plot(scatter(x = light_intensities, y = µµC, mode = "line"),  Layout(title = "µµC"))
# pl5 = plot(scatter(x = light_intensities, y = µµCG, mode = "line"),  Layout(title = "µµCG"))
display(pl1)
display(pl2)
display(pl3)


# a = mp.k * mp.sigma
# b = mp.tau * mp.sigma
# c = (mp.kd/mp.kr)*mp.tau*mp.sigma^2

# 4.9e-6
# 4.9e-6 == 0.0000049
# mp.k == 4.9e-6

# (mp.kd/mp.kr)*mp.tau*mp.sigma^2*168
# (mp.kd/mp.kr)*mp.tau*mp.sigma*mp.sigma*200
# mp.kd*mp.tau*mp.sigma^2*200/mp.kr

# I = 168
# R = mp.RD + (mp.RL - mp.RD) * (I^mp.n/(I^mp.n + mp.Ik^mp.n))
# a = mp.k*mp.sigma*I
# b = mp.tau*mp.sigma*I
# c = mp.kd*mp.tau*mp.sigma^2*I^2/mp.kr
# µg = a / (1 + b +c)
# µ = (µg - R)
# µ*86400/2