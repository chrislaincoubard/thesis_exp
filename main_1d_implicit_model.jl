#1D simulation of biofilm development (method : Finite Volume, implicit)
#Author : Chrislain Coubard (adapted from Patrick Perre)
#Date : 05/2024

using PlotlyJS
using LinearAlgebra
include("parameters.jl")
include("functions_implicits.jl")

sp = SpaceParams()
tp = TimeParams()
mp = ModelParams()


pop_save = zeros(tp.n_save)
time_save = zeros(tp.n_save)
glu_save = zeros(tp.n_save, sp.n_tot)
time = zeros(tp.n_save)

deltax_glu = fill(sp.L_glu/sp.n_glu, sp.n_glu)  # size of CVs regular mesh, full CVs at boundary, uniform by zone
deltax_pop = fill(sp.L_pop/sp.n_pop, sp.n_pop)
deltax = vcat(deltax_glu, deltax_pop)
x = zeros(sp.n_tot)
dx = zeros(sp.n_tot-1)
D = fill(mp.D_gl, sp.n_tot)

for i in 1:sp.n_glu
    x[i] = deltax_glu[1]*(i-0.5)
end
for i in 1:sp.n_pop
    x[i+sp.n_glu] = sp.L_glu + deltax_pop[1] * (i-0.5)
end
for i in 1:(sp.n_tot-1)
    dx[i] = (deltax[i] + deltax[i+1])/2
end

glu = zeros(sp.n_tot)
d_glu = zeros(sp.n_tot)
pop = zeros(sp.n_tot)
s_glu = zeros(sp.n_tot)

glu[1:sp.n_glu] .= 10
ind = [1,2,3,4,5] .+ sp.n_glu 
pop[ind] .= 1

# 1D simulation
diag = zeros(sp.n_tot)
low = zeros(sp.n_tot-1)
up = zeros(sp.n_tot-1)
B = zeros(sp.n_tot)


for time_steps in 1:tp.n_save
    for i_inner in 1:tp.n_inner
        UpdateValues!(sp.n_glu, sp.n_tot, mp.mu_0, mp.k_s, tp.dt, mp.ratio, mp.rho_pop, glu, pop, s_glu)
        SmoothArr!(sp.n_glu, sp.n_tot, pop)
        UpdateD!(sp.n_glu, sp.n_tot, mp.D_gl, D, pop)
        tridiag = MakeTriDiag(D, sp.n_tot, dx, tp.dt, deltax, up, diag, low)
        UpdateGlu!(B, glu, s_glu, tridiag)
    end
    glu_save[time_steps,:] = glu
    pop_save[time_steps] = sum(pop)
    time_save[time_steps] = tp.dt*time_steps*tp.n_inner/3600
end

plt1 = glucose_graph(tp.n_save, tp.n_inner, tp.dt, x, glu_save)
plt2 = population_graph(time_save, pop_save)
display(plt1)
display(plt2)