using LinearAlgebra
using PlotlyJS
using Distributions
include("model_parameters.jl")
include("functions_model.jl")

function O2source(mu, VO2x, mx, dz)
    source = zeros(ind)
    source .= 0
    return source
end

function compB(O2, O2sat, D, dz, dt)
    B = zeros(length(O2))
    for i in eachindex(O2)
        B[i] = O2[i]
    end
    B[1] += 2 * dt * D * O2sat / dz^2
    B[end] += 2 * dt * D * O2sat / dz^2
    return B
end

function diagonals(O2,D, dz, dt)
    coefDiff = D*dt/dz^2
    diag = fill(1 + 2 * coefDiff, length(O2))
    up = fill(-coefDiff, length(O2)-1)
    low = fill(-coefDiff, length(O2)-1)
    diag[1] = 1+3*coefDiff
    diag[end] = 1+3*coefDiff
    return low, diag, up
end

gp = GasesParams()
O2 = zeros(1000)
# O2[1:300] .= 0.5
# O2[600:700] .= 0.8
dz = 1e-6

t_tot = 600 #(s)
dt = 1e-4
n_iter::Int = round(Int, t_tot/dt)
println(n_iter)
n_save = 40
n_inner::Int = n_iter/n_save
println(n_inner)
S = zeros(1000)
fill!(S, 0.8)
for time_step in 1:n_save
    println("Start $time_step")
    for i_inner in 1:n_inner
        b = compB(O2, gp.O2sat, gp.D_oxygen,dz,dt)
        low, diag, up = diagonals(O2, gp.D_oxygen, dz, dt)
        LinearAlgebra.LAPACK.gtsv!(low, diag, up, b)
        O2 .= b
        if i_inner % 10000 == 0 #&& time_step == n_save
                # p = plot(scatter(x = eachindex(B), y = B, mode = "line" ),
                # Layout(title = "Source term"))
                # display(p)
                pp = plot(scatter(x = eachindex(O2), y = O2, mode = "line"), 
                Layout(title = "02 concentration profile $(i_inner/10000)",
                xaxis_title = "Depth (µm)",
                yaxis_title = "O2 concentration mol/m3"))
                display(pp)
            end
    end
end