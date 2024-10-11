using LinearAlgebra
using PlotlyJS
using Distributions
using Printf
include("model_parameters.jl")
include("functions_model.jl")

function grossmu!(µ_gross, I, k, sigma, tau, kd, kr, pop)
    ind = findfirst(x->x==0, pop)
    for i in 1:ind
        µ_gross[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2)
    end
end

function respiration!(R, I, RD, RL, Ik, n, pop)
    ind = findfirst(x->x==0, pop)
    for i in 1:ind
        R[i] = RD + (RL-RD)*(I[i]^n/(I[i]^n+Ik^n))
        # R[i] = 0.12/86400
    end
end


function O2source(mu, VO2x, mx, dz)
    source = zeros(ind)
    source[0:100] .= 0.3
    return source
end

function compB(O2,S, O2sat, D, dz, dt)
    B = zeros(length(O2))
    for i in eachindex(O2)
        B[i] = O2[i] + S[i] * dt
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


save_path = mkpath("C:/Users/Chrislain/Documents/Results/O2_nosource2")
gp = GasesParams()

# O2[600:700] .= 0.8
dz = 1e-6
dts = [1e-1,1e-2,1e-3,1e-4]
# divs = [10,100,1000,10000]
div = 1000
t_tot = 600 #(s)
deltaT = 1e-3
#a = 0.1, b = 0.01, c = 0.001, d = 0.0001
aa = []
bb = []
cc = []
dd = []

S = zeros(300)
# S[1:100] .= 0.03

println("Start")
# for (deltaT,div,arr) in zip(dts,divs,[aa,bb,cc,dd])
    O2 = zeros(300)
    # O2[1:300] .= 0.5
    # println("Start for $(deltaT)")
    # dt_good_format = @sprintf("%.e", deltaT)
    # save_path = mkpath("C:/Users/Chrislain/Documents/Results/O2_$(dt_good_format)s")
    n_iter::Int = round(Int, t_tot/deltaT)
    for i in 1:n_iter
        b = compB(O2,S, gp.O2sat, gp.D_oxygen,dz,deltaT)
        low, diag, up = diagonals(O2, gp.D_oxygen, dz, deltaT)
        LinearAlgebra.LAPACK.gtsv!(low, diag, up, b)
        O2 .= b
        if i % div == 0 #&& time_step == n_save
            time = Int(i/div)
            # p = plot(scatter(x = eachindex(B), y = B, mode = "line" ),
            # Layout(title = "Source term"))
            # display(p)
            
            pp = plot(scatter(x = eachindex(O2), y = O2, mode = "line"), 
            Layout(title = "02 concentration profile after $(time) s",
            xaxis_title = "Depth (µm)",
            yaxis_title = "O2 concentration mol/m3",
            yaxis_range = [-0.01, 0.55]))
            if time < 10
                file_name = "P_0$(time)_source.png"
            else
                file_name = "P_$(time)_source.png"
            end
            savefig(pp, joinpath(save_path, file_name))
            display(pp)
            # # #Test concat plot
            # push!(arr, scatter(x = eachindex(O2), y = O2, mode = "line"))
        end
    end
# end

# println("Start creating plot")
# for i in 1:600
#     traces = [a[i], b[i], c[i], d[i]]
#     pp = plot(traces, 
#     Layout(title = "02 concentration profile after $(i) s",
#     xaxis_title = "Depth (µm)",
#     yaxis_title = "O2 concentration mol/m3",
#     yaxis_range = [-0.01, 0.55]))
#     if i < 10
#         file_name = "P_0$i.png"
#     else
#         file_name = "P_$i.png"
#     end
#     savefig(pp, joinpath(save_path, file_name))
#     # display(pp)
# end

println("End")