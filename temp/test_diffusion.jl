using LinearAlgebra
using PlotlyJS
using Distributions
using Printf
using DataFrames
include("model_parameters.jl")
include("functions_model.jl")

function onegrossmu!(I, k, sigma, tau, kd, kr)
    return (k*sigma*I) / (1 + tau*sigma*I + (kd/kr)*tau*(sigma*I)^2)
    
end

function onerespiration!(I, RD, RL, Ik, n)
    return RD + (RL-RD)*(I^n/(I^n+Ik^n))
end


function oneO2source(mu,X, VO2x, mx, dz)
    return ((mu*(X/dz)*VO2x)/mx)
end

function compB(O2,S, O2sat, D, dz, dt)
    B = zeros(length(O2))
    for i in eachindex(O2)
        B[i] = O2[i] + S * dt
        # B[i] = O2[i]
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


save_path = mkpath("C:/Users/Chrislain/Documents/Results/O2_source2")
gp = GasesParams()
hmp = HanModelParams()



# O2[600:700] .= 0.8
dz = 1e-6
dts = [1e-1,1e-2,1e-3,1e-4]
X = 3e-4*1.4e2
Xdz = X/300
# divs = [10,100,1000,10000]

t_tot = 600 #(s)
deltaT = 1e-2
div = 1/deltaT
#a = 0.1, b = 0.01, c = 0.001, d = 0.0001
aa = []
bb = []
cc = []
dd = []
Light = 200
S = zeros(300)
# S[1:100] .= 0.03

mu_g = onegrossmu!(Light, hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr)
println(mu_g)
r = onerespiration!(Light, hmp.RD, hmp.RL, hmp.Ik,hmp.n)
println(r)

mu_net = mu_g - r
println(mu_net)

SO2 = oneO2source(mu_net,Xdz, gp.VO2_x, gp.Mx,dz)
println(SO2)
df = DataFrame()
println("Start")
# for (deltaT,div,arr) in zip(dts,divs,[aa,bb,cc,dd])
    O2 = zeros(300)
    # O2[1:300] .= 0.5
    # println("Start for $(deltaT)")
    # dt_good_format = @sprintf("%.e", deltaT)
    # save_path = mkpath("C:/Users/Chrislain/Documents/Results/O2_$(dt_good_format)s")
    n_iter::Int = round(Int, t_tot/deltaT)
    for i in 1:n_iter
        b = compB(O2,SO2, gp.O2sat, gp.D_oxygen,dz,deltaT)
        low, diag, up = diagonals(O2, gp.D_oxygen, dz, deltaT)
        LinearAlgebra.LAPACK.gtsv!(low, diag, up, b)
        O2 .= b
        if i % div == 0 #&& time_step == n_save
            time = Int(i/div)
            # p = plot(scatter(x = eachindex(B), y = B, mode = "line" ),
            # Layout(title = "Source term"))
            # display(p)
            df[!,"$(time)"] = O2
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
            # save_path = mkpath("C:/Users/Chrislain/Documents/Results/O2_sourcess")
            # savefig(pp, joinpath(save_path, file_name))
            # display(pp)
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

xtest = collect(0:1e-6:3e-4)

function analytical_source(x, S, L, cO2,DO2)
    a = -S/(2*DO2)
    b = -a*L
    @. return a*x^2 + b*x + cO2   
end

y_nosource = fill(gp.O2sat,300)


L = 3e-4
y = analytical_source(xtest, SO2, L, gp.O2sat,gp.D_oxygen)

p_ss_nosource = plot([scatter(x = xtest.*1e6, y = y_nosource, mode = "line", name = "analytical"),
            scatter(x = eachindex(df."600"), y = df."600", mode = "line", name = "FVM")],
            Layout(title = "Steady state no source",xaxis_title = "Depth (µm)",
            yaxis_title = "O2 concentration mol/m3"))
display(p_ss_nosource)
p_ss_source = plot([scatter(x = xtest.*1e6, y = y, mode = "line", name = "analytical"),
                            scatter(x = eachindex(df."600"), y = df."600", mode = "line", name = "FVM")],
                            Layout(title = "Steady State source"))
display(p_ss_source)
save_gpath = mkpath("C:/Users/Chrislain/Documents/Results/comparaison_FVM_Analytical")

# savefig(p_ss_nosource, joinpath(save_gpath,"FVMvsAnalytical_steady_state_nosource_zoom.png"))
savefig(p_ss_source, joinpath(save_gpath, "FVMvsAnalytical_steady_state_source.png"))
println("End")
