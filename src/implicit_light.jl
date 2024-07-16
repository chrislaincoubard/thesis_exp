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


function makeO2coefmatrix(DO2, dz, dt, pop)
    ind = findfirst(x -> x == 0, pop) -1
    A = zeros(ind, ind)
    coef = DO2*dt/dz^2
    for i in 2:ind-1
        A[i, i-1] = -coef
        A[i, i+1] = -coef
        A[i, i] = 1+ 2*coef
    end
    A[1,1] = 1 +3*coef
    A[1,2] = -coef
    A[ind, ind-1] = -coef
    A[ind, ind] = 1 + coef
    return A
end

function make02Bmatrix(O2, O2atm, D, dz, dt, pop,mu)
    ind = findfirst(x -> x == 0, pop) -1 
    B = zeros(ind)
    for i in eachindex(B)
        B[i] = O2[i] +(mu[i]*gp.VO2_x)/(pop[i]*gp.Mx)*dt
    end
    B[1] = O2[1] + 2*dt*D*O2atm/dz^2
    return B
end


tp = TimeParams()
hmp = HanModelParams() 
gp = GasesParams()
light_intensities = [100,200,300,500,1000]
# light_intensities = [200]
z = 1e-3
nz = 1000
dz = z/nz
zplt = 0:dz:z
matO2 = zeros(tp.n_save, nz)
time_save = zeros(tp.n_save)
df_mu = DataFrame()
df_height = DataFrame(Height = Float64[], Intensity=String[], Time = Float64[])
df = DataFrame()
df_O2 = DataFrame()

for I0 in light_intensities
    LI = zeros(nz)
    µ = zeros(nz)
    pop = zeros(nz)
    height = zeros(tp.n_save)
    X0 = hmp.rho * dz
    O2 = zeros(nz)
    pop[1:30] .= X0
    println("Start for $I0")
    for time_step in 1:tp.n_save
        for i_inner in 1:tp.n_inner
            computelight!(LI, I0, hmp.ke, dz, pop)
            updatemu!(µ, hmp.RD, hmp.RL, LI, hmp.n, hmp.Ik,hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, pop)
            pop .= solvematrix(µ, tp.dt, nz, pop)
            A = makeO2coefmatrix(gp.D_oxygen, dz, tp.dt, pop)
            B= make02Bmatrix(O2, gp.O2atm, gp.D_oxygen,dz, tp.dt, pop,µ)
            if i_inner in 1:2 && time_step == 1
                println(B)
            end
            new_O2 = A \ B
            O2[1:length(new_O2)] .= new_O2
            smootharray!(pop, nz, X0)
        end
        time = tp.dt*time_step*tp.n_inner/3600
        currheight = sum(pop) / hmp.rho
        time_save[time_step] = time
        height[time_step] = currheight
        push!(df_height, [currheight, "$I0", time])
        # df_O2[!,"$time_step"] = O2
        matO2[time_step,:] .= O2 
    end
    
    colname = "$I0"
    df_mu[!,colname] = µ
    df[!,colname] = height
    df_O2[!,colname] = O2
end
println("Done with calulation")

######## Plots #############
# println(first(df_mu,5))

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

# mu_traces = buildmutraces(df_mu, zplt)

O2_traces = buildO2traces(df_O2, zplt)

plot(O2_traces, Layout(title = "test"))

# plotmutraces(mu_traces, "C:/Users/LGPM Bamako/Documents/Results", "mu_plot.png")
# # writedlm( "height.csv",  df_height, ',')
# println(first(df,5))

# pl = plot([
#     scatter(x = time_save, y = df."100", mode = "markers", name = "100"), 
#     scatter(x = time_save, y = df."200", mode = "markers", name = "200"),
#     scatter(x = time_save, y = df."300", mode = "markers", name = "300"),
#     scatter(x = time_save, y = df."500", mode = "markers", name = "500"),
#     scatter(x = time_save, y = df."1000", mode = "markers", name = "1000"),
# ])
# display(pl)
# plpl = plot(O2_traces)
# display(plpl)
# println(first(df_O2, 50))

# plot(surface(
#     z=matO2,x=time_save,y=zplt, colorscale = "Earth"))

