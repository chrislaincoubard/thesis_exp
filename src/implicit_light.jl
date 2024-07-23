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

function computesource(mu, VO2x, mx, pop, dz,dt)
    ind = findfirst(x -> x == 0, pop) - 1
    source = zeros(ind)
    for i in eachindex(source)
        source[i] = ((mu[i]*pop[i])*VO2x*dt)/(mx*dz)
    end
    return source
end

function make02Bmatrix(O2, O2surf, D, dz, dt, pop, S)
    ind = findfirst(x -> x == 0, pop) -1 
    B = zeros(ind)
    for i in eachindex(B)[Not(1)]
        B[i] = O2[i] +S[i]
    end
    B[1] = O2[1] + (2*dt*D*O2surf)/dz^2 + S[1]
    return B
end


tp = TimeParams()
hmp = HanModelParams() 
gp = GasesParams()
# light_intensities = [100,200,300,500,1000]
light_intensities = [500]
z = 1e-3
nz = 1000
dz = z/nz
zplt = 0:dz:z
matO2 = zeros(tp.n_save, nz)
time_save = zeros(tp.n_save)
df_H = DataFrame()
df_height = DataFrame(Height = Float64[], Intensity=String[], Time = Float64[])
df = DataFrame()

for I0 in light_intensities
    LI = zeros(nz)
    µ_gross = zeros(nz)
    R = zeros(nz)
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
            grossmu(µ_gross, LI, hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, pop)
            respiration(R, LI, hmp.RD, hmp.RL, hmp.Ik, hmp.n, pop)
            updatemu!(µ, hmp.RD, hmp.RL, LI, hmp.n, hmp.Ik,hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, pop)
            pop .= solvematrix(µ, tp.dt, nz, pop)
            smootharray!(pop, X0)
            A = makeO2coefmatrix(gp.D_oxygen, dz, tp.dt, pop)
            S = computesource(µ, gp.VO2_x, gp.Mx, pop, dz, tp.dt)
            B= make02Bmatrix(O2, gp.O2surf, gp.D_oxygen,dz, tp.dt, pop,S)
            if i_inner == tp.n_inner && time_step == tp.n_save
                clean_source = filter(x-> x!=0, S)
                pip = plot(scatter(x = zplt, y = clean_source, mode = "markers"), Layout(title = "Source 02"))
                display(pip)
            end
            new_O2 = A \ B
            O2[1:length(new_O2)] .= new_O2
        end
        time = tp.dt*time_step*tp.n_inner/3600
        currheight = sum(pop) / hmp.rho
        time_save[time_step] = time
        height[time_step] = currheight
        push!(df_height, [currheight, "$I0", time])
        # df_O2[!,"$time_step"] = O2
        # matO2[time_step,:] .= O2 
    end
    colname = "$I0"
    df[!,"µ_net"] = µ
    df[!,"µ_gross"] = µ_gross
    df[!,"R"] = R
    df_H[!,"height"] = height
    df[!,"light"] = LI
    df[!,"Intensity"] .= I0
    df[!,"Oxygen"] = O2
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
clean_O2 = filter(x -> x!=0, df[!,"Oxygen"])
plott = plot(scatter(x = zplt, y = clean_O2, mode = "markers"),Layout(title = "O2 along depth"))
display(plott)

# plotmutraces(mu_traces, "C:/Users/LGPM Bamako/Documents/Results", "mu_plot.png")
# # writedlm( "height.csv",  df_height, ',')
# println(first(df,5))
name = "$(light_intensities[1])"

clean_R = filter(x -> x!= 0, df[!,"R"])
clean_µ = filter(x -> x!= 0, df[!,"µ_net"])
clean_grossµ = filter(x -> x!= 0, df[!,"µ_gross"])
clean_light = filter(x -> x!=0, df[!,"light"])



pl = plot([
    # scatter(x = zplt, y = df_mu."100".*86400, mode = "markers", name = "100"), 
    scatter(x = zplt, y = clean_grossµ*86400, mode = "markers", name = "µ gross"),
    scatter(x = zplt, y = clean_µ.*86400, mode = "markers", name = "µ net"),
    scatter(x = zplt, y = clean_R*86400, mode = "markers", name = "R")
    
    # scatter(x = zplt, y = df_mu."300".*86400, mode = "markers", name = "300"),
    # scatter(x = zplt, y = df_mu."500".*86400, mode = "markers", name = "500"),
    # scatter(x = zplt, y = df_mu."1000".*86400, mode = "markers", name = "1000"),
], Layout(title = "growth rate"))
display(pl)

pp = plot(scatter(x = zplt, y = clean_light, mode = "markers"), Layout(title = "Light"))
display(pp)

# plpl = plot(O2_traces)
# display(plpl)
# println(first(df_O2, 50))

# plot(surface(
#     z=matO2,x=time_save,y=zplt, colorscale = "Earth"))

