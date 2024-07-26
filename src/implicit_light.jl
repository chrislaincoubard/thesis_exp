using PlotlyJS
using LinearAlgebra
using ColorSchemes
using DataFrames
using ColorTypes
using CSV
using Statistics
using DelimitedFiles
using LinearAlgebra.LAPACK
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
    A[1,1] = 1 +3coef
    A[1,2] = -coef
    A[ind, ind-1] = -coef
    A[ind, ind] = 1 + coef
    return A
end

function getdiagonals(DO2, dz, dt, pop)
    ind = findfirst(x -> x == 0, pop) -1
    diag = zeros(ind)
    up = zeros(ind-1)
    low = zeros(ind-1)
    coef = DO2*dt/dz^2
    diag .= 1+2*coef
    up .= -coef
    low .= -coef
    diag[1] = 1+3coef
    diag[ind] = 1+coef
    return low, diag, up
end

function ComputeSource(mu, VO2x, mx, pop, dz, dt)
    ind = findfirst(x -> x == 0, pop) -1 
    source = zeros(ind)
    for i in eachindex(source)
        source[i] = ((mu[i]*(pop[i]/dz)*VO2x)/mx)
    end
    return source
end

function computeB(O2, O2atm, D, dz, dt, pop,S)
    ind = findfirst(x -> x == 0, pop) -1 
    B = zeros(ind)
    for i in eachindex(B)[Not(1)]
        B[i] = O2[i] + S[i]*dt
        # B[i] = O2[i]
    end
    # B[1] = O2[1] + 2*dt*D*O2atm/dz^2
    B[1] = O2[1] + 2dt*D*O2atm/dz^2 +S[1]
    return B
end


tp = TimeParams()
hmp = HanModelParams() 
gp = GasesParams()
# light_intensities = [100,200,300,500,1000]
light_intensities = [200]
z = 1e-3
nz = 100000
dz = z/nz
zplt = 0:dz:z
time_save = zeros(tp.n_save)
df_mu = DataFrame()
df_height = DataFrame()
df = DataFrame()
df_O2 = DataFrame()
df_S = DataFrame()

for I0 in light_intensities
    ##Initialize Array
    LI = zeros(nz)
    µ = zeros(nz)
    µ_gross = zeros(nz)
    R = zeros(nz)
    pop = zeros(nz)
    height = zeros(tp.n_save)
    X0 = hmp.rho * dz
    O2 = zeros(nz)
    O22 = zeros(nz)
    pop[1:3000] .= X0
    println("Start for $I0")
    #Starting time (2 loops to save data at specific time points)
    for time_step in 1:tp.n_save
        println("Start $time_step")
        for i_inner in 1:tp.n_inner
            computelight!(LI, I0, hmp.ke, dz, pop)
            grossmu!(µ_gross, LI, hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, pop)
            respiration!(R, LI, hmp.RD, hmp.RL, hmp.Ik, hmp.n, pop)
            updatemu!(µ, µ_gross, R, pop)
            pop .= solvematrix(µ, tp.dt, nz, pop)
            ind = findfirst(x -> x == 0, pop)
            smootharray!(pop, X0)
            S = ComputeSource(µ, gp.VO2_x, gp.Mx,pop, dz, tp.dt)
            low, diag, up = getdiagonals(gp.D_oxygen, dz, tp.dt, pop)
            B = computeB(O2, gp.O2atm, gp.D_oxygen, dz, tp.dt, pop,S)
            LinearAlgebra.LAPACK.gtsv!(low, diag, up, B)
            O2[1:length(B)] .= B
            if i_inner == tp.n_inner #&& time_step == tp.n_save
                pp = plot(scatter(x = zplt*10^6, y = O2, mode = "line"), 
                Layout(title = "02 concentration profile",
                xaxis_title = "Depth (µm)",
                yaxis_title = "O2 concentration mol/m3"))
                display(pp)
                # pp = plot(scatter(x = zplt*10^6, y = S, mode = "markers"), 
                # Layout(title = "O2 production by photosynthesis",
                # xaxis_title = "Depth (µm)",
                # yaxis_title = "O2 production mol/m3/s"))
                # display(pp)
                # savefig(p, "C:/Users/chris/OneDrive/Images/plots/O2 prod.png")
            end
        end
        time = tp.dt*time_step*tp.n_inner/3600
        currheight = sum(pop) / hmp.rho
        time_save[time_step] = time
        height[time_step] = currheight
        # pl = plot(scatter(x = zplt, y = S, mode = "markers"), Layout(title = "test_source"))
        # display(pl)
        # df_O2[!,"$time_step"] = O2
    end
    colname = "$I0"
    df[!,"µ_net"] = µ
    df[!,"gross_µ"] = µ_gross
    df[!,"R"] = R
    df[!,"light"] = LI
    df[!,"Intensity"] .= I0
    df_height[!,"Height"] = height
    df_height[!,"Intensity"] .= I0
    # df[!,"oxygen"] = O2
    # df[!,colname] = height
    # df_O2[!,colname] = O2
end
println("Done with calculation")

##################################################################
#                                                                #
#                            Plots                               #
#                                                                #
##################################################################
function cleanarr(arr)
    return filter(x->x!=0, arr)    
end

clean_R = cleanarr(df[!,"R"])
clean_µ = cleanarr(df[!,"µ_net"])
clean_grossµ = cleanarr(df[!,"gross_µ"])
clean_light = cleanarr(df[!,"light"])
# clean_oxygen = cleanarr(df[!,"oxygen"])

pltR = plot([
    scatter(x = zplt*10^6, y = clean_R.*86400, mode = "line", name = "Respiration"),
    scatter(x = zplt*10^6, y = clean_grossµ.*86400, mode = "line", name = "µ brut"),
    scatter(x = zplt*10^6, y = clean_µ.*86400, mode = "line", name = "µ net")],
    Layout(title = "growth rate", xaxis_title = "Depth of biofilm (µm)", yaxis_title = "s-1"))


pltLight = plot(scatter(x = zplt*10^6, y = clean_light, mode = "markers"), 
Layout(title = "Light", xaxis_title = "Depth (µm)", yaxis_title = "PPFD (µmol/m²/s)"))
pltheight = plot(scatter(x = time_save, y = df_height[!,"Height"].*10^6, mode = "line"), 
    Layout(title = "Biofilm Growth", xaxis_title = "Time (h)", yaxis_title = "Height (µm)"))
display(pltheight)
# pltOxygen = plot(scatter(x = zplt*10^6, y = clean_oxygen, mode = "markers"), 
# Layout(title = "Oxygen", xaxis_title = "Depth (µm)", yaxis_title = "O2 concentration in mM"))

display(pltR)
# savefig(pltheight, "C:/Users/chris/OneDrive/Images/plots/height_200_epsi2e4.png")
# display(pltLight)
# display(pltOxygen)
# println(first(df_mu,5))

# function buildO2traces(df, x)
#     traces = Vector{GenericTrace}(undef, length(names(df)))
#     for (index, n) in enumerate(names(df))
#     dfcol_clean = filter(x -> x!= 0, df[!,n])
#     sct = scatter(y = dfcol_clean, x = x, mode = "markers",
#      marker = attr(size = 3), name = n)
#     traces[index] = sct
#     end
#     return traces
# end

# mu_traces = buildmutraces(df_mu, zplt)

# O2_traces = buildO2traces(df_O2, zplt)

# m=plot(O2_traces, Layout(title = "test"))

# display(m)

# plotmutraces(mu_traces, "C:/Users/LGPM Bamako/Documents/Results", "mu_plot.png")
# # writedlm( "height.csv",  df_height, ',')
# println(first(df,5))

# pl = plot([
#     # scatter(x = time_save, y = df."100", mode = "markers", name = "100"), 
#     # scatter(x = time_save, y = df."200", mode = "markers", name = "200"),
#     scatter(x = time_save, y = df."300", mode = "markers", name = "300"),
#     # scatter(x = time_save, y = df."500", mode = "markers", name = "500"),
#     # scatter(x = time_save, y = df."1000", mode = "markers", name = "1000"),
# ], Layout(title = "Biofilm Height"))
# display(pl)

# pla = plot([
#     scatter(x = zplt, y = df_mu."100".*86400, mode = "markers", name = "100"), 
#     scatter(x = zplt, y = df_mu."200".*86400, mode = "markers", name = "200"),
#     scatter(x = zplt, y = df_mu."300".*86400, mode = "markers", name = "300"),
#     scatter(x = zplt, y = df_mu."500".*86400, mode = "markers", name = "500"),
#     scatter(x = zplt, y = df_mu."1000".*86400, mode = "markers", name = "1000"),
# ])
# display(pla)
# plpl = plot(O2_traces)
# display(plpl)
# println(first(df_O2, 50))

# plot(surface(
#     z=matO2,x=time_save,y=zplt, colorscale = "Earth"))
