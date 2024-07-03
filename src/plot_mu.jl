using PlotlyJS
using Statistics
using LinearAlgebra
using ColorSchemes
using DataFrames
include("model_parameters.jl")

function computelight!(arrLight, I0, ke, dz, nz)
    arrLight[1] = I0
    for i in 2:nz
        arrLight[i] = I0 * exp(-ke*(dz*i))
    end
    
end

function updatemu!(µ, RD, RL, I, n, Ik, k, sigma, tau, kd, kr, nz)
    for i in 1:nz
        R = RD + (RL - RD) * (I[i]^n/(I[i]^n + Ik^n))
        µ[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2) - R
    end
end

function buildmutraces(df, x)
    traces = Vector{GenericTrace}(undef, length(names(df)))
    for (index, n) in enumerate(names(df))
    dfcol_clean = filter(x -> x!= 0, df[!,n])
    sct = scatter(y = dfcol_clean.*86400, x = x, mode = "markers",
     marker = attr(size = 3), name = n)
    traces[index] = sct
    end
    return traces
end


mp = ModelParams() 
light_intensities = [100,200,300,500,1000]
z = 2e-4
nz = 10000
dz = z/nz
zplt = 0:dz:z
df_mu = DataFrame()


for I0 in light_intensities
    LI = zeros(nz)
    µ = zeros(nz)
    computelight!(LI, I0, mp.ke, dz, nz)
    updatemu!(µ, mp.RD, mp.RL, LI, mp.n, mp.Ik,mp.k, mp.sigma, mp.tau, mp.kd, mp.kr, nz)
    colname = "$I0"
    df_mu[!,colname] = µ
end

layout = Layout(
    title = "growth rate variation over depth of biofilm", 
    xaxis_title = "depth (m)", 
    yaxis_title = "growth rate (d-1)",
    legend_title_text="Incident light intensity",
    legend_itemsizing="constant"
    )

# traces = buildmutraces(df_mu, zplt)

# plttest = plot(traces, layout)
# display(plttest)

# savefig(plttest, "mu.png")

println(mean(df."200"))