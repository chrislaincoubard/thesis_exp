using PlotlyJS
using Statistics
using LinearAlgebra
using ColorSchemes
using DataFrames
include("parameters_light.jl")

function computelight!(arrLight, I0, ke, dz, nz,pop)
    arrLight[1] = I0
    for i in 2:nz
        if pop[i] != 0
        arrLight[i] = I0 * exp(-ke*(dz*i))
        end 
    end
    
end

function updatemu!(µ, RD, RL, I, n, Ik, k, sigma, tau, kd, kr, nz, pop)
    for i in 1:nz
        if pop[i] != 0
            R = RD + (RL - RD) * (I[i]^n/(I[i]^n + Ik^n))
            µ[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2) - R
        end
    end
end

function solvematrix(mu, deltaT, nz, b)
    diag = zeros(nz)
    for i in 1:nz
        diag[i] = 1- mu[i] * deltaT
    end
    mat_diag = Diagonal(diag)
    return mat_diag \ b
end

function smootharray!(pop, nz, X)
    for i in 1:nz-1
        if pop[i] > X
            pop[i+1] = pop[i+1] + pop[i] -X
            pop[i] = X
        end
    end
end

function barheight(height_arr)
    barH = zeros(length(height_arr))
    barH[1] = height_arr[1]
    for i in eachindex(height_arr)[Not(1)]
        barH[i] = height_arr[i] - height_arr[i-1]
    end
    return barH
end

tp = TimeParams()
mp = ModelParams() 
light_intensities = [100,200,300,500,1000]
z = 1e-3
nz = 10000
dz = z/nz
zplt = 0:dz:z

µ_save = zeros(tp.n_save, nz)
time_save = zeros(tp.n_save)
pop_save = zeros(tp.n_save, nz)

df_mu = DataFrame()
df_light = DataFrame()
df_height = DataFrame()


for (i,I0) in enumerate(light_intensities)
    LI = zeros(nz)
    µ = zeros(nz)
    pop = zeros(nz)
    X0 = mp.rho * dz
    pop[1:30] .= X0
    println("Start for $I0")
    height_save = zeros(tp.n_save)
    for time_step in 1:tp.n_save
        for i_inner in 1:tp.n_inner
            computelight!(LI, I0, mp.ke, dz, nz, pop)
            updatemu!(µ, mp.RD, mp.RL, LI, mp.n, mp.Ik,mp.k, mp.sigma, mp.tau, mp.kd, mp.kr, nz, pop)
            pop .= solvematrix(µ, tp.dt, nz, pop)
            smootharray!(pop, nz, X0)
        end
        time_save[time_step] = tp.dt*time_step*tp.n_inner/3600
        height_save[time_step] = sum(pop) / mp.rho
        µ_save[time_step,:] = µ
        pop_save[time_step,:] = pop
    end
    colname = "$I0"
    df_mu[!,colname] = µ
    df_height[!,colname] = barheight(height_save)
end

layout1 = Layout(
    title = "Light attenuation",
    xaxis_title = "depth ",
    yaxis_title = "light intensity (µmol*m-2*s-1)"
)

layout2 = Layout(
    title = "growth rate variation over depth of biofilm", 
    xaxis_title = "depth (m)", 
    yaxis_title = "growth rate (d-1)"
)

layout3 = Layout(
    title = "pop of biofilm",
    xaxis_title = "time (h)",
    yaxis_title = "pop (m)"
)

# plt = plot(
#     scatter(y = LI, x = zplt, mode = "line", line = attr(color = "blue")), layout1)
# display(plt)

# plt2 = plot(scatter(y = µ.*86400, x = zplt, mode = "markers"), layout2)
# display(plt2)

plt3 = plot(scatter(x = time_save, y = height_save, mode = "markers"), layout3)
display(plt3)

col = reverse(ColorSchemes.rainbow)
bar_traces = Vector{GenericTrace}(undef, tp.n_save)


barH = barheight(height_save)

for i in 1:tp.n_save
    t = time_save[i]
    scale_i = (i-1)/(tp.n_save - 1)
    bar_traces[i] = bar(name = "$t h", x = ["$(mp.I0)"], y = [barH[i]],
    width = [0.2], marker_color = col[scale_i])
end


plt_histo = plot(bar_traces, Layout(barmode = "stack"))
display(plt_histo)
print(first(df_mu,5))
cleanhisto = mapcols(col -> barheight(height_save),df_height)
df_height[!,"Time"] = time_save
println(df_height)
plot(cleanhisto,kind = "bar", x=:Time, y =:100, Layout(barmode = "stack"))