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
height_save = zeros(tp.n_save)

df_mu = DataFrame()
df_light = DataFrame()
df_height = DataFrame(Height = Float64[], Intensity=Int64[], Time = Int64[])
print(df_height)


for (i,I0) in enumerate(light_intensities)
    LI = zeros(nz)
    µ = zeros(nz)
    pop = zeros(nz)
    height = zeros(tp.n_save)
    X0 = mp.rho * dz
    pop[1:30] .= X0
    println("Start for $I0")
    for time_step in 1:tp.n_save
        for i_inner in 1:tp.n_inner
            computelight!(LI, I0, mp.ke, dz, nz, pop)
            updatemu!(µ, mp.RD, mp.RL, LI, mp.n, mp.Ik,mp.k, mp.sigma, mp.tau, mp.kd, mp.kr, nz, pop)
            pop .= solvematrix(µ, tp.dt, nz, pop)
            smootharray!(pop, nz, X0)
        end
        time = tp.dt*time_step*tp.n_inner/3600
        currheight = sum(pop) / mp.rho
        time_save[time_step] = time
        height_save[time_step] = currheight
        height[time_step] = currheight
        µ_save[time_step,:] = µ
        pop_save[time_step,:] = pop
        push!(df_height, [currheight, I0, time])
    end
    
    colname = "$I0"
    df_mu[!,colname] = µ
end

layout1 = Layout(
    title = "Light attenuation",
    xaxis_title = "depth ",
    yaxis_title = "light intensity (µmol*m-2*s-1)"
)



layout3 = Layout(
    title = "pop of biofilm",
    xaxis_title = "time (h)",
    yaxis_title = "pop (m)"
)

# plt = plot(
#     scatter(y = LI, x = zplt, mode = "line", line = attr(color = "blue")), layout1)
# display(plt)


layout2 = Layout(
    title = "growth rate variation over depth of biofilm", 
    xaxis_title = "depth (m)", 
    yaxis_title = "growth rate (d-1)"
    )


# plt2 = plot(scatter(y = df_mu."200".*86400, x = zplt, mode = "markers"), layout2)
# display(plt2)

plt3 = plot(scatter(x = time_save, y = height_save, mode = "markers"), layout3)
display(plt3)

col = reverse(ColorSchemes.rainbow)
# bar_traces = Vector{GenericTrace}(undef, tp.n_save)




# plt_histo = plot(traces_bar, Layout(barmode = "stack"))
# display(plt_histo)
# print(first(df_mu,5))
# cleanhisto = mapcols(col -> barheight(height_save),df_height)
# df_height[!,"Time"] = time_save

# plot(cleanhisto,kind = "bar", x=:Time, y =:100, Layout(barmode = "stack"))

function barheight(height_arr)
    barH = zeros(length(height_arr))
    barH[1] = height_arr[1]
    for i in eachindex(height_arr)[Not(1)]
        barH[i] = height_arr[i] - height_arr[i-1]
    end
    return barH
end

function bar_df!(df)
    for n in names(df)
        df[!,n] = barheight(df[!,n])
    end
end

function trace1bar!(traces, data, time, I, colors)
    for i in eachindex(time)
        t = time[i]
        scale_i = (i-1)/(length(time) - 1)
        traces[i] = bar(name = "$t h", x = ["$I"], y = [data[i]],
        width = [0.2], marker_color = colors[scale_i])
    end
    return traces_bar
end
# bar_df!(df_height)

print(first(df_height,5))



# traces_bar = Vector{GenericTrace}(undef, tp.n_save)
# tb = Vector{GenericTrace}(undef, tp.n_save)
# barh = barheight(df_height."100")
# barh2 = barheight(df_height."200")
# trrr = trace1bar!(traces_bar, barh, time_save, 100, col)
# trrr2 = trace1bar!(traces_bar, barh2, time_save, 200, col)
# pp =plot(trrr, Layout(barmode = "stack", title = "this one"))
# display(pp)