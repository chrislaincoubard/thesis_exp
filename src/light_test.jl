using PlotlyJS
include("parameters_light.jl")
using Statistics
using DelimitedFiles

function computelight!(arrLight, I0, ke, dx, nx,pop)
    arrLight[1] = I0
    for i in 1:nx-1
        if pop[i] != 0
            arrLight[i+1] = arrLight[i]*((1-(ke*dx*i)))
        end 
    end
end

function computerespiration(RD, RL, I, n, Ik)
    return RD + (RL - RD) * (I^n/(I^n + Ik^n))
end

function updatemu!(µ, RD, RL, I, n, Ik, k, sigma, tau, kd, kr, nx, pop)
    for i in 1:nx
        if pop[i] != 0
            R = RD + (RL - RD) * (I[i]^n/(I[i]^n + Ik^n))
            µ[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2) - R
        end
    end
end

function updateheight!(pop,µ,dt)
    for i in eachindex(pop)
        pop[i] = pop[i] * µ[i] * dt + pop[i]
    end
end

function smootharray!(pop, nx, dx)
    for i in 1:nx-1
        if pop[i] > dx
            pop[i+1] = pop[i+1] + pop[i] -dx
            pop[i] = dx
        end
    end
end



tp = TimeParams()
mp = ModelParams()

x = 1e-4 
nx = 1000
dx = x/nx
xplt = 0:dx:x

DLI = zeros(nx)
light_save = zeros(tp.n_save, nx)
µ = zeros(nx)

height = zeros(nx)
height[1:30] .= dx
mumean_save = zeros(tp.n_save)
µ_save = zeros(tp.n_save, nx)
time_save = zeros(tp.n_save)
pop_save = zeros(tp.n_save, nx)
height_save = zeros(tp.n_save)
height_mean = zeros(tp.n_save)
height_mean[1] = 3e-6
#Random testing

for time_step in 1:tp.n_save
    for i_inner in 1:tp.n_inner
        computelight!(DLI, mp.I0, mp.ke, dx, nx, height)
        updatemu!(µ, mp.RD, mp.RL, DLI, mp.n, mp.Ik,mp.k, mp.sigma, mp.tau, mp.kd, mp.kr, nx, height)
        updateheight!(height, µ, tp.dt)
        smootharray!(height, nx, dx)

    end
    if time_step != 1
        mumean = mean(µ)
        height_mean[time_step] = mumean * height_mean[time_step-1] * tp.n_inner + height_mean[time_step-1]
    end
    time_save[time_step] = tp.dt*time_step*tp.n_inner/3600
    height_save[time_step] = sum(height)
    µ_save[time_step,:] = µ
    light_save[time_step,:] = DLI
    pop_save[time_step,:] = height
end


########## Plots ##############


layout1 = Layout(
    title = "Light attenuation",
    xaxis_title = "depth (m)",
    yaxis_title = "light intensity (µmol*m-2*s-1)"
)

layout2 = Layout(
    title = "growth rate variation over height of biofilm", 
    xaxis_title = "depth (m)", 
    yaxis_title = "growth rate (d-1)"
)

layout3 = Layout(
    title = "Height of biofilm",
    xaxis_title = "time (h)",
    yaxis_title = "height (m)"
)


traces = Vector{GenericTrace}(undef, tp.n_save)
for i in 1:tp.n_save
    t = i*tp.dt*tp.n_inner/3600
    traces[i] = scatter(x = xplt, y = light_save[i,:], 
        mode = "line", name = "$t h",
        line = attr(color = "red", width = 1))
end

traces_µ = Vector{GenericTrace}(undef, tp.n_save)
for i in 1:tp.n_save
    t = i*tp.dt*tp.n_inner/3600
    traces_µ[i] = scatter(x = xplt, y = µ_save[i,:], 
        mode = "line", name = "$t h",
        line = attr(color = "red", width = 1))
end

plt = plot(
    scatter(y = DLI, x = xplt*1e6, mode = "line", line = attr(color = "blue")), layout1)
display(plt)

plt2 = plot(scatter(y = µ, x = xplt, mode = "markers"), layout2)
display(plt2)

plt3 = plot(scatter(x = time_save, y = height_save, mode = "markers"), layout3)
display(plt3)

plt_light = plot(traces)
display(plt_light)

plt_mu = plot(traces_µ)
display(plt_mu)
# print(height[1:30])
muuu = filter(x -> x!=0, µ)
println(mean(muuu)*86400)
# writedlm("mu_save.csv", µ_save, ",")
# writedlm("height_save.csv", height_save, ",")
# writedlm("pop_save.csv",pop_save,",")