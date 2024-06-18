using PlotlyJS
include("parameters_light.jl")
using Statistics

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

function computegrowthrate!(µ,i, I, R, k, sigma, tau, kd, kr)
    µ[i] = (k*sigma*I) / (1 + tau*sigma*I + (kd/kr)*tau*(sigma*I)^2) - R 
end

function updatemu!(µ, RD, RL, I, n, Ik, k, sigma, tau, kd, kr, nx, pop)
    for i in 1:nx
        if pop[i] != 0
            R = RD + (RL - RD) * (I[i]^n/(I[i]^n + Ik^n))
            µ[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2) - R
        end
    end
end

function updateheight!(pop,µ)
    for i in eachindex(pop)
        pop[i] = pop[i] * µ[i] + pop[i]
    end
end

function smootharray!(pop, nx)
    for i in 1:nx-1
        if pop[i] > dx
            pop[i+1] = pop[i+1] + pop[i] -dx
            pop[i] = dx
        end
    end
end

# sp = SpaceParams()
tp = TimeParams()
mp = ModelParams()

x = 5e-4 
nx = 1000
dx = x/nx
xplt = 0:dx:x


space = zeros(nx)
space[1:100] .= mp.ke
DLI = zeros(nx)
µ = zeros(nx)
height = zeros(nx)
height[1:100] .= dx

µ_save = zeros(tp.n_save, nx)
time_save = zeros(tp.n_save)
height_save = zeros(tp.n_save)
#Random testing
dt = 150
nt = tp.t_tot / dt
height2 = zeros(nx)

for time_step in 1:tp.n_save
    for i_inner in 1:tp.n_inner
        computelight!(DLI, mp.I0, mp.ke, dx, nx, height)
        updatemu!(µ, mp.RD, mp.RL, DLI, mp.n, mp.Ik,mp.k, mp.sigma, mp.tau, mp.kd, mp.kr, nx, height)
        updateheight!(height, µ)
        smootharray!(height, nx)

    end
    if time_step == 2
        for i in 1:nx
            if height[i] != 0
                height[i] != 0
                println("height[i] : ", height[i], " -- i : ",i)
            end
        end
    end
    time_save[time_step] = tp.dt*time_step*tp.n_inner/3600
    height_save[time_step] = sum(height)
    µ_save[time_step,:] = µ
end


########## Plots ##############


layout1 = Layout(
    title = "Light attenuation",
    xaxis_title = "Height (mm)",
    yaxis_title = "light intensity (µmol*m-2*s-1)"
)

layout2 = Layout(
    title = "growth rate variation over height of biofilm", 
    xaxis_title = "height (mm)", 
    yaxis_title = "growth rate (d-1)"
)

# layout3 = Layout(
#     title = "Height of biofilm",
#     xaxis_title = "time (h)",
#     yaxis_title = "height (m)"
# )
plt = plot(
    scatter(y = DLI, x = xplt*1e3, mode = "line", line = attr(color = "blue")), layout1)
display(plt)

plt2 = plot(scatter(y = µ*86400, x = xplt*1e4, mode = "markers"), layout2)
display(plt2)

plt3 = plot(scatter(x = 0:40, y = height_save, mode = "markers"))
display(plt3)
# ppp = plot(scatter(y=height, x = 0:ttest, mode = "line"), layout3)
# display(ppp)