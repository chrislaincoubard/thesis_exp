using PlotlyJS
include("parameters_light.jl")
using Statistics

function computelight!(arrLight, I0, ke, dx, nx)
    arrLight[1] = I0
    for i in 1:nx-1
        arrLight[i+1] = arrLight[i]*((1-(ke*dx*i)))
    end
end

function computerespiration(RD, RL, I, n, Ik)
    return RD + (RL - RD) * (I^n/(I^n + Ik^n))
end

function computegrowthrate!(µ,i, I, R, k, sigma, tau, kd, kr)
    µ[i] = (k*sigma*I) / (1 + tau*sigma*I + (kd/kr)*tau*(sigma*I)^2) - R 
end

# sp = SpaceParams()
tp = TimeParams()
mp = ModelParams()


x = 1e-4 
nx = 1000
dx = x/nx

DLI = zeros(nx)
xplt = 0:dx:x
µ = zeros(nx)
pop = zeros(nx)



# DLI[1] = mp.I0
for i in 1:nx
    R = computerespiration(mp.RD, mp.RL, DLI[i], mp.n, mp.Ik)
    mu = computegrowthrate!(µ, i, DLI[i], R, mp.k, mp.sigma, mp.tau, mp.kd, mp.kr)
end


# for i in 1:nx-1
#     DLI[i+1] = DLI[i]*((1-(mp.ke*dx*i)))
# end

computelight!(DLI, mp.I0, mp.ke, dx, nx)

#
mumean = mean(µ)*3600

#Random testing
x0 = 1e-6
ttest = 180
biom = zeros(ttest)
biom[1] = x0
cb = 1.4e5
height = zeros(ttest)
height[1] = 1e-6

for i in 1:(ttest-1)
    biom[i+1] = biom[i] * mumean + biom[i]
    height[i+1] = height[i] *mumean + height[i]
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

layout3 = Layout(
    title = "dry biomass per unit of support area (surface biomass)",
    xaxis_title = "Time (h)",
    yaxis_title = "surface biomass"
)

layout4 = Layout(
    title = "Height of biofilm",
    xaxis_title = "time (h)",
    yaxis_title = "height (m)"
)
plt = plot(
    scatter(y = DLI, x = xplt*1e3, mode = "line", line = attr(color = "blue")), layout1)
display(plt)


plt2 = plot(scatter(y = µ*86400, x = xplt*1e4, mode = "markers"), layout2)
display(plt2)


pp = plot(scatter(y=biom, x = 0:ttest, mode = "markers"), layout3)
display(pp)
ppp = plot(scatter(y=height, x = 0:ttest, mode = "line"))
display(ppp)