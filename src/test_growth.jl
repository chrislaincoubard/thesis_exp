using PlotlyJS
using DataFrames
include("model_parameters.jl")

function grossmu!(µ_gross, I, k, sigma, tau, kd, kr)
    for i in eachindex(µ_gross)
        µ_gross[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2)
    end
end

function respiration!(R, I, RD, RL, Ik, n)
    for i in eachindex(R)
        R[i] = RD + (RL-RD)*(I[i]^n/(I[i]^n+Ik^n))
        # R[i] = 0.12/86400
    end
end

function updatemu!(µ,gross_µ, R)
    for i in eachindex(µ)
        µ[i] =  gross_µ[i] - R[i]
    end
end

function computelight!(arrLight, I0, ke, dz)
    arrLight[1] = I0
    for i in eachindex(arrLight)[Not(1)]
        arrLight[i] = I0 * exp(-ke*(dz*i))
    end 
end

function cleanarr(arr)
    return filter(x->x!=0, arr)    
end

function po2!(po2, light, a, b, c)
    for i in eachindex(po2)
        I = light[i]
        po2[i] = I / (a*I^2 + b*I + c)
    end
end

mp = HanModelParams()
z = 5e-4
nz = 500
dz = z/nz
zplt = 0:dz:z
over_light = collect(0:5:2500)

k = 3.6467e-4
kr = 4.8e-4
kd = 2.99e-4
tau = 6.849
sigma = 1.9e-3

PO2 = zeros(length(over_light))

a = 3.42e-7
b = 9.30e-4
c = 2.90e-1
mu = zeros(nz)
mu_gross = zeros(nz)
mu_gross2 = zeros(nz)
mu2 = zeros(nz)
resp = zeros(nz)
light = zeros(nz)
light_200 = zeros(nz)
light_100 = zeros(nz)
light_500 = zeros(nz)
light_1000 = zeros(nz)

computelight!(light, mp.I0, mp.ke, dz)
po2!(PO2, over_light, a, b, c)
# computelight!(light_100, 100, mp.ke, dz)
# computelight!(light_500, 500, mp.ke, dz)
# computelight!(light_1000, 1000, mp.ke, dz)
# grossmu!(mu_gross, over_light, mp.k, mp.sigma, mp.tau, mp.kd, mp.kr)
# grossmu!(mu_gross2, over_light, k, sigma, tau, kd, kr)
# respiration!(resp, over_light, mp.RD, mp.RL, mp.Ik, mp.n)
# updatemu!(mu, mu_gross, resp)
# updatemu!(mu2, mu_gross2, resp)

grossmu!(mu_gross, light, mp.k, mp.sigma, mp.tau, mp.kd, mp.kr)
respiration!(resp, light, mp.RD, mp.RL, mp.Ik, mp.n)
updatemu!(mu, mu_gross, resp)

pp = plot(scatter(x = over_light, y = PO2, mode = "line"), Layout(title = "TEST"))
# clean_oxygen = cleanarr(df[!,"oxygen"])

# pltR = plot([
#     scatter(x = zplt, y = resp.*86400, mode = "line", name = "Respiration", line = attr(color = "orange")),
#     scatter(x = zplt, y = mu_gross.*86400, mode = "line", name = "µ brut", line =attr(color = "blue")),
#     scatter(x = zplt, y = mu.*86400, mode = "line", name = "µ net"),
#     scatter(x = zplt, y = light, mode = "line", name = "Light intensity", yaxis="y2", line = attr(color = "firebrick",dash = "dot")) ],
#     Layout(title = "Growth rate in biofilm",
#      xaxis_title = "Depth (µm)",
#       yaxis_title = "µ (d-1)",
#       yaxis2 = attr(title ="Light intensity (µmol/m2/s)", overlaying = "y", side ="right")))
# display(pltR)

# pltR2 = plot([
#     scatter(x = over_light, y = mu_gross.*86400, mode = "line", name = "µ brut(params Yan)"),
#     scatter(x = over_light, y = mu_gross2.*86400, mode = "line", name = "µ brut (params Wu)")],
#     Layout(title = "Gross growth rate", xaxis_title = "Depth (µm)", yaxis_title = "µ (d-1)"))
# display(pltR2)
# savefig(pltR2,"C:/Users/chris/OneDrive/Images/plots/comparaison_mu_params.png" )
# lightpl = plot(
#     scatter(x = zplt*10^6, y=light, mode = "markers", name = "200"),
#     Layout(title = "light attenuation", xaxis_title="depth (µm)",
#      yaxis_title = "Light intensity (µmol-1/m²/s)",
#      legend_title_text="Incident light intensity" 
#      ))

# display(lightpl)
# savefig(pltR, "C:/Users/chris/OneDrive/Images/plots/growth_200_Yan_plus_light.png")
