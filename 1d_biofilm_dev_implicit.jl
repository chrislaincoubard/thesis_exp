#1D simulation of biofilm development (method : Finite Volume, implicit)
#Author : Chrislain Coubard (adapted from Patrick Perre)
#Date : 05/2024

using PlotlyJS
using LinearAlgebra

function UpdateValues!(n_glu, n_tot, mu_0, k_s, dt, ratio, rho_pop,  glu, pop, s_glu)
    for  i_x in (1+n_glu):n_tot
        mu = (mu_0 * glu[i_x]) / (k_s + glu[i_x])
        d_pop = pop[i_x] * mu * dt
        pop[i_x] = pop[i_x] + d_pop
        s_glu[i_x] = -d_pop / ratio * rho_pop
    end
end

function SmoothArr!(n_glu, n_tot, pop)
    for i_x in (1+n_glu):n_tot
        if pop[i_x] > 1
            pop[i_x+1] = pop[i_x+1] + pop[i_x] - 1
            pop[i_x] = 1
        end
    end
end

function UpdateD!(n_glu, n_tot, D_gl, D, pop)
    for i_x in (1+n_glu):n_tot
        D[i_x] = D_gl * pop[i_x]
    end
end

function MakeTriDiag(D, n_tot, dx, dt, deltax,up, diag, low)

    up[1] = -D[1]*dt/dx[1]/deltax[1] #Flux nul
    diag[1] = 1 - up[1]
    low[n_tot-1] = -D[n_tot-1] * dt/dx[n_tot-1]/deltax[n_tot]
    diag[n_tot] = 1 - low[n_tot-1]
    
    for i_x in 2:(n_tot-1)
        low[i_x-1] = -D[i_x-1] * dt/dx[i_x-1]/deltax[i_x]
        up[i_x] = -D[i_x] * dt/dx[i_x]/deltax[i_x]
        diag[i_x] = 1 - up[i_x] - low[i_x-1]    
    end
    mat = Tridiagonal(low, diag, up)
    return mat
end

function UpdateGlu!(B, glu, s_glu, tridiag)
    @. B = glu + s_glu
    glu .= tridiag \ B
end
#Init values

n_glu = 100
n_pop = 200
L_glu = 1e-2 #gelose
L_pop = 2e-3 #culture

n_tot = n_glu + n_pop

mu_0 = 0.29/(24*3600)
k_s = 0.2
dt = 150
t_tot = 30*24*3600
D_gl = 6.5e-11
ratio = 0.08
rho_pop = 1000 * 0.68 * 0.08

n_iter = Int64(t_tot/dt)
n_save = 40
n_inner = n_iter/n_save

pop_save = zeros(n_save)
time_save = zeros(n_save)
glu_save = zeros(n_save, n_tot)
time = zeros(n_save)

deltax_glu = fill(L_glu/n_glu, n_glu)  # size of CVs regular mesh, full CVs at boundary, uniform by zone
deltax_pop = fill(L_pop/n_pop, n_pop)
deltax = vcat(deltax_glu, deltax_pop)
x = zeros(n_tot)
dx = zeros(n_tot-1)
D = fill(D_gl, n_tot)

for i in 1:n_glu
    x[i] = deltax_glu[1]*(i-0.5)
end
for i in 1:n_pop
    x[i+n_glu] = L_glu + deltax_pop[1] * (i-0.5)
end
for i in 1:(n_tot-1)
    dx[i] = (deltax[i] + deltax[i+1])/2
end

glu = zeros(n_tot)
d_glu = zeros(n_tot)
pop = zeros(n_tot)
s_glu = zeros(n_tot)

glu[1:n_glu] .= 10
ind = [1,2,3,4,5] .+ n_glu 
pop[ind] .= 1

# 1D simulation
diag = zeros(n_tot)
low = zeros(n_tot-1)
up = zeros(n_tot-1)
B = zeros(n_tot)


for time_steps in 1:n_save
    for i_inner in 1:n_inner
        UpdateValues!(n_glu, n_tot, mu_0, k_s, dt, ratio, rho_pop, glu, pop, s_glu)
        SmoothArr!(n_glu, n_tot, pop)
        UpdateD!(n_glu, n_tot, D_gl, D, pop)
        tridiag = MakeTriDiag(D, n_tot, dx, dt, deltax,up, diag, low)
        UpdateGlu!(B, glu, s_glu, tridiag)
    end
    glu_save[time_steps,:] = glu
    pop_save[time_steps] = sum(pop)
    time_save[time_steps] = dt*time_steps*n_inner/3600
end

layout1 = Layout(
    title = "Bacterial Growth",
    xaxis_title = "Temps (h)",
    yaxis_title = "Height"
    )

layout2 = Layout(
    title = "Glucose concentration",
    xaxis_title = "X * 1000",
    yaxis_title = "Glucose Concentration"
)
plt1 = plot(scatter(;x = time_save, y = pop_save, mode = "markers"), layout1)

traces = Vector{GenericTrace}(undef, n_save)
for i in 1:n_save
    t = i*dt*n_inner/3600
    traces[i] = scatter(x = x*1e3, y = glu_save[i,:], 
        mode = "line", name = "$t h",
        line = attr(color = "red", width = 1))
end
plt2 = plot(traces, layout2)
display(plt1)
display(plt2)