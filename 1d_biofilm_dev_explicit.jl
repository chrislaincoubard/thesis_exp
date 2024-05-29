#1D simulation of biofilm development (method : Finite Volume, explicit)
#Author : Chrislain Coubard (adapted from Patrick Perre)
#Date : 05/2024

#Init values
using PlotlyJS

n_glu = 100
n_pop = 100
L_glu = 1e-2
L_pop = 1e-3

n_tot = n_glu + n_pop

mu_0 = 0.29/(24*3600)
k_s = 0.2
dt = 1
t_tot = 15*24*3600
D_gl = 5e-11
dx = zeros(n_tot-1)
ratio = 0.08
rho_pop = 1000 * 0.68 *0.08

n_iter = Int64(t_tot/dt)
n_save = 40
n_inner = n_iter/n_save

pop_save = zeros(n_save)
time_save = zeros(n_save)
times = zeros(n_save)

deltax_glu = fill(L_glu/n_glu,n_glu)
deltax_pop = fill(L_pop/n_pop, n_pop)
deltax = vcat(deltax_glu, deltax_pop)
D = fill(D_gl, n_tot)



x = zeros(n_tot)
for j in 1:n_glu
    x[j] = deltax_glu[1] * (j-0.5)
end
for h in 1:n_pop
    x[h+n_glu] = L_glu + deltax_pop[1]*(h-0.5)
end
println(deltax)
for i in 1:n_tot-1
    dx[i] = (deltax[i] - deltax[i+1])/2
end
println(dx)
glu = zeros(n_iter,n_tot)
d_glu = zeros(n_iter, n_tot)
pop = zeros(n_iter, n_tot)
s_glu = zeros(n_iter, n_tot)

glu[1,1:n_glu] .= 10
glu
ind = [1,2,3,4,5] .+ n_glu 
pop[1,ind] .= 1

# traces_arr = []
# p = plot(scatter(;x=x*1e3, y=glu, mode="line"))
# push!(traces_arr,p)
# addtraces(p1, scatter(x = (x2*1e3), y = glu), mode = "markers")
for i in 1:n_iter
    for j in (1+n_glu):n_tot
        mu = mu_0 * glu[i,j]/(k_s + glu[i,j])
        d_pop = pop[i,j] * mu * dt
        pop[i,j] = pop[i,j] + d_pop
        s_glu[i,j] = -d_pop / ratio * rho_pop
    end
    for j in (1+n_glu):n_tot
        if pop[i,j] > 1
            pop[i, j+1] = pop[i,j+1] + pop[i, j] -1
            pop[i,j] = 1
        end
    end
    for j in (n_glu+1):n_tot
        D[i,j] = D_gl * pop[i,j]
    end
    for j in 2:(n_tot-1)
        d_glu[i,j] = ((glu[i,j-1]-glu[i,j])*D[i,j-1]/dx[i,j-1]-(glu[i,j]-glu[i,j+1])*D[i,j]/dx[i,j])*dt/deltax[i,j]
    end
    d_glu[i,1] = (-(glu[i,1]-glu[i,2])/dx[1])*D[i,1]*dt/deltax[1]
    d_glu[n_tot] = 0
    @. glu = glu + d_glu + s_glu
end


