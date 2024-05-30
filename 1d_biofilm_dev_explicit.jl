#1D simulation of biofilm development (method : Finite Volume, explicit)
#Author : Chrislain Coubard (adapted from Patrick Perre)
#Date : 05/2024


using PlotlyJS

function UpdateValues!(n_glu, n_tot, mu_0, k_s, dt, ratio, rho_pop,  glu, pop, s_glu)
    for  i_x in (1+n_glu):n_tot
        mu = mu_0 * glu[i_x] / (k_s + glu[i_x])
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

function UpdateDGlu!(n_tot, d_glu,dt, glu, D, dx, deltax)
    for i_x in 2:n_tot-1      
        d_glu[i_x] = ((glu[i_x-1]-glu[i_x])*D[i_x-1]/dx[i_x-1]-(glu[i_x]-glu[i_x+1])*D[i_x]/dx[i_x])*dt/deltax[i_x]
    end
    d_glu[1] = (-(glu[1]-glu[2])/dx[1])*D[1]*dt/deltax[1]
    d_glu[n_tot] = 0
end

function buildScatter(x, y, n_save)
    arr_traces = Vector{GenericTrace}
    for i in 1:n_save
        println("y[i,:]", y[i,:])
        trace = scatter(;x=x, y=y[i,:],mode ="line")
        push!(arr_traces, trace)
    end
    return plot(arr_traces)
end

#Init values
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
glu_save = zeros(n_save, n_tot)
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
for i in 1:n_tot-1
    dx[i] = (deltax[i] + deltax[i+1])/2
end
glu = zeros(n_tot)
d_glu = zeros( n_tot)
pop = zeros(n_tot)
s_glu = zeros( n_tot)

glu[1:n_glu] .= 10
ind = [1,2,3,4,5] .+ n_glu 
pop[ind] .= 1



for time_step in 1:n_save #save points
    for i_inner in 1:n_inner #actual computation
        UpdateValues!(n_glu, n_tot, mu_0, k_s, dt, ratio, rho_pop, glu, pop, s_glu)
        SmoothArr!(n_glu, n_tot, pop)
        UpdateD!(n_glu, n_tot, D_gl, D, pop)
        UpdateDGlu!(n_tot, d_glu, dt, glu, D, dx, deltax)
        @. glu = glu + d_glu + s_glu
    end
    glu_save[time_step,:] = glu
    pop_save[time_step] = sum(pop)
    time_save[time_step] = dt*time_step*n_inner/3600
end



plot(scatter(;x = time_save,y= pop_save, mode = "markers"))


# traces = Vector{GenericTrace}(undef, n_save)
# for i in 1:n_save
#     traces[i] = scatter(x = x*1e3, y = glu_save[i,:], mode = "line", name = "$i",line_color = "red")
# end
# plot(traces)