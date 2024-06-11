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

function glucose_graph(n_save, n_inner, dt, x, glu_save)
    layout = Layout(
    title = "Glucose concentration",
    xaxis_title = "X * 1000",
    yaxis_title = "Glucose Concentration"
    )

    traces = Vector{GenericTrace}(undef, n_save)
    for i in 1:n_save
        t = i*dt*n_inner/3600
        traces[i] = scatter(x = x*1e3, y = glu_save[i,:], 
            mode = "line", name = "$t h",
            line = attr(color = "red", width = 1))
    end
    plt = plot(traces, layout)
    return plt  
end

function population_graph(time_save, pop_save)
    layout = Layout(
    title = "Bacterial Growth",
    xaxis_title = "Temps (h)",
    yaxis_title = "Height (m)"
    )
    plt = plot(scatter(;x = time_save, y = pop_save, mode = "markers"), layout)
    return plt

end