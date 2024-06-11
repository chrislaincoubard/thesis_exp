using Parameters

@with_kw struct TimeParams @deftype Float64
    t_tot = 30*24*3600
    dt = 150
    n_iter = t_tot/dt
    n_save::Int64 = 40
    n_inner = n_iter/n_save
end

@with_kw struct SpaceParams @deftype Float64
    n_glu::Int64 = 100
    n_pop::Int64 = 200
    L_glu = 1e-2
    L_pop = 2e-3
    n_tot::Int64 = n_glu + n_pop
end


@with_kw struct ModelParams @deftype Float64
    mu_0 = 0.29/(34*3600)
    k_s = 0.2
    D_gl = 6.5e-11
    ratio = 0.08
    rho_pop = 1000 * 0.68 * 0.08
end
