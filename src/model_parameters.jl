using Parameters

@with_kw struct TimeParams @deftype Number
    t_tot = 7*86400 #(s)
    dt = 150
    n_iter::Int = round(Int, t_tot/dt)
    n_save = 32
    n_inner::Int = n_iter/n_save
end

@with_kw struct SpaceParams @deftype Number
    x = 1e-4 #total height of biofilm
    nx = 500
    dx = x/nx
end

@with_kw struct HanModelParams @deftype Number
    I0 = 200 #Initial Light intensity (µmol*m-2*s-1)
    ke = 11650 #biofilm light extinction coefficient (m-1)
    k = 4.9e-6 #scaling factor, linking the rate of photons capture (s*s-1) 
    kr = 1.16e-3 #repair rate (s-1)
    kd = 0.15#damage rate (dimentionless)
    tau = 1.06e-3 #turnover time electron transfer chain (s-1)
    sigma = 1.59e-2 #cross section of PSII (m2*µmol-1)
    n = 3 #switching rate between RL and RD (dimentionless)
    Ik = 50 #half saturating light for light respiration (µmol*m-2*s-1)
    RD = 9.15e-8 #respiration in dark phase (s-1)
    RL = 1.78e-6 #respiration in light phase (s-1)
    rho = 1.4e2  #biomass volumic mass (kg*m-3)
end

@with_kw struct GasesParams @deftype Number
    D_oxygen = 0.68*1.99e-9 #coef de diffusion de l'O2 dans le biofilm = 68% du coef de diffusion de l'O2 dans l'eau (m/s)
    O2atm = 0.274 #Valeur d'un milieu solide saturé en 02, assimilé au biofilm (mol/m3)
    VO2_x = 1.125 #Coefficient stoechiométrique de production d'O2 par rapport à la biomasse
    Mx = 0.024 # C-masse molaire de la biomasse (kgx/C-mol)
end