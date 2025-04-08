using Parameters

@with_kw struct TimeParams @deftype Number
    t_tot = 7*24 #(h)
    dt = 1 #(s)
    n_inner::Int = 3600/dt
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
    kd = 0.15 #damage rate (dimentionless)
    tau = 1.06e-3 #turnover time electron transfer chain (s-1)
    sigma = 1.59e-2 #cross section of PSII (m2*µmol-1)
    n = 3 #switching rate between RL and RD (dimentionless)
    Ik = 50 #half saturating light for light respiration (µmol*m-2*s-1)
    RD = 9.15e-8 #respiration in dark phase (s-1)
    RL = 1.78e-6 #respiration in light phase (s-1)
    rho = 1.4e2  #biomass volumic mass (kg*m-3)
end

@with_kw struct GasesParams @deftype Number
    D_oxygen = 0.68*1.99e-9 #coef de diffusion de l'O2 dans le biofilm = 68% du coef de diffusion de l'O2 dans l'eau (m2/s)
    D_CO2 = 1.9e-9#coefficient de diffusion du CO2 assimilé à celui des ions bicarbonate
    O2sat = 0.274 #Valeur d'un milieu solide saturé en 02, assimilé au biofilm (mol/m3)
    H_constant_CO2 = 3.4e-2 #Henry's constant for CO2 in water at 25°C (mol/L/atm)
    PP_CO2 = 4.2e-4 #partial pressure of CO2 in air (atm)
    CO2_surf = H_constant_CO2 * PP_CO2 * 1000 #(mol/m3)
    VO2_x = 1.125 #Coefficient stoechiométrique de production d'O2 par rapport à la biomasse
    VCO2_x = -1 #Coefficient stoechiométrique de consommation du CO2 par rapport à la biomasse
    Mx = 0.024 # C-masse molaire de la biomasse (kgx/C-mol)
    KI = 4.47e-7 #dissociation coefficient of CO2
    KII = 4.67e-11 #dissociation coefficient of HCO3-
    
end

@with_kw struct Ions_params @deftype Number
    D_NO3 = 1.701e-9 #Diffusion coefficient of nitrate in water (m2/s)
    D_H2PO4 = 7.159e-10 #Diffusion coefficient of phosphate ion in water (m2/s)
    C_NO3 = 8.8 #Concentration de NO3-dans le milieu (mol/m3)
    C_H2PO4 = 1.71 #Concentration de H2PO4- dans le milieu de culture (mol/m3)
    C_Na = 8.8 #Concentration de Na+ dans le milieu (mol/m3)
    C_K = 2.14 #Concentration de K+ le milieu de culture (mol/m3)
    QI = 6.2e-8 #dissociation coefficient of H2PO4- into HPO42-
    QII = 4.8e-13 #dissociation coefficient of HPO42- into PO43-
    kw = 1e-14 #dissociation coefficient of H20
    VN_X = -0.16 #Coefficient stoechiométrique de consommation de NO3- par rapport à la biomasse
    VP_X = -0.007 #Coefficient stoechiométrique de consommation de H2PO42- par rapport à la biomasse

end