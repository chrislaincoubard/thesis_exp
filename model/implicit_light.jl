using PlotlyJS
using LinearAlgebra
using ColorSchemes
using DataFrames
using ColorTypes
using CSV
using Glob
using Statistics
using DelimitedFiles
using LinearAlgebra.LAPACK
include("model_parameters.jl")
include("functions_model.jl")
include("utils.jl")

## Values of parameters in the "model_parameters.jl" file
## Functions to run the model in "functions_model.jl" file
tp = TimeParams()
hmp = HanModelParams() 
gp = GasesParams()
ip = Ions_params()
light_intensities = [200]
z = 1e-3
nz = Int64(1e3)
dz = z/nz
zplt = 0:dz:z
time_save = zeros(tp.t_tot)

#To change depending of the machine
save_path_data = mkpath(raw"C:\Users\chris\Documents\Results\model_data")


@time begin
for I0 in light_intensities
    ##Initialize Arrays for computation and Dataframe for data export
    df_height = DataFrame()
    df = DataFrame()
    df_mu, df_mu_gross, df_R = DataFrame(), DataFrame(), DataFrame() 
    df_O2, df_CO2, df_NO3, df_PO2 = DataFrame(), DataFrame(), DataFrame(), DataFrame()
    filename_O2 = "model_O2_$I0.csv"
    filename_CO2 = "model_CO2_$I0.csv"
    filename_N = "model_N_$I0.csv"
    filename_P = "model_P_$I0.csv"
    filename_µ = "model_mu_$I0.csv"
    filename_µgross = "model_mu_gross_$I0.csv"
    filename_R = "model_R_$I0.csv"
    light = zeros(nz)
    µ, µ_gross, R = zeros(nz), zeros(nz), zeros(nz)
    mean_mu, height = zeros(tp.t_tot), zeros(tp.t_tot)
    pop = zeros(nz)
    X0 = hmp.rho * dz
    O2, CO2, NO3, H2PO4 = zeros(nz), zeros(nz), zeros(nz), zeros(nz)
    pop[1:25] .= X0
    println("Start for $I0")
    #Starting time (2 loops to save data at specific time points)
    @time for time_step in 1:tp.t_tot
        # println("Start $time_step")
        for i_inner in 1:tp.n_inner
            ## Compute growth ##
            computelight!(light, I0, hmp.ke, dz, pop)
            grossmu!(µ_gross, light, hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, pop)
            respiration!(R, light, hmp.RD, hmp.RL, hmp.Ik, hmp.n, pop)
            updatemu!(µ, µ_gross, R, pop)
            pop .= solvematrix(µ, tp.dt, pop)
            smootharray!(pop, X0)
            ## Compute gases and ions##
            SO2 = computeSource(µ, gp.VO2_x, gp.Mx, pop, dz)
            SCO2 = computeSource(µ, gp.VCO2_x, gp.Mx, pop, dz) 
            SNO3 = computeSource(µ, ip.VN_X, gp.Mx, pop, dz)
            ##Get the diagonals from the tridiagonal matrix for the computation
            lowO2, diagO2, upO2 = getdiagonals(gp.D_oxygen, dz, tp.dt, pop)
            lowCO2, diagCO2, upCO2 = getdiagonals(gp.D_CO2, dz, tp.dt, pop)
            lowNO3, diagNO3, upNO3 = getdiagonals_ions(ip.D_NO3, dz, tp.dt, pop)
            #lowH2PO4, diagH2PO4, upH2PO4 = getdiagonals(gp.D_PO4, dz, tp.dt, pop)
            ##Compute the B matrix for the computation
            B = computeB_gases(O2, gp.O2sat, gp.D_oxygen, dz, tp.dt, pop, SO2)
            BCO2 = computeB_gases(CO2, gp.CO2_surf, gp.D_CO2, dz, tp.dt, pop, SCO2)
            BNO3 = computeB_ions(NO3, ip.C_NO3, ip.D_NO3, dz, tp.dt, pop, SNO3)
            #BH2PO4 = computeB(H2PO4, C_PO4, gp.D_H2PO4, dz, tp.dt, SH2PO4)
            LinearAlgebra.LAPACK.gtsv!(lowO2, diagO2, upO2, B)
            LinearAlgebra.LAPACK.gtsv!(lowCO2, diagCO2, upCO2, BCO2)
            LinearAlgebra.LAPACK.gtsv!(lowNO3, diagNO3, upNO3, BNO3)
            #LinearAlgebra.LAPACK.gtsv!(lowH2PO4, diagH2PO4, upH2PO4, BH2PO4)
            O2[1:length(B)] .= B
            CO2[1:length(BCO2)] .= BCO2
            NO3[1:length(BNO3)] .= BNO3
            
        end
        #Append dataframes for each hour of the simulation
        df_O2[!, "$time_step"] = copy(O2)
        df_CO2[!,"$time_step"] = copy(CO2)
        df_NO3[!, "$time_step"] = copy(NO3)
        df_mu[!,"$time_step"] = copy(µ)
        df_mu_gross[!, "$time_step"] = copy(µ_gross)
        df_R[!,"$time_step"] = copy(R)       
    end
    #Export data in CSV format
    
    CSV.write(joinpath(save_path_data,filename_O2), df_O2)
    CSV.write(joinpath(save_path_data,filename_CO2), df_CO2)
    CSV.write(joinpath(save_path_data,filename_N), df_NO3)
    CSV.write(joinpath(save_path_data, filename_µ),df_mu)
    CSV.write(joinpath(save_path_data, filename_µgross), df_mu_gross)
    CSV.write(joinpath(save_path_data, filename_R), df_R)
end
end
println("Done with calculation")
println("Data is saved at $(save_path_data)")



