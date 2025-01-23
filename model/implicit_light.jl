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
save_path_data = mkpath(raw"C:\Users\Chrislain\Documents\Results\model_data_test")


@time begin
for I0 in light_intensities
    ##Initialize Arrays for computation and Dataframe for data export
    df_height = DataFrame()
    df = DataFrame()
    df_mu, df_mu_gross, df_R = DataFrame(), DataFrame(), DataFrame() 
    df_O2, df_CO2, df_NO3, df_H2PO4 = DataFrame(), DataFrame(), DataFrame(), DataFrame()
    df_H = DataFrame()
    df_pH = DataFrame()
    filename_O2 = "model_O2_$I0.csv"
    filename_CO2 = "model_CO2_$I0.csv"
    filename_N = "model_N_$I0.csv"
    filename_P = "model_P_$I0.csv"
    filename_µ = "model_mu_$I0.csv"
    filename_µgross = "model_mu_gross_$I0.csv"
    filename_R = "model_R_$I0.csv"
    filename_H = "model_H_$I0.csv"
    filename_pH = "model_pH_$I0.csv"
    light = zeros(nz)
    µ, µ_gross, R = zeros(nz), zeros(nz), zeros(nz)
    mean_mu, height = zeros(tp.t_tot), zeros(tp.t_tot)
    pop = zeros(nz)
    X0 = hmp.rho * dz
    O2, CO2, NO3, H2PO4 = zeros(nz), zeros(nz), zeros(nz), zeros(nz)
    H1, H2 = zeros(nz), zeros(nz)
    pH1, pH2 = zeros(nz), zeros(nz)
    pop[1:25] .= X0
    println("Start for $I0")
    #Starting time (2 loops to save data at specific time points)
    @time for time_step in 1:tp.t_tot
    # @time for time_step in 1:2
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
            SH2PO4 = computeSource(µ, ip.VP_X, gp.Mx, pop, dz)
            ##Get the diagonals from the tridiagonal matrix for the computation
            lowO2, diagO2, upO2 = getdiagonals(gp.D_oxygen, dz, tp.dt, pop)
            lowCO2, diagCO2, upCO2 = getdiagonals(gp.D_CO2, dz, tp.dt, pop)
            lowNO3, diagNO3, upNO3 = getdiagonals_ions(ip.D_NO3, dz, tp.dt, pop)
            lowH2PO4, diagH2PO4, upH2PO4 = getdiagonals(ip.D_H2PO4, dz, tp.dt, pop)
            ##Compute the B matrix for the computation
            B = computeB_gases(O2, gp.O2sat, gp.D_oxygen, dz, tp.dt, pop, SO2)
            BCO2 = computeB_gases(CO2, gp.CO2_surf, gp.D_CO2, dz, tp.dt, pop, SCO2)
            BNO3 = computeB_ions(NO3, ip.C_NO3, ip.D_NO3, dz, tp.dt, pop, SNO3)
            BH2PO4 = computeB_ions(H2PO4, ip.C_H2PO4, ip.D_H2PO4, dz, tp.dt,pop ,SH2PO4)
            LinearAlgebra.LAPACK.gtsv!(lowO2, diagO2, upO2, B)
            LinearAlgebra.LAPACK.gtsv!(lowCO2, diagCO2, upCO2, BCO2)
            LinearAlgebra.LAPACK.gtsv!(lowNO3, diagNO3, upNO3, BNO3)
            LinearAlgebra.LAPACK.gtsv!(lowH2PO4, diagH2PO4, upH2PO4, BH2PO4)
            
            O2[1:length(B)] .= B
            CO2[1:length(BCO2)] .= BCO2
            NO3[1:length(BNO3)] .= BNO3
            H2PO4[1:length(BH2PO4)] .= BH2PO4
            B = make_B(NO3, H2PO4)
            C = make_C(H2PO4, CO2, gp.ka1_co2, gp.ka2_co2, ip.ka1_h3po4, ip.ka2_h3po4, ip.kw)
            for (i,(b, c)) in enumerate(zip(B,C))
                h1,h2 = solve_quadratic(1,b,c)
                H1[i] = h1
                H2[i] = h2
                # pH1[i] = -log(h1/1e-3)
                # pH2[i] = -log(h2/1e-3)
            end
            
        end
        #Append dataframes for each hour of the simulation
        df_H[!,"H1_$time_step"] = copy(H1)
        df_H[!,"H2_$time_step"] = copy(H2)
        df_pH[!, "A_$time_step"] = copy(pH1)
        df_pH[!,"B_$time_step"] = copy(pH2)
        df_O2[!, "$time_step"] = copy(O2)
        df_CO2[!,"$time_step"] = copy(CO2)
        df_NO3[!, "$time_step"] = copy(NO3)
        df_H2PO4[!, "$time_step"] = copy(H2PO4)
        df_mu[!,"$time_step"] = copy(µ)
        df_mu_gross[!, "$time_step"] = copy(µ_gross)
        df_R[!,"$time_step"] = copy(R)
    end
    #Export data in CSV format
    CSV.write(joinpath(save_path_data,filename_O2), df_O2)
    CSV.write(joinpath(save_path_data,filename_CO2), df_CO2)
    CSV.write(joinpath(save_path_data,filename_N), df_NO3)
    CSV.write(joinpath(save_path_data, filename_P), df_H2PO4)
    CSV.write(joinpath(save_path_data, filename_µ),df_mu)
    CSV.write(joinpath(save_path_data, filename_µgross), df_mu_gross)
    CSV.write(joinpath(save_path_data, filename_R), df_R)
    CSV.write(joinpath(save_path_data, filename_H), df_H)
    CSV.write(joinpath(save_path_data, filename_pH), df_pH)
end
end
println("Done with calculation")
println("Data is saved at $(save_path_data)")



