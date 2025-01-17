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
save_path_data = mkpath(raw"C:\Users\Chrislain\Documents\Results\model_data")
save_path = mkpath(raw"C:\Users\Chrislain\Documents\Results\model_O2_clean_test")
save_path_CO2 = mkpath(raw"C:\Users\Chrislain\Documents\Results\model_CO2_clean_test")
save_path_N = mkpath(raw"C:\Users\Chrislain\Documents\Results\model_N03_test")


@time begin
for I0 in light_intensities
    ##Initialize Array
    df_height = DataFrame()
    df = DataFrame()
    df_mu, df_mu_gross, df_R = DataFrame(), DataFrame(), DataFrame() 
    df_O2, df_CO2, df_NO3, df_PO2 = DataFrame(), DataFrame(), DataFrame(), DataFrame()
    filename_O2 = "model_O2_$I0.csv"
    filename_CO2 = "model_CO2_$I0.csv"
    filename_N = "model_N_$I0.csv"
    filename_P = "model_P_$I0.csv"
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
            # println("Start $i_inner for $time_step" )
            computelight!(light, I0, hmp.ke, dz, pop)
            grossmu!(µ_gross, light, hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, pop)
            respiration!(R, light, hmp.RD, hmp.RL, hmp.Ik, hmp.n, pop)
            updatemu!(µ, µ_gross, R, pop)
            pop .= solvematrix(µ, tp.dt, pop)
            smootharray!(pop, X0)
            ## Compute gases ##
            SO2 = computeSource(µ, gp.VO2_x, gp.Mx, pop, dz)
            SCO2 = computeSource(µ, gp.VCO2_x, gp.Mx, pop, dz) 
            SNO3 = computeSource(µ, ip.VN_X, gp.Mx, pop, dz)
            lowO2, diagO2, upO2 = getdiagonals(gp.D_oxygen, dz, tp.dt, pop)
            lowCO2, diagCO2, upCO2 = getdiagonals(gp.D_CO2, dz, tp.dt, pop)
            lowNO3, diagNO3, upNO3 = getdiagonals_ions(ip.D_NO3, dz, tp.dt, pop)
            #lowH2PO4, diagH2PO4, upH2PO4 = getdiagonals(gp.D_PO4, dz, tp.dt, pop)
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
            # Optionnal plots (comment for speed) ##
            
            # # Optionnal plots (comment for speed) ##
            # if i_inner in tp.n_inner #&& time_step == tp.n_save
            #     pp = plot(scatter(x = zplt*10^6, y = CO2, mode = "line"), 
            #     Layout(title = "C02 concentration profile $time_step",
            #     xaxis_title = "Depth (µm)",
            #     yaxis_title = "CO2 concentration mol/m3"))
            #     display(pp)
            # end
            
        end
        # p = plot(scatter(x = eachindex(SO2), y = SO2, mode = "line" ))
        # display(p)
        
    
        df_O2[!, "$time_step"] = copy(O2)
        df_CO2[!,"$time_step"] = copy(CO2)
        df_NO3[!, "$time_step"] = copy(NO3)
        cleanO2 = removezeros(O2)
        cleanmu = removezeros(µ)
        cleanCO2 = removezeros(CO2)
        cleanN = removezeros(NO3)
        file_name_O2 = "O2_profile_$(time_step)_h.png"
        file_name_CO2 = "CO2_profile_$(time_step)_h.png"
        file_name_NO3 = "NO3_profile_$(time_step)_h.png"
        # pp = plot([scatter(x = zplt*10^6, y = cleanO2, mode = "line", name = "O2"),
        # scatter(x = zplt*10^6, y = cleanmu.*86400, mode = "line", name = "mu", yaxis = "y2")], 
        # Layout(title = "02 concentration profile $(time_step*i_inner/72) h",
        # xaxis_title = "Depth (µm)",
        # yaxis_title = "O2 concentration mol/m3",
        # xaxis_range = [0,325],
        # yaxis_range = [0.25,0.45],
        # yaxis2 = attr(title = "Growth rate", overlaying = "y", side = "right")))

        # pp = plot([scatter(x = zplt*10^6, y = cleanO2, mode = "line", name = "O2")], 
        # Layout(title = "02 concentration profile for 200 µmol<sub>photons</sup>/m<sup>2</sup>/s after $(time_step) hours",
        # xaxis_title = "Depth (µm)",
        # yaxis_title = "O2 concentration mol/m<sup>3</sup>",
        # xaxis_range = [0,325],
        # yaxis_range = [0.25,0.45]))
        # display(pp)

        # p = plot(scatter(x = zplt*10^6, y= cleanCO2, mode = "line"),
        # Layout(title = "CO2 concentration profile $(time_step) h",
        # xaxis_title = "Depth (µm)",
        # yaxis_title = "CO2 concentration mol/m3",
        # # xaxis_range = [0,325],
        # # yaxis_range = [0, 0.28]
        # ))
        # # display(p)
        # p_n = plot(scatter(x = zplt*10^6, y= cleanN, mode = "line"),
        # Layout(title = "NO3 concentration profile $(time_step) h",
        # xaxis_title = "Depth (µm)",
        # yaxis_title = "NO3 concentration mol/m3",
        # # xaxis_range = [0,325],
        # # yaxis_range = [0, 0.28]
        # ))
        # display(pp)
        # savefig(p, joinpath(save_path_CO2, file_name_CO2))
        # savefig(pp, joinpath(save_path, file_name_O2))
        # savefig(p_n, joinpath(save_path_N, file_name_NO3))
        time = time_step
        currheight = sum(pop) / hmp.rho
        time_save[time_step] = time
        height[time_step] = currheight
        mu_no_zeros = removezeros(µ)
        mean_mu[time_step] = mean(mu_no_zeros)        
    end
    #save data at some time points in DataFrame
    
    CSV.write(joinpath(save_path_data,filename_O2), df_O2)
    CSV.write(joinpath(save_path_data,filename_CO2), df_CO2)
    CSV.write(joinpath(save_path_data,filename_N), df_NO3)
    clean_R = removezeros(R)
    clean_µ = removezeros(µ)
    clean_grossµ = removezeros(µ_gross)
    clean_light = removezeros(light)
    df[!,"mu_net"] = clean_µ
    df[!,"gross_mu"] = clean_grossµ
    df[!,"R"] = clean_R
    df[!,"light"] = clean_light
    df[!,"Intensity"] .= I0
    df_height[!,"Height"] = height
    df_height[!,"Intensity"] .= I0
    df_height[!, "time"] .= time_save
    df_height[!,"mean_mu"] .= mean_mu

    #Export dataframe
    filename = "data_model_$(I0)_clean.csv" 
    filename2 = "height_model_$(I0)_clean.csv"
    path = joinpath(raw"C:\Users\Chrislain\Documents\Results\result_model",filename)
    path2 = joinpath(raw"C:\Users\Chrislain\Documents\Results\result_model",filename2)
    CSV.write(path, df)
    CSV.write(path2, df_height)
end
end
println("Done with calculation")
println("Data is saved at $(save_path_data)")



