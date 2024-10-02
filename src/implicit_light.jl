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
light_intensities = [300]
z = 1e-3
nz = Int64(1e3)
dz = z/nz
zplt = 0:dz:z
time_save = zeros(tp.n_save)

@time begin
for I0 in light_intensities
    ##Initialize Array
    df_height = DataFrame()
    df = DataFrame()
    light = zeros(nz)
    µ = zeros(nz)
    µ_gross = zeros(nz)
    mean_mu = zeros(tp.n_save)
    R = zeros(nz)
    pop = zeros(nz)
    height = zeros(tp.n_save)
    X0 = hmp.rho * dz
    O2 = zeros(nz)
    CO2 = zeros(nz)
    CO2sat = gp.PCo2 * gp.HCo2 * (1+gp.ka1/10.0^(-gp.phs)+(gp.ka1*gp.ka2/(10.0^(-gp.phs))^2))
    pop[1:25] .= X0
    print("sum pop start = $(sum(pop))")
    println("Start for $I0")
    #Starting time (2 loops to save data at specific time points)
    @time for time_step in 1:tp.n_save
        println("Start $time_step")
        for i_inner in 1:tp.n_inner
            ## Compute growth ##
            computelight!(light, I0, hmp.ke, dz, pop)
            grossmu!(µ_gross, light, hmp.k, hmp.sigma, hmp.tau, hmp.kd, hmp.kr, pop)
            respiration!(R, light, hmp.RD, hmp.RL, hmp.Ik, hmp.n, pop)
            updatemu!(µ, µ_gross, R, pop)
            pop .= solvematrix(µ, tp.dt, pop)
            smootharray!(pop, X0)
            ## Compute gases ##
            SO2 = computeO2source(µ, gp.VO2_x, gp.Mx,pop, dz)
            SCO2 = .-SO2
            low, diag, up = getdiagonals(gp.D_oxygen, dz, tp.dt, pop)
            lowCO2, diagCO2, upCO2 = getdiagonals(gp.D_CO2, dz, tp.dt, pop)
            B = computeB(O2, gp.O2sat, gp.D_oxygen, dz, tp.dt, pop, SO2)
            
            BCO2 = computeB(CO2, CO2sat, gp.D_CO2, dz, tp.dt, pop, SCO2)
            LinearAlgebra.LAPACK.gtsv!(low, diag, up, B)
            LinearAlgebra.LAPACK.gtsv!(lowCO2, diagCO2, upCO2, BCO2)
            O2[1:length(B)] .= B
            CO2[1:length(BCO2)] .= BCO2
            # Optionnal plots (comment for speed) ##
            if i_inner in tp.n_inner #&& time_step == tp.n_save
                p = plot(scatter(x = eachindex(B), y = B, mode = "line" ))
                display(p)
                pp = plot(scatter(x = zplt*10^6, y = O2, mode = "line"), 
                Layout(title = "02 concentration profile $time_step",
                xaxis_title = "Depth (µm)",
                yaxis_title = "O2 concentration mol/m3"))
                display(pp)
            end
            # # Optionnal plots (comment for speed) ##
            # if i_inner in tp.n_inner #&& time_step == tp.n_save
            #     pp = plot(scatter(x = zplt*10^6, y = CO2, mode = "line"), 
            #     Layout(title = "C02 concentration profile $time_step",
            #     xaxis_title = "Depth (µm)",
            #     yaxis_title = "CO2 concentration mol/m3"))
            #     display(pp)
            # end
        end
        time = tp.dt*time_step*tp.n_inner/3600
        currheight = sum(pop) / hmp.rho
        time_save[time_step] = time
        height[time_step] = currheight
        mu_no_zeros = removezeros(µ)
        mean_mu[time_step] = mean(mu_no_zeros)
        
    end
    #save data at some time points in DataFrame
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
    println(length(mean_mu))
    println(length(height))
    df_height[!,"mean_mu"] .= mean_mu

    #Export dataframe
    filename = "data_model_$I0.csv" 
    filename2 = "height_model_$I0.csv"
    path = joinpath(raw"C:\Users\Chrislain\Documents\Results\result_model",filename)
    path2 = joinpath(raw"C:\Users\Chrislain\Documents\Results\result_model",filename2)
    CSV.write(path, df)
    CSV.write(path2, df_height)
end
end
println("Done with calculation")
println(raw"Data is saved at C:\Users\Chrislain\Documents\Results\result_model")



