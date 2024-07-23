function computelight!(arrLight, I0, ke, dz,pop)
    arrLight[1] = I0
    ind = findfirst(x -> x == 0, pop)
    for i in 2:ind
        arrLight[i] = I0 * exp(-ke*(dz*i))
    end 
end


function updateheight!(pop,µ,dt)
    for i in eachindex(pop)
        pop[i] = pop[i] * µ[i] * dt + pop[i]
    end
end

function grossmu(µ_gross,I,k, sigma, tau, kd, kr, pop)
    ind = findfirst(x -> x == 0, pop)
    for i in 1:ind
        µ_gross[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2)
    end
end

function respiration(R, I, RD, RL, Ik, n, pop)
    ind = findfirst(x -> x == 0, pop)
    for i in 1:ind
        R[i] = RD + (RL - RD) * (I[i]^n/(I[i]^n + Ik^n))
    end
end

function updatemu!(µ, RD, RL, I, n, Ik, k, sigma, tau, kd, kr, pop)
    ind = findfirst(x -> x == 0, pop)
    for i in 1:ind
        R = RD + (RL - RD) * (I[i]^n/(I[i]^n + Ik^n))
        µ[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2) - R
    end
end


function solvematrix(mu, deltaT, nz, b)
    diag = zeros(nz)
    for i in 1:nz
        diag[i] = 1- mu[i] * deltaT
    end
    mat_diag = Diagonal(diag)
    return mat_diag \ b
end

function smootharray!(arr, X)
    for i in eachindex(arr)
        if arr[i] > X
            arr[i+1] = arr[i+1] + arr[i] -X
            arr[i] = X
        end
    end
end

function updateO2!(O, dt, dz, D,pop, Oatm)
    ind = findfirst(x -> x == 0, pop)
    for i in 2:ind
        O[i] = ((O[i-1]-O[i])*D/dz-(O[i]-O[i+1]*D/dz))*dt/dz
    end
    O[1] = ((O[2]-Oatm)*2*D/dz+((O[2]-O[1])*D/dz))*dt/dz
    O[ind] = ((O[ind-1]-O[ind])*D*2/dz)*dt/dz

end