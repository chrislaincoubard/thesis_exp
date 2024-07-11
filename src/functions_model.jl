function computelight!(arrLight, I0, ke, dz, nz,pop)
    arrLight[1] = I0
    for i in 2:nz
        if pop[i] != 0
        arrLight[i] = I0 * exp(-ke*(dz*i))
        end 
    end
    
end

function updateheight!(pop,µ,dt)
    for i in eachindex(pop)
        pop[i] = pop[i] * µ[i] * dt + pop[i]
    end
end

function updatemu!(µ, RD, RL, I, n, Ik, k, sigma, tau, kd, kr, nz, pop)
    for i in 1:nz
        if pop[i] != 0
            R = RD + (RL - RD) * (I[i]^n/(I[i]^n + Ik^n))
            µ[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2) - R
        end
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

function smootharray!(pop, nz, X)
    for i in 1:nz-1
        if pop[i] > X
            pop[i+1] = pop[i+1] + pop[i] -X
            pop[i] = X
        end
    end
end

function updateO2!(O, dt, dz, D)
    if pop != 0
        for i in eachindex(O)[Not 1]
            O[i] = ((O[i-1]-O[i])*D/dz-(O[i]-O[i+1]*D/dz))*dt/dz
        end
        O[1] = ((glu[2]-glu[1])/dz)*D*dt/dz
    end
end