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

function grossmu!(µ_gross, I, k, sigma, tau, kd, kr, pop)
    ind = findfirst(x->x==0, pop)
    for i in 1:ind
        µ_gross[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2)
    end
end

function respiration!(R, I, RD, RL, Ik, n, pop)
    ind = findfirst(x->x==0, pop)
    for i in 1:ind
        R[i] = RD + (RL-RD)*(I[i]^n/(I[i]^n+Ik^n))
        # R[i] = 0.12/86400
    end
end

function updatemu!(µ,gross_µ, R, pop)
    ind = findfirst(x -> x == 0, pop)
    for i in 1:ind
        µ[i] =  gross_µ[i] - R[i]
    end
end


function solvematrix(mu, deltaT, b)
    diag = 1 .- mu * deltaT
    mat_diag = Diagonal(diag)
    return mat_diag \ b
end

function smootharray!(arr, X)
    n = length(arr)
    for i in 1:n-1
        if arr[i] > X
            excess = arr[i] - X
            arr[i] = X
            arr[i+1] += excess
        end
    end
end

function getdiagonals(D, dz, dt, pop)
    ind = findfirst(x -> x == 0, pop) -1
    coefDiff = D*dt/dz^2
    diag = fill(1 + 2 * coefDiff, ind)
    up = fill(-coefDiff, ind-1)
    low = fill(-coefDiff, ind-1)
    diag[1] = 1+3*coefDiff
    diag[ind] = 1+coefDiff
    return low, diag, up
end



function computeO2source(mu, VO2x, mx, pop, dz)
    ind = findfirst(x -> x == 0, pop)-1 
    source = zeros(ind)
    for i in 1:ind
        source[i] = ((mu[i]*(pop[i]/dz)*VO2x)/mx)
    end
    return source
end

function computeB(O2, O2surf, D, dz, dt, pop,S)
    ind = findfirst(x -> x == 0, pop) -1 
    B = O2[1:ind] .+ S[1:ind] .* dt
    # B = O2[1:ind]
    B[1] += 2 * dt * D * O2surf / dz^2
    return B
end

function calcexcess(arr, X)
    n = length(arr)
    total_excess = 0
    for i in 1:n-1
        if arr[i] > X
        total_excess += arr[i] - X
        end
    end
    return total_excess
end
