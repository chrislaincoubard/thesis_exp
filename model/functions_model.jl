function computelight!(arrLight, I0, ke, dz, ind)
    arrLight[1] = I0
    for i in 2:ind
        arrLight[i] = I0 * exp(-ke*(dz*i))
    end 
end


function updateheight!(pop,µ,dt)
    for i in eachindex(pop)
        pop[i] = pop[i] * µ[i] * dt + pop[i]
    end
end

function grossmu!(µ_gross, I, k, sigma, tau, kd, kr, ind)
    for i in 1:ind
        µ_gross[i] = (k*sigma*I[i]) / (1 + tau*sigma*I[i] + (kd/kr)*tau*(sigma*I[i])^2)
    end
end

function respiration!(R, I, RD, RL, Ik, n, ind)
    for i in 1:ind
        R[i] = RD + (RL-RD)*(I[i]^n/(I[i]^n+Ik^n))
        # R[i] = 0.12/86400
    end
end

function updatemu!(µ,gross_µ, R, ind)
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

function getdiagonals(D, dz, dt, ind)
    coefDiff = D*dt/dz^2
    diag = fill(1 + 2 * coefDiff, ind)
    up = fill(-coefDiff, ind-1)
    low = fill(-coefDiff, ind-1)
    diag[1] = 1+3*coefDiff
    diag[ind] = 1+3*coefDiff
    return low, diag, up
end

function getdiagonals_ions(D, dz, dt, ind)
    coefDiff = D*dt/dz^2
    diag = fill(1 + 2 * coefDiff, ind)
    up = fill(-coefDiff, ind-1)
    low = fill(-coefDiff, ind-1)
    diag[1] = 1+coefDiff
    diag[ind] = 1+3*coefDiff
    return low, diag, up
end


function computeSource(mu, Stoech, mx, dz,pop, ind)
    source = zeros(ind)
    for i in 1:ind
        source[i] = ((mu[i]*(pop[i]/dz)*Stoech)/mx)
    end
    return source
end


function computeB_gases(Phi, Phi_surf, D, dz, dt, S, ind)
    B = zeros(ind)
    for i in 1:ind
    B[i] = Phi[i] + S[i] * dt
    end
    B[1] += 2 * dt * D * Phi_surf / dz^2
    B[ind] += 2 * dt * D * Phi_surf / dz^2
    return B
end

function computeB_ions(Phi, Phi_surf, D, dz, dt, S, ind)
    B = zeros(ind)
    for i in 1:ind
    B[i] = Phi[i] + S[i] * dt
    # B[i] = Phi[i]
    end
    B[ind] += 2 * dt * D * Phi_surf / dz^2
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

function make_B(NO3, H2PO4)
    NO3_clean = removezeros(NO3)
    H2PO4_clean = removezeros(H2PO4)
    B = []
    for (N, P) in zip(NO3_clean, H2PO4_clean)
        b = -(N + P)
        push!(B, b)
    end
    return B
end

function make_C(H2PO4, CO2, KI, KII, QI, QII, Kw)
    H2PO4_clean = removezeros(H2PO4)
    CO2_clean = removezeros(CO2)
    Carr = []
    for (P, C) in zip(H2PO4_clean, CO2_clean)
        c = -(KI*C + 2*KI*KII*C + 2*QI*P + 3*QI*QII*P + Kw)
        push!(Carr, c)
    end
    return Carr
end