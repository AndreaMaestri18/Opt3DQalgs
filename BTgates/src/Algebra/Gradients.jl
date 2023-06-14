using LinearAlgebra
include("../Amplitudes/Chebyshev.jl")

function gradient(H_drift, H_ctrl_c, H_ctrl_q, T, δt, coefficients, unitary)

    h = 10e-7
    g=[]
    dim = size(H_drift,1)
    propagator = Matrix{ComplexF64}(I,dim,dim)
    amplitude_c(t) = chebyshev_amplitude(coefficients[1:Int(length(coefficients)/2)], T, t)
    amplitude_q(t) = chebyshev_amplitude(coefficients[Int(length(coefficients)/2) + 1:length(coefficients)], T, t)

    # time ordered product of the single exponential matrices
    for l in 0:δt:T
        H = H_drift + amplitude_c(l) * H_ctrl_c(l) + amplitude_q(l) * H_ctrl_q(l)
        infinitesimal_propagator  = exp(-1im * H * δt)
        propagator = infinitesimal_propagator * propagator
    end

    c = 1 - norm(tr(propagator*unitary')/dim)
    
    for j=1:1:length(coefficients)
        temp = deepcopy(coefficients)
        temp[j] = coefficients[j] + h
        amplitude_c_temp(t) = chebyshev_amplitude(temp[1:Int(length(temp)/2)], T, t)
        amplitude_q_temp(t) = chebyshev_amplitude(temp[Int(length(temp)/2) + 1:length(temp)], T, t)

        # initialising the propagator at x + h
        propagator_plus_dx = Matrix{ComplexF64}(I,dim,dim)
        # time ordered product of the single exponential matrices
        for l in 0:δt:T
            H = H_drift + amplitude_c_temp(l) * H_ctrl_c(l) + amplitude_q_temp(l) * H_ctrl_q(l)
            infinitesimal_propagator  = exp(-1im * H * δt)
            propagator_plus_dx = infinitesimal_propagator * propagator_plus_dx
        end
        c_plus_dx = 1 - norm(tr(propagator_plus_dx*unitary')/dim)
        push!(g,(c_plus_dx-c)/h)
    end

    return g, c
end
