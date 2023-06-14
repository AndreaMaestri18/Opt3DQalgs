include("../Amplitudes/Chebyshev.jl")

function propagator_dispersive(H_drift, H_ctrl_c, H_ctrl_q, T, δt, coefficients)

    # initialising the propagator
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

    return propagator
end