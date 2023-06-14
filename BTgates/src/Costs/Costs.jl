using LinearAlgebra
include("../Amplitudes/Chebyshev.jl")

function cost(unitary, propagator)
    # computes the value of the cost function given the unitary and the propagator
    c=tr(propagator*unitary')/size(unitary,1)
    1-norm(c)
end

function cost_from_0_dispersive(H_drift, H_ctrl_c, H_ctrl_q, T, δt, coefficients, unitary)

    # initialising the propagator
    dim = size(H_drift,1)
    propagator = Matrix{ComplexF64}(I,dim,dim)
    amplitude_c(t) = chebyshev_amplitude(coefficients[1:Int(length(coefficients)/2)], T, t)
    amplitude_q(t) = chebyshev_amplitude(coefficients[Int(length(coefficients)/2) + 1:end], T, t)

    # time ordered product of the single exponential matrices
    for l in 0:δt:T
        H = H_drift + amplitude_c(l) * H_ctrl_c(l) + amplitude_q(l) * H_ctrl_q(l)
        infinitesimal_propagator  = exp(-1im * H * δt)
        propagator = infinitesimal_propagator * propagator
    end

    c=tr(propagator*unitary')/size(unitary,1)
    1-norm(c)
    return 1-norm(c)
    
end


function cost_from_0_stepwise(H_drift, H_ctrl_c, H_ctrl_q, coefficients, unitary)

    # initialising the propagator
    dim = size(H_drift,1)
    propagator = Matrix{ComplexF64}(I,dim,dim)
    len = Int(length(coefficients)/2)
    # time ordered product of the single exponential matrices
    for l in 1:1:len
        H = H_drift + coefficients[l] * H_ctrl_c(l / len) + coefficients[len + l] * H_ctrl_q(l / len)
        infinitesimal_propagator  = exp(-1im * H / len)
        propagator = infinitesimal_propagator * propagator
    end

    c=tr(propagator*unitary')/size(unitary,1)
    1-norm(c)
    return 1-norm(c)

end 

function cost_dispersive_norm(H_drift, H_ctrl_c, H_ctrl_q, T, δt, coefficients, unitary)

    # initialising the propagator
    dim = size(H_drift,1)
    propagator = Matrix{ComplexF64}(I,dim,dim)
    amplitude_c(t) = chebyshev_amplitude(coefficients[1:Int(length(coefficients)/2)], T, t)
    amplitude_q(t) = chebyshev_amplitude(coefficients[Int(length(coefficients)/2) + 1:end], T, t)

    # time ordered product of the single exponential matrices
    for l in 0:δt:T
        H = H_drift + amplitude_c(l) * H_ctrl_c(l) + amplitude_q(l) * H_ctrl_q(l)
        infinitesimal_propagator  = exp(-1im * H * δt)
        propagator = infinitesimal_propagator * propagator
    end

    c = norm(propagator*unitary' - Matrix{ComplexF64}(I,dim,dim))
    return c
    
end