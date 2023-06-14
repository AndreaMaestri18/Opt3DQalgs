using Polynomials

function chebyshev_amplitude(coefficients, T, t)
    return ChebyshevT(coefficients)((2t - T)/T)
end

