using LinearAlgebra

function x_gate(fock_space_dimension, indeces)
    """
    generate a generalized fock_space_dimension dimensional x_gate, swapping the states |indeces[1]> and |indeces[2]>
    """
    x_gate = Matrix{Float64}(I, fock_space_dimension, fock_space_dimension)
    for index in indeces
        x_gate[index[1],index[1]] = x_gate[index[2],index[2]] = 0
        x_gate[index[1],index[2]] = x_gate[index[2],index[1]] = 1
    end
    
    return x_gate
    
end

function Rz(fock_space_dimension, indeces, θ)
    """
    generate a generalized fock_space_dimension dimensional x_gate, swapping the states |indeces[1]> and |indeces[2]>
    """
    Rz = Matrix{ComplexF64}(I, fock_space_dimension, fock_space_dimension)
    for index in indeces
        Rz[index[1],index[1]] = exp(-im*θ/2)
        Rz[index[2],index[2]] = exp(-im*θ/2)
    end

    
    return Rz
    
end
function Rx(fock_space_dimension, indeces, θ)
    """
    generate a generalized fock_space_dimension dimensional x_gate, swapping the states |indeces[1]> and |indeces[2]>
    """
    Rx = Matrix{ComplexF64}(I, fock_space_dimension, fock_space_dimension)
    for index in indeces
        Rx[index[1],index[1]] = Rx[index[2],index[2]] = cos(θ/2)
        Rx[index[1],index[2]] = Rx[index[2],index[1]] = -im*sin(θ/2)
    end
    
    return Rx
    
end
function U2(fock_space_dimension, indeces, θs)
    """
    generate a generalized fock_space_dimension dimensional x_gate, swapping the states |indeces[1]> and |indeces[2]>
    """
    U = Rz(fock_space_dimension, indeces, θs[1]) *  Rx(fock_space_dimension, indeces, θs[2]) *  Rz(fock_space_dimension, indeces, θs[3])
    
    
    
    return U
    
end

