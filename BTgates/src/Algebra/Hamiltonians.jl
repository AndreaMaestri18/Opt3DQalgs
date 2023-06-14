using LinearAlgebra

include("Matrices.jl")

function Jaynes_Cummings_H(fock_space_dimension, ω_c, ω_q, g)

    # generating matrices
    a,adag,sp,sm,sz = generate_matrices(fock_space_dimension)

    # jaynes cummings hamiltonian in non-dispersive regime
    return ω_c * adag * a + g * (sp * a + sm * adag) + ω_q * sz
end

function Jaynes_Cummings_Dispersive(fock_space_dimension, ω_c, ω_q, χ)

    # generating matrices
    a,adag,sp,sm,sz = generate_matrices(fock_space_dimension)

    I_qubit = Matrix(I,2,2)
    I_photon = Matrix(I,fock_space_dimension,fock_space_dimension)
    I_big = kron(I_photon, I_qubit)

    # jaynes cummings hamiltonian in non-dispersive regime
    return ω_c * (adag * a + 1/2 * I_big) + χ * (adag * a + 1/2 * I_big) * sz + ( ω_q / 2 ) * sz 
end

function control_hamiltonian(fock_space_dimension, ω_control, t, amp)

    # generating matrices
    a,adag,sp,sm,sz = generate_matrices(fock_space_dimension)

    return cos(ω_control * t) * (sp + sm)
end

function control_hamiltonian_cavity(fock_space_dimension, ω_c, t, amp)

    # generating matrices
    a,adag,sp,sm,sz = generate_matrices(fock_space_dimension)

    return amp * (exp(1im * ω_c * t ) * adag + exp(-1im * ω_c * t ) * a)
    
end

function control_hamiltonian_qubit(fock_space_dimension, ω_q, t, amp)

    # generating matrices
    a,adag,sp,sm,sz = generate_matrices(fock_space_dimension)

    return amp * (exp(1im * ω_q * t ) * sp + exp(-1im * ω_q * t ) * sm)
end


