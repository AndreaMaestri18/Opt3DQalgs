using LinearAlgebra

function annihilation_operator(N::Int)
    A = zeros(Float64, N, N)
    for i in 1:N-1
        A[i, i+1] = sqrt(i)
    end
    return A
end

function generate_matrices(fock_space_dimension)
    """
    generates the matrices for the optimisation problem

    INPUT:
    dimensions of the fock space as input

    RETURNS:
    annihilation and creation operator as a tensor product of the 2x2 identity and the fock_space_dimension x fock_space_dimension a and a_dagger
    σp σz σm as a tensor product of the 2x2 σp σz σm and the fock_space_dimension x fock_space_dimension Identity matrix
    
    """
    annihilation = annihilation_operator(fock_space_dimension)
    sp_initial = [0 1; 0 0]
    sm_initial = [0 0; 1 0]
    sz_initial = [1 0; 0 -1]
    I_qubit = Matrix(I,2,2)
    I_photon = Matrix(I,fock_space_dimension,fock_space_dimension)
    a = kron(annihilation,I_qubit)
    adag = kron(annihilation',I_qubit)
    sp = kron(I_photon,sp_initial)
    sm = kron(I_photon,sm_initial)
    sz = kron(I_photon,sz_initial)
    
    return a,adag,sp,sm,sz
end
