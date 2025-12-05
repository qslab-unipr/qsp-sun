# The optimized version of the Sun et al.'s QSP algorithm is based on a new decomposition
# which prepares the real part of the desired state separately from the complex one.
# For state with real positive coefficients there is an even more reduced version.
# This module contains the definition of the optimized QSP circuit
# for the general complex case and for its reductions.

from circuit_classes import *

def qspCircuitReduced(N : int, M : int, alphas : list, ucg_0_angle : float) -> 'QuantumCircuit':

    # QUANTUM STATE PREPARATION CIRCUIT FOR REAL CASE
    circuit = QuantumCircuit(N + M)

    for n in range(1, N+1): 

        match n:
            case 1:

                ry = RotY(N, 1, ucg_0_angle, 0)
                circuit.addStage(ry)

            case _:

                basis_vectors = compBasis(n)[1:]
                alpha_vector = {basis_vectors[idx]: coeff for idx, coeff in enumerate(alphas[n-2])}
                    
                # UNITATRY OPERATORS (S+, H)
                pre = UnitPre(n, n, n-1)
                circuit.addStage(pre)

                # LAMBDA STAGES
                s1 = PrefixCopyStage(n, n)
                s2 = GrayInitialStage(n, n, alpha_vector)
                s3 = SuffixCopyStage(n, n)
                s4 = GrayPathStage(n, n, alpha_vector)
                s5 = InverseStage(n, n)

                circuit.addStage(s1)
                circuit.addStage(s2)
                circuit.addStage(s3)
                circuit.addStage(s4)
                circuit.addStage(s5)

                # UNITATRY OPERATORS (H, S)
                post = UnitPost(n, n, n-1)
                circuit.addStage(post)

    return circuit


def qspCircuitReduced_mod(N : int, M : int, alphas_mag : list, alphas_pha : list, ucg_0_angle : float) -> 'QuantumCircuit':

    # QUANTUM STATE PREPARATION CIRCUIT WITH LAMBDA2 CONSTRUCTORS
    circuit = QuantumCircuit(N + M)

    for n in range(1, N+1): 

        match n:
            case 1:

                ry = RotY(N, 1, ucg_0_angle, 0)
                circuit.addStage(ry)

            case _:

                basis_vectors = compBasis(n)[1:]
                alpha_vector1 = {basis_vectors[idx]: coeff for idx, coeff in enumerate(alphas_mag[n-2])}
                    
                # UNITARY OPERATORS (S+, H)
                pre = UnitPre(n, n, n-1)
                circuit.addStage(pre)

                # LAMBDA STAGES
                s1 = PrefixCopyStage(n, n)
                s2 = GrayInitialStage(n, n, alpha_vector1)
                s3 = SuffixCopyStage(n, n)
                s4 = GrayPathStage(n, n, alpha_vector1)
                s5 = InverseStage(n, n)

                circuit.addStage(s1)
                circuit.addStage(s2)
                circuit.addStage(s3)
                circuit.addStage(s4)
                circuit.addStage(s5)

                # UNITARY OPERATORS (H, S)
                post = UnitPost(n, n, n-1)
                circuit.addStage(post)
    
    basis_vectors = compBasis(N)[1:]
    alpha_vector2 = {basis_vectors[idx]: coeff for idx, coeff in enumerate(alphas_pha)}

    # FINAL LAMBDA FOR COMPLEX PHASES
    s1 = PrefixCopyStage(N, N)
    s2 = GrayInitialStage(N, N, alpha_vector2)
    s3 = SuffixCopyStage(N, N)
    s4 = GrayPathStage(N, N, alpha_vector2)
    s5 = InverseStage(N, N)

    circuit.addStage(s1)
    circuit.addStage(s2)
    circuit.addStage(s3)
    circuit.addStage(s4)
    circuit.addStage(s5)

    return circuit