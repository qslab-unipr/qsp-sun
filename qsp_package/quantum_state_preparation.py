from .circuit_classes import *

def qspCircuit(N : int, M : int, alphas : list, ucg_0_angles : list) -> 'QuantumCircuit':

    # QUANTUM STATE PREPARATION CIRCUIT
    circuit = QuantumCircuit(N + M)

    for n in range(1, N+1): 

        match n:
            case 1:
                
                rz1 = RotZ(N, 1, ucg_0_angles[2], 0)
                circuit.addStage(rz1)

                ry = RotY(N, 1, ucg_0_angles[1], 0)
                circuit.addStage(ry)

                rz2 = RotZ(N, 1, ucg_0_angles[0], 0)
                circuit.addStage(rz2)

            case _:

                for m, a in enumerate(reversed(alphas)):
                    basis_vectors = compBasis(n)[1:]
                    alpha_vector = {basis_vectors[idx]: coeff for idx, coeff in enumerate(a[n-2])}
                    
                    if m == 1:
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
                    
                    if m == 1:
                        # UNITATRY OPERATORS (H, S)
                        post = UnitPost(n, n, n-1)
                        circuit.addStage(post)

    return circuit

