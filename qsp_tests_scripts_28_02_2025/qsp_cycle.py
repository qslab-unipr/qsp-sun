from circuit_classes_cycle import *
import sys

sys.stdout = open('n2_sparse_complex_30', 'w') # Comment this to write on console (comment also stdout.close() at the end)

np.set_printoptions(precision=3, suppress=True, linewidth=100)

# NUMER OF QUBITS
N = 2
M = 2 * N

# INPUT STATES
initial_state = np.array([0]*(N+M))
initial_state_vector = stateToVector(initial_state)

for _ in range(30):
    # DENSE COEFFICIENTS
    complex_numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
    complex_numbers /= np.linalg.norm(complex_numbers, ord=2)
    coefficients = complex_numbers.tolist() # Modify this to insert the desired vector

    # LAMBDA PARAMETERS
    prefix_alphas, middle_alphas, suffix_alphas, ucg_0_angles, global_phases, lambda_diagonals = qspParameters(coefficients, N)
    alphas = [prefix_alphas, middle_alphas, suffix_alphas]

    # TENSOR PRODUCT WITH ANCILLARIES
    new_qubits = [0] * M
    state_vector = stateToVector(new_qubits)
    coefficients = np.kron(coefficients, state_vector)

    # QUANTUM STATE PREPARATION CIRCUIT
    circuit = QuantumCircuit(N + M, initial_state_vector)

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

    output_state = circuit.computeQuantumState(modulo=False)

    circuit.printMetrics(coefficients, output_state)
    print("-"*50)

sys.stdout.close()