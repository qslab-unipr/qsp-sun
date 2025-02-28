from circuit_classes import *
import sys

sys.stdout = open('results/console_output', 'w') # Comment this to write on console (comment also stdout.close() at the end)

np.set_printoptions(precision=3, suppress=True, linewidth=100)

# NUMER OF QUBITS
N = 4
M = 2 * N

# INPUT STATES
initial_state = np.array([0]*(N+M))
initial_state_vector = stateToVector(initial_state)

# NEW COEFFICIENTS GENERATION PART
complex_numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
complex_numbers /= np.linalg.norm(complex_numbers, ord=2)
coefficients = complex_numbers.tolist() # Modify this to insert the desired vector

print(f"Coefficients: {coefficients}")
print("_"*50)

# LAMBDA PARAMETERS
prefix_alphas, middle_alphas, suffix_alphas, ucg_0_angles, global_phases, lambda_diagonals = qspParameters(coefficients, N)
alphas = [prefix_alphas, middle_alphas, suffix_alphas]

lambda_check = np.diag(np.array(lambda_diagonals[1]))

# TENSOR PRODUCT WITH ANCILLARIES
new_qubits = [0] * M
state_vector = stateToVector(new_qubits)
coefficients = np.kron(coefficients, state_vector)

# QUANTUM STATE PREPARATION CIRCUIT
circuit = QuantumCircuit(N + M, initial_state_vector)

for m, a in enumerate(reversed(alphas)):    
    if m == 1:
        basis_vectors = compBasis(N)[1:]
        alpha_vector = {basis_vectors[idx]: coeff for idx, coeff in enumerate(a[N-2])}

        # LAMBDA STAGES
        s1 = PrefixCopyStage(N, N)
        s2 = GrayInitialStage(N, N, alpha_vector)
        s3 = SuffixCopyStage(N, N)
        s4 = GrayPathStage(N, N, alpha_vector)
        s5 = InverseStage(N, N)

        circuit.addStage(s1)
        circuit.addStage(s2)
        circuit.addStage(s3)
        circuit.addStage(s4)
        circuit.addStage(s5)

random_state_vector = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
random_state_vector /= np.linalg.norm(random_state_vector, ord=2)

new_qubits_sv = stateToVector(new_qubits)

th_state = lambda_check @ random_state_vector
th_state = np.kron(th_state, new_qubits_sv)

ac_state = circuit.computeQuantumState(modulo=False, state_vector=np.kron(random_state_vector, new_qubits_sv))

circuit.printMetrics(th_state, ac_state)

#circuit.printMatrixComparison(lambda_check, circuitMatrixReduction(circuit.computeMatrix(modulo=False), N, M)) # Takes very long time for N >= 4

#print(circuit.getInfos())

sys.stdout.close()