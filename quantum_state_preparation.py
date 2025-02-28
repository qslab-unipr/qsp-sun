import math
from circuit_classes import *
import sys

sys.stdout = open('results/console_output', 'w') # Comment this to write on console (comment also stdout.close() at the end)

np.set_printoptions(precision=3, suppress=True, linewidth=100)

# NUMER OF QUBITS
N = 2
M = 2 * N

# INPUT STATES
initial_state = np.array([0]*(N+M))
initial_state_vector = stateToVector(initial_state)

# DENSE DOMPLEX COEFFICIENTS
complex_numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
complex_numbers /= np.linalg.norm(complex_numbers, ord=2)
coefficients = complex_numbers.tolist() # Modify this to insert the desired vector

# DENSE REAL COEFFICIENTS
#real_numbers = np.random.normal(size=2**N)
#real_numbers /= np.linalg.norm(real_numbers, ord=2)
#coefficients = real_numbers.tolist() # Modify this to insert the desired vector

# SPARSE COMPLEX COEFFICIENTS
#complex_numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
#num_zeros = np.random.randint(1, 2**N)
#zero_indices = np.random.choice(2**N, num_zeros, replace=False)
#complex_numbers[zero_indices] = 0
#complex_numbers /= np.linalg.norm(complex_numbers, ord=2)
#coefficients = complex_numbers.tolist() # Modify this to insert the desired vector

# SPARSE REAL COEFFICIENTS
#num_zeros = np.random.randint(1, 2**N)
#zero_indices = np.random.choice(2**N, num_zeros, replace=False)
#real_numbers = np.random.normal(size=2**N)
#real_numbers[zero_indices] = 0
#real_numbers /= np.linalg.norm(real_numbers, ord=2)
#coefficients = real_numbers.tolist() # Modify this to insert the desired vector

print(f"Coefficients: {coefficients}")
print(f"Modulo: {[np.sqrt(np.real(x)**2 + np.imag(x)**2) for x in coefficients]}")
print("_"*50)

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


#print(f"\nGlobal phases: {global_phases}")

output_state = circuit.computeQuantumState(modulo=False)

print(f"\nQSP quantum state with global phases:\n{toKet([math.prod(global_phases) * x for x in output_state], N+M)}")

#print(f"\nOutput state (MODULO):")
#circuit.printCircuit(mode="ket", modulo=True)

#print(f"\nOutput state (without global phase):")
#circuit.printCircuit(mode="ket", modulo=False)

#circuit.printCircuit(mode="matrix", modulo=False)
#circuit.printCircuit(mode="matrix", modulo=True)
#circuit.printCircuit(mode="figure", modulo=False)
#circuit.printCircuit(mode="console", modulo=False)

circuit.printMetrics(coefficients, output_state)

#print(circuit.getCircuitInfos())

sys.stdout.close()