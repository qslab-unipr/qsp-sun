import os
import sys
import math
import time
import numpy as np

from qsp_package import __doc__, qspCircuit, qspParameters, toKet, stateToVector
print(__doc__)

np.set_printoptions(precision=3, suppress=True, linewidth=100)

result_path = 'results/'

if not os.path.exists(result_path):
    os.makedirs(result_path)

try:
    match input("Specific quantum state [s] or random vector [r]? > ").strip().lower():
        case 's':
            match int(input("Select a known quantum state or insert a custom one:\n" \
            "1) Bell (n=2)\n" \
            "2) GHZ (n=3)\n" \
            "3) GHZ (n=4)\n" \
            "4) W (n=3)\n" \
            "5) W (n=4)\n" \
            "6) Dicke (n=3)\n" \
            "7) Dicke (n=4)\n" \
            "8) Custom\n" \
            " > ").strip()):
                case 1:
                    N = 2
                    coefficients = [math.sqrt(1/2), 0, 0, math.sqrt(1/2)]
                case 2:
                    N = 3
                    coefficients = [math.sqrt(1/2), 0, 0, 0, 0, 0, 0, math.sqrt(1/2)]
                case 3:
                    N = 4
                    coefficients = [math.sqrt(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, math.sqrt(1/2)]
                case 4:
                    N = 3
                    coefficients = [0, math.sqrt(1/3), math.sqrt(1/3), 0, math.sqrt(1/3), 0, 0, 0]
                case 5:
                    N = 4
                    coefficients = [0, math.sqrt(1/4), math.sqrt(1/4), 0, math.sqrt(1/4), 0, 0, 0, math.sqrt(1/4), 0, 0, 0, 0, 0, 0, 0]
                case 6:
                    N = 3
                    coefficients = [0, 0, 0, math.sqrt(1/3), 0, math.sqrt(1/3), math.sqrt(1/3), 0]
                case 7:
                    N = 4
                    coefficients = [0, 0, 0, math.sqrt(1/6), 0, math.sqrt(1/6), math.sqrt(1/6), 0, 0, math.sqrt(1/6), math.sqrt(1/6), 0, math.sqrt(1/6), 0, 0, 0]
                case 8:
                    N = int(input("Number of qubits? > ").strip())
                    coefficients = list(map(float, input("Insert real coefficients (floats) separated by spaces > ").strip().split()))
                case _:
                    raise Exception("Invalid input!")

        case 'r':
            N = int(input("Number of qubits? > ").strip())

            match input("Dense [d] or sparse [s] random state? > ").strip().lower():
                case 'd':
                    match input("Real [r] or complex [c] coefficients? > ").strip().lower():
                        case 'r':
                            real_numbers = np.random.normal(size=2**N)
                            real_numbers /= np.linalg.norm(real_numbers, ord=2)
                            coefficients = real_numbers.tolist()
                        case 'c':
                            complex_numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
                            complex_numbers /= np.linalg.norm(complex_numbers, ord=2)
                            coefficients = complex_numbers.tolist()
                        case _:
                            raise Exception("Invalid input!")
                case 's':
                    match input("Real [r] or complex [c] coefficients? > ").strip().lower():
                        case 'r':
                            num_zeros = np.random.randint(1, 2**N)
                            zero_indices = np.random.choice(2**N, num_zeros, replace=False)
                            real_numbers = np.random.normal(size=2**N)
                            real_numbers[zero_indices] = 0
                            real_numbers /= np.linalg.norm(real_numbers, ord=2)
                            coefficients = real_numbers.tolist()
                        case 'c':
                            complex_numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
                            num_zeros = np.random.randint(1, 2**N)
                            zero_indices = np.random.choice(2**N, num_zeros, replace=False)
                            complex_numbers[zero_indices] = 0
                            complex_numbers /= np.linalg.norm(complex_numbers, ord=2)
                            coefficients = complex_numbers.tolist()
                        case _:
                            raise Exception("Invalid input!")
                case _:
                    raise Exception("Invalid input!")
        case _:
            raise Exception("Invalid input!")
        
except Exception as e:
    print(e)
    sys.exit(1)

sys.stdout = open(result_path + 'qsp_output.txt', 'w', encoding='utf-8') # Comment this to write on console (comment also stdout.close() at the end)

print(f"Coefficients: {coefficients}")
print(f"Modulo: {[np.sqrt(np.real(x)**2 + np.imag(x)**2) for x in coefficients]}")
print("_"*50)

M = 2 * N

param_start_time = time.perf_counter() # START PARAMETERS GENERATION

# LAMBDA PARAMETERS
prefix_alphas, middle_alphas, suffix_alphas, ucg_0_angles, global_phases, lambda_diagonals = qspParameters(coefficients, N)
alphas = [prefix_alphas, middle_alphas, suffix_alphas]

param_end_time = time.perf_counter() # END PARAMETERS GENERATION

# QSP CIRCUIT
circuit = qspCircuit(N, M, alphas, ucg_0_angles)

circ_start_time = time.perf_counter() # START QUANTUM COMPUTATION

# OUTPUT STATE
output_state = circuit.computeQuantumState(modulo=False)

circ_end_time = time.perf_counter() # END QUANTUM COMPUTATION

print(f"\nQSP quantum state with global phases:\n{toKet([math.prod(global_phases) * x for x in output_state], N+M)}")

# METRICS
new_qubits = [0] * M
state_vector = stateToVector(new_qubits)
coefficients = np.kron(coefficients, state_vector)
circuit.printMetrics(coefficients, output_state)

print(f"\nTime taken to generate parameters: {param_end_time - param_start_time}s.")
print(f"Time taken to execute the circuit: {circ_end_time - circ_start_time}s.")

# CIRCUIT INFOS
#specs = circuit.getCircuitInfos()
#print(f"\n{specs['resources']}")

sys.stdout.close()
sys.__stdout__.write(f"\nTask completed! Results saved in {result_path} directory.")