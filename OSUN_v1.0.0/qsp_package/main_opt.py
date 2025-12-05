import os
import sys
import math
import time
import numpy as np
from quantum_state_preparation_mod  import qspCircuitReduced, qspCircuitReduced_mod
from utils_mod import toKet, stateToVector, qspParameters_mod

np.set_printoptions(precision=3, suppress=True, linewidth=100)

result_path = 'results/'

if not os.path.exists(result_path):
    os.makedirs(result_path)

try:
    match input("Specific quantum state [s] or random vector [r]? > ").strip().lower():
        case 's':
            match int(input("Select a known quantum state or insert a custom one:\n" \
            "1) Bell Phi+ (n=2)\n" \
            "2) Bell Phi- (n=2)\n" \
            "3) Bell Psi+ (n=2)\n" \
            "4) Bell Psi- (n=2)\n" \
            "5) GHZ (n=3)\n" \
            "6) GHZ (n=4)\n" \
            "7) W (n=3)\n" \
            "8) W (n=4)\n" \
            "9) Dicke (n=3)\n" \
            "10) Dicke (n=4)\n" \
            "11) Custom\n" \
            " > ").strip()):
                case 1:
                    N = 2
                    coefficients = [math.sqrt(1/2), 0, 0, math.sqrt(1/2)]
                case 2:
                    N = 2
                    coefficients = [math.sqrt(1/2), 0, 0, -math.sqrt(1/2)]
                case 3:
                    N = 2
                    coefficients = [0, math.sqrt(1/2), math.sqrt(1/2), 0]
                case 4:
                    N = 2
                    coefficients = [0, math.sqrt(1/2), -math.sqrt(1/2), 0]
                case 5:
                    N = 3
                    coefficients = [math.sqrt(1/2), 0, 0, 0, 0, 0, 0, math.sqrt(1/2)]
                case 6:
                    N = 4
                    coefficients = [math.sqrt(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, math.sqrt(1/2)]
                case 7:
                    N = 3
                    coefficients = [0, math.sqrt(1/3), math.sqrt(1/3), 0, math.sqrt(1/3), 0, 0, 0]
                case 8:
                    N = 4
                    coefficients = [0, math.sqrt(1/4), math.sqrt(1/4), 0, math.sqrt(1/4), 0, 0, 0, math.sqrt(1/4), 0, 0, 0, 0, 0, 0, 0]
                case 9:
                    N = 3
                    coefficients = [0, 0, 0, math.sqrt(1/3), 0, math.sqrt(1/3), math.sqrt(1/3), 0]
                case 10:
                    N = 4
                    coefficients = [0, 0, 0, math.sqrt(1/6), 0, math.sqrt(1/6), math.sqrt(1/6), 0, 0, math.sqrt(1/6), math.sqrt(1/6), 0, math.sqrt(1/6), 0, 0, 0]
                case 11:
                    N = int(input("Number of qubits? > ").strip())
                    coefficients = list(map(float, input("Insert real coefficients (floats) separated by spaces > ").strip().split()))
                case _:
                    raise Exception("Invalid input!")

        case 'r':
            N = int(input("Number of qubits? > ").strip())

            match input("Dense [d] or sparse [s] random state? > ").strip().lower():
                case 'd':
                    match input("Real positive [p], real negative [n] or complex [c] coefficients? > ").strip().lower():
                        case 'p':
                            numbers = np.abs(np.random.normal(size=2**N))
                        case 'n':
                            numbers = np.random.normal(size=2**N)
                        case 'c':
                            numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
                        case _:
                            raise Exception("Invalid input!")
                        
                    numbers /= np.linalg.norm(numbers, ord=2)
                    coefficients = numbers.tolist()
                    
                case 's':
                    match input("Real positive [p], real negative [n] or complex [c] coefficients? > ").strip().lower():
                        case 'p':
                            numbers = np.abs(np.random.normal(size=2**N))
                        case 'n':
                            numbers = np.random.normal(size=2**N)
                        case 'c':
                            numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
                        case _:
                            raise Exception("Invalid input!")
                    
                    zero_indices = np.random.choice(2**N, np.random.randint(1, 2**N), replace=False)
                    numbers[zero_indices] = 0
                    numbers /= np.linalg.norm(numbers, ord=2)
                    coefficients = numbers.tolist()

                case _:
                    raise Exception("Invalid input!")
        case _:
            raise Exception("Invalid input!")
        
except Exception as e:
    print(e)
    sys.exit(1)

sys.stdout = open(result_path + 'qsp_output.txt', 'w', encoding='utf-8') # Comment this to write on console (comment also stdout.close() at the end)

print(f"Coefficients: {coefficients}")
#print(f"Modulo: {[np.sqrt(np.real(x)**2 + np.imag(x)**2) for x in coefficients]}")
print(f"Modulo: {[abs(x) for x in coefficients]}")
print("_"*50)

M = 2 * N

param_start_time = time.perf_counter() # START PARAMETERS GENERATION

# LAMBDA_2 + LAMBDA_FINAL PARAMETERS
middle_alphas, phases_alphas, gamma0, global_phases, lambda_diagonals = qspParameters_mod(coefficients, N)

param_end_time = time.perf_counter() # END PARAMETERS GENERATION

# QSP CIRCUIT (REDUCED)
circuit = qspCircuitReduced(N, M, middle_alphas, gamma0) if all([True if x == 0 else False for x in np.angle(coefficients)]) else qspCircuitReduced_mod(N, M, middle_alphas, phases_alphas, gamma0)
#circuit = qspCircuitReduced_mod(N, M, middle_alphas, phases_alphas, gamma0)

circ_start_time = time.perf_counter() # START QUANTUM COMPUTATION

circuit.printCircuit('figure', modulo=False)


# OUTPUT STATE
output_state = circuit.computeQuantumState(modulo=False)

circ_end_time = time.perf_counter() # END QUANTUM COMPUTATION

print(f"\nDesired quantum state with global phases:\n{toKet([math.prod(global_phases) * x for x in output_state], N+M)}")


# METRICS
new_qubits = [0] * M
state_vector = stateToVector(new_qubits)
coefficients = np.kron(coefficients, state_vector)
circuit.printMetrics(coefficients, output_state)

print(f"\nTime taken to generate parameters: {param_end_time - param_start_time}s.")
print(f"Time taken to execute the circuit: {circ_end_time - circ_start_time}s.")

# CIRCUIT INFOS
specs = circuit.getCircuitInfos()
print(f"\n{specs['resources']}")

sys.stdout.close()
sys.__stdout__.write(f"\nTask completed! Results saved in {result_path} directory.")

