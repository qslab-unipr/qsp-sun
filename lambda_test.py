import os
import sys
import numpy as np

from qsp_package import __doc__, lambdaCircuit, lambdaTest, qspParameters
print(__doc__)

np.set_printoptions(precision=3, suppress=True, linewidth=100)

result_path = 'results/'

if not os.path.exists(result_path):
    os.makedirs(result_path)

N = int(input("Number of qubits of the input register? (for n >= 4 fast hardware and large memory are needed) > ").strip())
M = 2 * N

sys.stdout = open(result_path + 'lambda_test.txt', 'w', encoding='utf-8')

# RANDOM DENSE COMPLEX COEFFICIENTS
complex_numbers = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
complex_numbers /= np.linalg.norm(complex_numbers, ord=2)
coefficients = complex_numbers.tolist()

# LAMBDA PARAMETERS
prefix_alphas, middle_alphas, suffix_alphas, ucg_0_angles, global_phases, lambda_diagonals = qspParameters(coefficients, N)
alphas = [prefix_alphas, middle_alphas, suffix_alphas]

# QSP CIRCUIT
circuit = lambdaCircuit(N, M, alphas)

lambdaTest(N, M, circuit, lambda_diagonals)

# CIRCUIT INFOS
specs = circuit.getCircuitInfos()
print(f"\n{specs['resources']}")

sys.stdout.close()
sys.__stdout__.write(f"\nTask completed! Results saved in {result_path} directory.")