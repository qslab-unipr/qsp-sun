# This module computes and analyzes the parameters for the last $\Lambda$-type constructor,
# the one dedicated to the preparation of complex phases. 

import numpy as np
import sys
from utils_mod import cartesian_to_polar, generateInnProdMatrix, generateAlphaVect

def matrixLambdaFinal(state_vector):

    # Step 1: Print the desired quantum state
    print("\n[Step 1] Input quantum state:\n")
    print(state_vector)
    assert np.isclose(np.linalg.norm(state_vector), 1, atol = 1e-6), f"Coefficient vector norm is not 1: {np.linalg.norm(state_vector)} != 1."

    # Step 2: Cartesian coordinates (real and imaginary parts)
    real_parts = np.real(state_vector)
    imag_parts = np.imag(state_vector)

    print("\n[Step 2] Cartesian coordinates for each component:\n")
    print("Real part: ", real_parts)
    print("Immaginary part: ", imag_parts)
 
    # Step 3: Cartesian to polar
    magnitudes, angles = cartesian_to_polar(state_vector)
    print("\n[Step 3] Polar coordinates for each component:\n")
    print("Magnitudes: ", magnitudes)
    print("Phase angles [rad]: ", angles)

    # Step 4: Magnitudes and complex phases
    magnitudes = np.abs(state_vector)
    phases = np.array([z / abs(z) if abs(z) > 1e-10 else 1.0 for z in state_vector])
    print("\n[Step 4] Euler representation for each component:\n")
    print("Magnitudes: ", magnitudes)
    print("Phases: ", phases)

    # Step 5: Collect the first entry in the array "phases"
    phase0 = phases[0]
    new_phases = phases / phase0
    new_phase_angles = np.angle(new_phases)
    print("\n[Step 5] New complex phases and related arguments:\n")
    print("New complex phases: ", new_phases)
    print("New phase angles [rad]:", new_phase_angles)
    print("Global phase: ", phase0)

    # Step 6: Diagonal matrix with array "new_phases"
    diagonal_matrix = np.diag(new_phases)
    print("\n[Step 6] Diagonal matrix from new complex phases:\n")
    print(diagonal_matrix)

    return {
        "Magnitudes": magnitudes,
        "New complex phases": new_phases,
        "New phase angles": new_phase_angles,
        "GLobal phase": phase0,
        "Diagonal matrix": diagonal_matrix}

# === TEST ===
N = int(input("Number of qubits? > ").strip())
try:
    match input("Specific quantum state [s] or random vector [r]? > ").strip().lower():
        case 's':
            #coefficients = [math.sqrt(1/2)+0j, 0+1j, 0-1j, -math.sqrt(1/2)+0j]
            coefficients = [0.+0.5j, 0.5+0.j, 0.-0.5j, -0.5+0.j]
            coefficients /= np.linalg.norm(coefficients, ord=2)
            assert np.isclose(np.linalg.norm(coefficients), 1, atol = 1e-6), f"Norm is not 1: {np.linalg.norm(coefficients)} != 1."
        case 'r':
            random_coeff = np.random.normal(size=2**N) + 1j *  np.random.normal(size=2**N)
            random_coeff /= np.linalg.norm(random_coeff, ord=2)
            assert np.isclose(np.linalg.norm(random_coeff), 1, atol = 1e-6), f"Norm is not 1: {np.linalg.norm(random_coeff)} != 1."
            coefficients = random_coeff.tolist()
        case _:
            raise Exception("Invalid input!")
except Exception as e:
    print(e)
    sys.exit(1)

risultati = matrixLambdaFinal(coefficients)


# Step 7: Solve the linear system for parallelized phases
phases_alphas = []

comb = risultati["New phase angles"]
print("\n[Step 7a] Checking the new phase angles that enter the linear system:\n")
print(comb)

binary_matrix = generateInnProdMatrix(N)
print("\n[Step 7b] Checking the coefficient matrix:\n")
print(binary_matrix)
phases_alphas = generateAlphaVect(binary_matrix, comb[1:])
alpha_check = (pow(2, 1 - N) * (2 * binary_matrix - np.ones((pow(2, N) - 1,pow(2, N) - 1))) @ comb[1:])
assert np.allclose(alpha_check, phases_alphas, atol = 1e-6), "Alpha values don't match the check: \n" + str(phases_alphas) + " -> " + str(alpha_check)
print("\n[Step 7c] Checking the alphas:\n")
print(phases_alphas)
