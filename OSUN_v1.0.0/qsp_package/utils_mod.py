# This module contains all the functions used to operate with vectors and compute parameters,
# according to the parallelization introduced by Sun et al. (arXiv:2108.06150)

from math import acos
import numpy as np

def binaryInnProd(x : int, y : int, n : int) -> str: # Compute the inner product in F^2 between two states of 'n' qubits
    v = []
    for i in range(n):
        v.append(int(x[i]) * int(y[i]))

    return bin(sum(v))[-1]

def compBasis(n: int) -> str: # Generate the computational basis vectors for 'n' qubits
    strings = []

    for k in range(pow(2, n)):
        strings.append(bin(k)[2:].zfill(n))

    return strings

def cartesian_to_polar(vector):
    magnitudes = np.abs(vector)
    angles = np.angle(vector)
    return magnitudes, angles

def toKet(vector : list, n : int) -> str: # Show the KET notation of an input vector of 'n' qubits
    basis = compBasis(n)
    state = ''

    for i, q in enumerate(vector):
        if q != 0:
            state += (str(q) + '|' + str(basis[i]) + '>  ')

    return state

def grayCode(n_gray : int, n : int):
    up_prefix = np.array([[0]])
    down_prefix = np.array([[1]])
    matrix = np.concatenate((up_prefix, down_prefix))

    match n_gray:
        case 1:
            for _ in range(n - 1):
                up_prefix = np.concatenate((up_prefix, up_prefix))
                down_prefix = np.concatenate((down_prefix, down_prefix))
                matrix = np.concatenate((np.concatenate((matrix, np.flip(matrix, 0))), np.concatenate((up_prefix, down_prefix))), axis = 1)
        case 2:
            for _ in range(n - 1):
                up_prefix = np.concatenate((up_prefix, up_prefix))
                down_prefix = np.concatenate((down_prefix, down_prefix))
                matrix = np.concatenate((np.concatenate((up_prefix, down_prefix)), np.concatenate((matrix, np.flip(matrix, 0)))), axis = 1)
            
    return matrix

def stateToVector(state : np.ndarray):
    state_vector = 1
    for ks in state:
        vector = np.array([1 - ks, ks]) # Vector notation of ket state
        state_vector = np.kron(state_vector, vector)
    return state_vector

def circuitMatrixReduction(diag_matrix, n, m):
    dim_full = 2 ** (n + m)
    selected_indices = [i for i in range(dim_full) if i & ((1 << m) - 1) == 0]
    reduced_diag = np.diag([diag_matrix[i, i] for i in selected_indices])
    return reduced_diag

class Node:
    def __init__(self, value, arc):
        self._value = value # Value of the node
        self._arc_val = arc # Value of the arc from node to parent
        self._children = [] # Children nodes

    def addChild(self, node: 'Node'): # Add a brach to the tree
        self._children.append(node)
    
    def nodeVal(self) -> float: # Value of the node
        return self._value
    
    def arcVal(self) -> str: # Value of the arc
        return self._arc_val
    
    def explore(self, string : str) -> float: # Explore recursively the tree until the desire state
        if not string:
            return self.nodeVal()
        
        token = string[0]
        for childe in self._children:
            if childe.arcVal() == token:
                return childe.explore(string[1:])
    
    def printTree(self, depth : int = 0):
        for childe in self._children:
            childe.printTree(depth + 1)
        ind = '-' * depth
        print(ind + str(self.nodeVal()))

def generateBinTree(vector_v : list, bit : str, k : int): # Generate a binary search tree of depth 'k' from a list of coeffiecients
    value = np.linalg.norm(vector_v, ord = 2) # Norm-2
    node = Node(value, bit)

    if k != 0:
        vector0 = vector_v[:len(vector_v)//2]
        vector1 = vector_v[len(vector_v)//2:]

        node.addChild(generateBinTree(vector0, '0', k - 1))
        node.addChild(generateBinTree(vector1, '1', k - 1))

    return node

def generateThetaVect(unit_vector : np.ndarray, n : int): # Generate the vector of theta from a multiplexor rappresentation in binary strings
    theta = np.array([])

    binary_strings = encodeMultiplexor(np.array(['0']), n)
    binary_strings = sorted(binary_strings, key=len)

    bst = generateBinTree(unit_vector, None, n)

    for child_string in binary_strings:
        
        parent_string = child_string[:-1]
        child = bst.explore(child_string)
        parent = bst.explore(parent_string)

        if parent == 0:
            theta = np.append(theta, acos(0))
        else:
            theta = np.append(theta, acos(child / parent))

    return theta

def generateInnProdMatrix(n: int) -> np.ndarray: # Generate the matrix of inner products between every computational basis vectors except |00...00>
    basis_vectors = compBasis(n)[1:]
    coefficients = np.array([])

    for i in basis_vectors:
        row = np.array([])
        for j in basis_vectors:
            row = np.append(row, float(binaryInnProd(i, j, n)))
        coefficients = np.append(coefficients, row)
    coefficients = np.reshape(coefficients, (len(basis_vectors), len(basis_vectors)))

    return coefficients

def generateAlphaVect(coefficients : np.ndarray, theta_vector : np.ndarray) -> np.ndarray: # Generate the vector of alpha by solving the linear system
    alpha_vector = np.linalg.solve(coefficients, theta_vector)
    return alpha_vector

def encodeMultiplexor(binary_vector : np.ndarray, n : int) -> np.ndarray: # Generate all the binary strings to use for the exploration of the binary tree
    n = n - 1

    if n == 0:
        return binary_vector
    
    elem_0 = '0' + binary_vector[-1]
    branch_0 = np.append(binary_vector, elem_0)
    branch_0 = encodeMultiplexor(branch_0, n)

    elem_1 = '1' + binary_vector[-1]
    branch_1 = np.append(binary_vector, elem_1)
    branch_1 = encodeMultiplexor(branch_1, n)

    binary_vector = np.union1d(branch_0, branch_1)

    return binary_vector

def matrixLambdaFinal(state_vector):

    # Step 1: Print the desired quantum state
    #print("\n[Step 1] Input quantum state:\n")
    #print(state_vector)
    #assert np.isclose(np.linalg.norm(state_vector), 1, atol = 1e-6), f"Coefficient vector norm is not 1: {np.linalg.norm(state_vector)} != 1."

    # Step 2: Cartesian coordinates (real and imaginary parts)
    #real_parts = np.real(state_vector)
    #imag_parts = np.imag(state_vector)

    #print("\n[Step 2] Cartesian coordinates for each component:\n")
    #print("Real part: ", real_parts)
    #print("Immaginary part: ", imag_parts)
 
    # Step 3: Cartesian to polar
    magnitudes, angles = cartesian_to_polar(state_vector)
    #print("\n[Step 3] Polar coordinates for each component:\n")
    #print("Magnitudes: ", magnitudes)
    #print("Phase angles [rad]: ", angles)

    # Step 4: Magnitudes and complex phases
    magnitudes = np.abs(state_vector)
    phases = np.array([z / abs(z) if abs(z) > 1e-10 else 1.0 for z in state_vector])
    #print("\n[Step 4] Euler representation for each component:\n")
    #print("Magnitudes: ", magnitudes)
    #print("Phases: ", phases)

    # Step 5: Collect the first entry in the array "phases"
    phase0 = phases[0]
    new_phases = phases / phase0
    new_phase_angles = np.angle(new_phases)
    #print("\n[Step 5] New complex phases and related arguments:\n")
    #print("New complex phases: ", new_phases)
    #print("New phase angles [rad]:", new_phase_angles)
    #print("Global phase: ", phase0)

    # Step 6: Diagonal matrix with array "new_phases"
    diagonal_matrix = np.diag(new_phases)
    #print("\n[Step 6] Diagonal matrix from new complex phases:\n")
    #print(diagonal_matrix)

    return {
        "Magnitudes": magnitudes,
        "New complex phases": new_phases,
        "New phase angles": new_phase_angles,
        "Global phase": phase0,
        "Diagonal matrix": diagonal_matrix}


def qspParameters_mod(unit_vector : np.ndarray, n : int) -> tuple[list, list, list, list, list]:
    
    # PARAMETERS FOR LAMBDA2

    #assert np.isclose(np.linalg.norm(unit_vector), 1, atol = 1e-6), f"Coefficient vector norm is not 1: {np.linalg.norm(unit_vector)} != 1."

    #assert (len(unit_vector) == pow(2, n)), f"Wrong number of coefficients 'V' or wrong value of 'n': {len(unit_vector)} != {pow(2, n)}."

    #coefficient_matrix, angles, solution, variables, unused_variables = solve_system(unit_vector, n) #commented for opt_qsp_lambda2

    #a = solution[0::3] #commented for opt_qsp_lambda2
    #b = solution[1::3] #commented for opt_qsp_lambda2
    #d = solution[2::3] #commented for opt_qsp_lambda2

    #prefix_vector, suffix_vector = (a, b), d #commented for opt_qsp_lambda2

    theta_vector = generateThetaVect(unit_vector, n)  # [(2^n)-1] theta angles
    
    #print(f"Thetas generated from BST: {theta_vector}")
    #print(f"Prefix lambda angles (alphas, betas): {prefix_vector}")
    #print(f"Suffix lambda angles (deltas): {suffix_vector}")
    #print("_"*50)

    #params_vector = [prefix_vector, theta_vector, suffix_vector] #commented for opt_qsp_lambda2, replaced by the next line
    params_vector = theta_vector

    #prefix_alphas, middle_alphas, suffix_alphas = [], [], [] # Cannot use numpy because the resulting array has different shapes
    #previous line commented for opt_qsp_lambda2, replaced by the next line
    middle_alphas = []
    
    #ucg_0_angles = [] #commented for opt_qsp_lambda2
    global_phases = []

    for k in range(n): # Multiplexor decomposition
        lambda_diagonals = []
        match k:
            case 0:
                #alpha0 = params_vector[0][0][0] #commented for opt_qsp_lambda2
                #global_phases.append(np.exp(1j*(alpha0))) #commented for opt_qsp_lambda2

                #beta0 = params_vector[0][1][0] #commented for opt_qsp_lambda2
                #print(f"Angle for Rz(beta0) of UCG {k+1}: {beta0}")
                #ucg_0_angles.append(beta0) #commented for opt_qsp_lambda2

                gamma0 = params_vector[0]
                #print(f"Angle for Ry(gamma0) of UCG {k+1}: {gamma0}")
                #ucg_0_angles.append(gamma0) #commented for opt_qsp_lambda2

                #delta0 = params_vector[2][0] #commented for opt_qsp_lambda2
                #print(f"Angle for Rz(delta0) of UCG {k+1}: {delta0}")
                #ucg_0_angles.append(delta0) #commented for opt_qsp_lambda2

                #print("_"*50)
            case _:

                #for i, pv in enumerate(params_vector): #commented for opt_qsp_lambda2

                    diag_rz = np.array([])

                    #match i:
                        #case 0:
                            #for p in range(pow(2, k)):
                                #a_angle = pv[0][pow(2, k) + p - 1]
                                #b_angle = pv[1][pow(2, k) + p - 1]
                                #diag = np.array([np.exp(1j*(a_angle - b_angle)), np.exp(1j*(a_angle + b_angle))])
                                #diag_rz = np.append(diag_rz, diag)
                        #case _:
                    for p in range(pow(2, k)):
                        angle = params_vector[pow(2, k) + p - 1]
                        diag = np.array([np.exp(-1j*(angle)), np.exp(1j*(angle))])
                        diag_rz = np.append(diag_rz, diag)
                    
                    #print(f"Diagonal matrix {i+1}/3 of UCG {k+1}: {diag_rz}")
                    #print(f"Angles for diagonal of matrix {i}/3 of UCG {k+1}: {[np.angle(x) for x in diag_rz]}")

                    lamb = np.array([coeff / diag_rz[0] for coeff in diag_rz])
                    lambda_diagonals.append(lamb)
                    #print(f"Exponential diagonal of Lambda {k+1}.{i+1}: {lamb}")

                    global_phases.append(diag_rz[0])

                    comb = [np.angle(x) for x in lamb]
                    #print(f"Angle combination for UCG {k+1}: {comb}")

                    binary_matrix = generateInnProdMatrix(k+1)
                    alphas = generateAlphaVect(binary_matrix, comb[1:])
                    
                    #print(f"Aphas for Lambda {k+1}.{i+1}: {alphas}")

                    #alpha_check = (pow(2, 1 - (k + 1)) * (2 * binary_matrix - np.ones((pow(2, (k + 1)) - 1,pow(2, (k + 1)) - 1))) @ comb[1:]) # (2^(1-n))*(2*A - J)*thetas

                    #assert np.allclose(alpha_check, alphas, atol = 1e-6), "Alpha values don't match the check: \n" + str(alphas) + " -> " + str(alpha_check)

                    #match i:
                        #case 0:
                            #prefix_alphas.append(alphas)
                        #case 1:
                    middle_alphas.append(alphas)
                        #case 2:
                            #suffix_alphas.append(alphas)

                    #print("-"*50)
                #print("_"*50)

    # PARAMETERS FOR FINAL LAMBDA
    risultati = matrixLambdaFinal(unit_vector)


    # Step 7: Solve the linear system for parallelized phases
    phases_alphas = []

    comb = risultati["New phase angles"]
    #print("\n[Step 7a] Checking the new phase angles that enter the linear system:\n")
    #print(comb)

    binary_matrix = generateInnProdMatrix(n)
    #print("\n[Step 7b] Checking the coefficient matrix:\n")
    #print(binary_matrix)
    phases_alphas = generateAlphaVect(binary_matrix, comb[1:])
    #alpha_check = (pow(2, 1 - n) * (2 * binary_matrix - np.ones((pow(2, n) - 1,pow(2, n) - 1))) @ comb[1:])
    #assert np.allclose(alpha_check, phases_alphas, atol = 1e-6), "Alpha values don't match the check: \n" + str(phases_alphas) + " -> " + str(alpha_check)
    #print("\n[Step 7c] Checking the alphas:\n")
    #print(phases_alphas)
    #print(type(phases_alphas))

    last_global_phase = risultati["Global phase"]
    global_phases.append(last_global_phase)
    last_diagonal = risultati["New complex phases"]
    lambda_diagonals.append(last_diagonal)

    #return prefix_alphas, middle_alphas, suffix_alphas, ucg_0_angles, global_phases, lambda_diagonals #commented for opt_qsp_lambda2
    return middle_alphas, phases_alphas, gamma0, global_phases, lambda_diagonals