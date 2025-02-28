from math import acos
import numpy as np
from phases_linear_system import *

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

def generateAlphaVect(coeffiecients : np.ndarray, theta_vector : np.ndarray) -> np.ndarray: # Generate the vector of alpha by solving the linear system
    alpha_vector = np.linalg.solve(coeffiecients, theta_vector)
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

def qspParameters(unit_vector : np.ndarray, n : int) -> tuple[list, list, list, list, list, list]:
    
    assert np.isclose(np.linalg.norm(unit_vector), 1, atol = 1e-6), f"Coefficient vector norm is not 1: {np.linalg.norm(unit_vector)} != 1."

    assert (len(unit_vector) == pow(2, n)), f"Wrong number of coefficients 'V' or wrong value of 'n': {len(unit_vector)} != {pow(2, n)}."

    coefficient_matrix, angles, solution, variables, unused_variables = solve_system(unit_vector, n)

    a = solution[0::3]
    b = solution[1::3]
    d = solution[2::3]

    prefix_vector, suffix_vector = (a, b), d

    theta_vector = generateThetaVect(unit_vector, n)  # [(2^n)-1] theta angles

    params_vector = [prefix_vector, theta_vector, suffix_vector]

    prefix_alphas, middle_alphas, suffix_alphas = [], [], [] # Cannot use numpy because the resulting array has different shapes
    ucg_0_angles = []
    global_phases = []

    for k in range(n): # Multiplexor decomposition
        lambda_diagonals = []
        match k:
            case 0:
                alpha0 = params_vector[0][0][0]
                global_phases.append(np.exp(1j*(alpha0)))

                beta0 = params_vector[0][1][0]
                ucg_0_angles.append(beta0)

                gamma0 = params_vector[1][0]
                ucg_0_angles.append(gamma0)

                delta0 = params_vector[2][0]
                ucg_0_angles.append(delta0)

            case _:

                for i, pv in enumerate(params_vector):

                    diag_rz = np.array([])

                    match i:
                        case 0:
                            for p in range(pow(2, k)):
                                a_angle = pv[0][pow(2, k) + p - 1]
                                b_angle = pv[1][pow(2, k) + p - 1]
                                diag = np.array([np.exp(1j*(a_angle - b_angle)), np.exp(1j*(a_angle + b_angle))])
                                diag_rz = np.append(diag_rz, diag)
                        case _:
                            for p in range(pow(2, k)):
                                angle = pv[pow(2, k) + p - 1]
                                diag = np.array([np.exp(-1j*(angle)), np.exp(1j*(angle))])
                                diag_rz = np.append(diag_rz, diag)
                    
                    lamb = np.array([coeff / diag_rz[0] for coeff in diag_rz])
                    lambda_diagonals.append(lamb)

                    global_phases.append(diag_rz[0])

                    comb = [np.angle(x) for x in lamb]

                    binary_matrix = generateInnProdMatrix(k+1)
                    alphas = generateAlphaVect(binary_matrix, comb[1:])
                    
                    alpha_check = (pow(2, 1 - (k + 1)) * (2 * binary_matrix - np.ones((pow(2, (k + 1)) - 1,pow(2, (k + 1)) - 1))) @ comb[1:]) # (2^(1-n))*(2*A - J)*thetas

                    assert np.allclose(alpha_check, alphas, atol = 1e-6), "Alpha values don't match the check: \n" + str(alphas) + " -> " + str(alpha_check)

                    match i:
                        case 0:
                            prefix_alphas.append(alphas)
                        case 1:
                            middle_alphas.append(alphas)
                        case 2:
                            suffix_alphas.append(alphas)

    return prefix_alphas, middle_alphas, suffix_alphas, ucg_0_angles, global_phases, lambda_diagonals