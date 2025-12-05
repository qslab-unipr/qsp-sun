import numpy as np

def cartesian_to_polar(vector):
    magnitudes = np.abs(vector)
    angles = np.angle(vector)
    return magnitudes, angles

def polar_to_cartesian(magnitudes, angles):
    return magnitudes * np.exp(1j * angles)

def verify_conversion(original, reconstructed):
    if np.allclose(original, reconstructed):
        return "ok"
    else:
        return "error"

def save_to_file(filename, complex_vector, magnitudes, angles, verification_result):
    with open(filename, 'w') as f:
        f.write("Vettore Complesso:\n")
        for item in complex_vector:
            f.write(f"{item}\n")
        
        f.write("\nModuli:\n")
        for item in magnitudes:
            f.write(f"{item}\n")
        
        f.write("\nAngoli di Fase:\n")
        for item in angles:
            f.write(f"{item}\n")
        
        f.write("\nVerifica Conversione: ")
        f.write(verification_result + "\n")

def generate_symbolic_matrix(n):
    """ Matrice simbolica con 'a' e 'b' """
    rows = 2**n
    cols = n
    matrix = np.empty((rows, cols), dtype=object)
    
    for j in range(cols):
        block_size = 2**(n - j - 1)
        repeat_count = 2**j
        symbols = np.tile(np.repeat(["a", "b"], block_size), repeat_count)
        matrix[:, j] = symbols
    
    return matrix

def generate_variable_matrix(n):
    """ Matrice dei coefficienti e array delle incognite """
    rows = 2**n
    cols = 3 * (2**n - 1)
    coefficient_matrix = np.zeros((rows, cols))
    variables = []
    
    for i in range(2**n - 1):
        variables.extend([f"x_{i}", f"y_{i}", f"z_{i}"])
    
    index = 0
    for j in range(n):
        block_size = 2**(n - j - 1)
        repeat_count = 2**j
        symbols = np.tile(np.repeat(["a", "b"], block_size), repeat_count)

        
        for i in range(rows):
            col_idx = 3 * index
            if col_idx + 2 < cols:
                if symbols[i] == "a":
                    coefficient_matrix[i, col_idx] = 1
                    coefficient_matrix[i, col_idx+1] = -1
                    coefficient_matrix[i, col_idx+2] = -1
                else:
                    coefficient_matrix[i, col_idx] = 1
                    coefficient_matrix[i, col_idx+1] = 1
                    coefficient_matrix[i, col_idx+2] = -1
            
            if (i + 1) % (2*block_size) == 0:
                index += 1
    
    return coefficient_matrix, variables

def check_variable_usage(coefficient_matrix, variables):
    """ Verifica che tutte le variabili siano state utilizzate nella matrice dei coefficienti """
    used_variables = set()
    for col in range(coefficient_matrix.shape[1]):
        if np.any(coefficient_matrix[:, col] != 0):
            used_variables.add(variables[col])
    
    unused_variables = set(variables) - used_variables
    return unused_variables

def solve_system(complex_vector, n):
    """ Genera e risolve il sistema lineare con una matrice dei coefficienti e un array random di termini noti (solo fasi) """
    coefficient_matrix, variables = generate_variable_matrix(n)

    magnitudes, angles = cartesian_to_polar(complex_vector)

    try:
        solution = np.linalg.lstsq(coefficient_matrix, angles, rcond=None)[0]
    except np.linalg.LinAlgError:
        solution = "Il sistema non ha soluzioni uniche."
    
    unused_variables = check_variable_usage(coefficient_matrix, variables)
    
    return coefficient_matrix, angles, solution, variables, unused_variables
