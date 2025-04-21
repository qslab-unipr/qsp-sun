from .circuit_classes import *

def lambdaCircuit(N : int, M : int, alphas : list) -> 'QuantumCircuit':
    
    # LAMBDA CIRCUIT
    circuit = QuantumCircuit(N + M)

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

    return circuit

def lambdaTest(N : int, M : int, circuit : 'QuantumCircuit', lambda_diagonals : list) -> None:

    lambda_check = np.diag(np.array(lambda_diagonals[1]))

    #input_state_vector = [1] + [0] * (2**N - 1)
    #ancillary_state_vector = [1] + [0] * (2**M - 1)

    #th_state = lambda_check @ input_state_vector
    #th_state = np.kron(th_state, ancillary_state_vector)

    #ac_state = circuit.computeQuantumState(modulo=False)

    #circuit.printMetrics(th_state, ac_state)

    circuit.printMatrixComparison(lambda_check, circuitMatrixReduction(circuit.computeMatrix(modulo=False), N, M)) # Takes very long time for N >= 4