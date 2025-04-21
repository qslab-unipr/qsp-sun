from .utils import *
from operator import add
from math import floor, log2
import matplotlib.pyplot as plt
import pennylane as qml
import pennylane.math as qmath

class Stage:
    def __init__(self, qsp_qubits : int, n_lambda : int):
        # QSP CIRCUIT QUIBITS NUMBER
        self._qsp_qubits = qsp_qubits

        # REGISTERS INITIALIZZATION WITH n QUBITS EACH
        self._input_register = n_lambda
        self._copy_register = n_lambda
        self._phase_register = n_lambda

        self._ancillaries = self._copy_register + self._phase_register

        if self._ancillaries != 0 : self._t = int(floor(log2(self._ancillaries / 2)))

    def circuit(self) -> list:
        raise "Abstract method!"
    
class PrefixCopyStage(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int):
        Stage.__init__(self, qsp_qubits, n_lambda)

    def circuit(self) -> list:
        # COPY REGISTER INITIALIZATION ("k" COPIES OF THE FIRST "t" QUBITS)
        k = int(floor(self._ancillaries / (2 * self._t)))
        
        for x in range(self._t):
            control, target = x, self._qsp_qubits + x
            qml.CNOT(wires = [control, target])

        for i in range(k - 1):
            for x in range(self._t):
                control, target = self._qsp_qubits + (self._t * i) + x, self._qsp_qubits + self._t + (self._t * i) + x
                qml.CNOT(wires = [control, target])

class GrayInitialStage(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int, alphas : list):
        Stage.__init__(self, qsp_qubits, n_lambda)
        self._alphas = alphas

    def circuit(self) -> list:
        binary = ['0' * (self._input_register - self._t)] * pow(2, self._t) # Is the same of calculating the first element of Gray code 1 and 2 (00...00)
        basis = [s[::-1] for s in compBasis(self._t)]
        strings = list(map(add, basis, binary)) # Strings in the first column have the last "(n − t)" bits at 0, and strings in each row share the same first "t" bits

        # PHASE REGISTER INITIALIZATION (NP!)
        for i, string in enumerate(strings):
            for j, bit in enumerate(string):
                if bit == '1':
                    control, target = self._qsp_qubits + j, self._qsp_qubits + self._copy_register + i
                    qml.CNOT(wires = [control, target])

        # ALPHA ROTATIONS
        for i, id in enumerate(strings):
            if '1' in id:
                wire = self._qsp_qubits + self._copy_register + i
                qml.PhaseShift(self._alphas[id], wires = [wire])

class SuffixCopyStage(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int):
        Stage.__init__(self, qsp_qubits, n_lambda)

    def circuit(self) -> list:
        # COPY REGISTER RESET (INVERSE OF PREFIX COPY STAGE)
        k = int(floor(self._ancillaries / (2 * self._t)))

        for i in range(k-1)[::-1]:
            for x in range(self._t)[::-1]:
                control, target = self._qsp_qubits + (self._t * i) + x, self._qsp_qubits + self._t + (self._t * i) + x
                qml.CNOT(wires = [control, target])

        for x in range(self._t)[::-1]:
            control, target = x, self._qsp_qubits + x
            qml.CNOT(wires = [control, target])

        # COPY REGISTER INITIALIZATION ("k" COPIES OF THE LAST "n - t" QUBITS)
        k = int(floor(self._ancillaries / (2 * (self._input_register - self._t))))

        for x in range(self._input_register - self._t):
            control, target = self._t + x, self._qsp_qubits + x
            qml.CNOT(wires = [control, target])

        for i in range(k - 1):
            for x in range(self._input_register - self._t):
                control, target = ((self._input_register - self._t) * i) + self._qsp_qubits + x, (self._input_register - self._t) + ((self._input_register - self._t) * i) + self._qsp_qubits + x
                qml.CNOT(wires = [control, target])

class GrayPathStage(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int, alphas : list):
        Stage.__init__(self, qsp_qubits, n_lambda)
        self._alphas = alphas

    def circuit(self) -> list:
        gray_1 = grayCode(1, self._input_register - self._t)
        gray_2 = grayCode(2, self._input_register - self._t)

        basis = [s[::-1] for s in compBasis(self._t)]
        binary = ['0' * (self._input_register - self._t)] * pow(2, self._t)

        prev = list(map(add, basis, binary)) # Vector of Gray Initial Stage binary strings

        for k in range(int(pow(2, self._input_register) / pow(2, self._t) - 1)):

            curr = [] # Current phase string vector

            for i, id in enumerate(basis):
                g1 = str("".join(map(str, gray_1[k + 1])))
                g2 = str("".join(map(str, gray_2[k + 1])))
                curr.append(id + g1 if i < (len(basis) / 2) else id + g2)

            # PHASE REGISTER INITIALIZATION (NP!)
            for i, id in enumerate(basis):
                curr_str = curr[i]
                prev_str = prev[i]
                change_vector = [True if b1 != b2 else False for b1, b2 in zip(prev_str, curr_str)]

                for bit, change in enumerate(change_vector[self._t:]):
                    if change:
                        control, target = self._qsp_qubits + bit, self._qsp_qubits +  self._copy_register + i
                        qml.CNOT(wires = [control, target])

            # ALPHA ROTATIONS
            for i, id in enumerate(basis):
                string = curr[i]
                wire = self._qsp_qubits + self._copy_register + i
                qml.PhaseShift(self._alphas[string], wires = [wire])

            prev = curr # Replace previous phase with current phase
        
class InverseStage(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int):
        Stage.__init__(self, qsp_qubits, n_lambda)

    def circuit(self) -> list:
        # INVERSE GRAY PATH STAGE
        gray_1 = grayCode(1, self._input_register - self._t)
        gray_2 = grayCode(2, self._input_register - self._t)

        basis = [s[::-1] for s in compBasis(self._t)]
        binary = ['0' * (self._input_register - self._t)] * pow(2, self._t)
        strings = [list(map(add, basis, binary))]

        for k in range(int(pow(2, self._input_register) / pow(2, self._t) - 1)):

            s = []

            for i, id in enumerate(basis):
                g1 = str("".join(map(str, gray_1[k + 1])))
                g2 = str("".join(map(str, gray_2[k + 1])))
                s.append(id + g1 if i < (len(basis) / 2) else id + g2)

            strings.append(s)
        
        strings = strings[::-1]

        for k in range(int(pow(2, self._input_register) / pow(2, self._t) - 1))[::-1]:

            curr = strings[k]
            post = strings[k + 1]

            for i in range(len(basis))[::-1]:
                curr_str = curr[i]
                post_str = post[i]
                change_vector = [True if b1 != b2 else False for b1, b2 in zip(post_str, curr_str)]

                for bit, change in enumerate(change_vector[self._t:]):
                    if change:
                        control, target = self._qsp_qubits + bit, self._qsp_qubits + self._copy_register + i
                        qml.CNOT(wires = [control, target])

        # INVERSE SUFFIX COPY STAGE
        k = int(floor(self._ancillaries / (2 * (self._input_register - self._t))))

        for i in range(k - 1):
            for x in range(self._input_register - self._t)[::-1]:
                control, target = ((self._input_register - self._t) * i) + self._qsp_qubits + x, (self._input_register - self._t) + ((self._input_register - self._t) * i) + self._qsp_qubits + x
                qml.CNOT(wires = [control, target])
                
        for x in range(self._input_register - self._t)[::-1]:
            control, target = self._t + x, self._qsp_qubits + x
            qml.CNOT(wires = [control, target])

        k = int(floor(self._ancillaries / (2 * self._t)))
        
        for x in range(self._t):
            control, target = x, self._qsp_qubits + x
            qml.CNOT(wires = [control, target])

        for i in range(k - 1):
            for x in range(self._t):
                control, target = (self._t * i) + self._qsp_qubits + x, self._t + (self._t * i) + self._qsp_qubits + x
                qml.CNOT(wires = [control, target])

        # INVERSE GRAY INITIAL STAGE
        strings = list(map(add, basis, binary))

        for i, string in enumerate(strings[::-1]):
            for j, bit in enumerate(string[::-1]):
                if bit == '1':
                    control, target = self._qsp_qubits + len(string) - 1 - j, self._qsp_qubits + self._copy_register + len(strings) - 1 - i
                    qml.CNOT(wires = [control, target])

        # INVERSE PREFIX COPY STAGE
        for i in range(k-1)[::-1]:
            for x in range(self._t)[::-1]:
                control, target = (self._t * i) + self._qsp_qubits + x, self._t + (self._t * i) + self._qsp_qubits + x
                qml.CNOT(wires = [control, target])

        for x in range(self._t)[::-1]:
            control, target = x, self._qsp_qubits + x
            qml.CNOT(wires = [control, target])

class RotY(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int, angle : float, wire : int):
        Stage.__init__(self, qsp_qubits, n_lambda)
        self._angle = angle
        self._wire = wire

    def circuit(self) -> list:
        qml.RY(2*self._angle, self._wire)

class RotZ(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int, angle : float, wire : int):
        Stage.__init__(self, qsp_qubits, n_lambda)
        self._angle = angle
        self._wire = wire

    def circuit(self) -> list:
        qml.RZ(2*self._angle, self._wire)

class UnitPre(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int,  wire : int):
        Stage.__init__(self, qsp_qubits, n_lambda)
        self._wire = wire
        
    def circuit(self) -> list:
        qml.adjoint(qml.S)(self._wire)
        qml.Hadamard(self._wire)

class UnitPost(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int,  wire : int):
        Stage.__init__(self, qsp_qubits, n_lambda)
        self._wire = wire

    def circuit(self) -> list:
        qml.Hadamard(self._wire)
        qml.S(self._wire)

class ControlGate(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int, angle : float, wire : int, cw):
        Stage.__init__(self, qsp_qubits, 0)
        self._angle = angle
        self._wire = wire
        self._cw = cw

    def circuit(self) -> list:
        U = np.matrix([[np.cos(self._angle),-np.sin(self._angle)],[np.sin(self._angle),np.cos(self._angle)]])
        qml.ControlledQubitUnitary(U, self._wire-1, self._wire, self._cw)

class QuantumCircuit():
    def __init__(self, wires : int):
        self._wires = wires
        self._stages = []

        self._device = qml.device('default.qubit', wires)
    
    def addStage(self, stage):
        self._stages.append(stage)
    
    def circuit(self):
        @qml.qnode(self._device)
        def _circuit():
            for s in self._stages:
                s.circuit()
            return qml.state(), qml.expval(qml.PauliZ(0))
        return _circuit

    def printCircuit(self, mode : str, modulo : bool = False) -> None:
        
        circuit = self.circuit()

        match mode:
            case "console":
                print(qml.draw(qnode = circuit, decimals = 4)())

            case "figure":
                qml.draw_mpl(qnode = circuit, decimals = 4)()
                #plt.savefig(f"results/circuit_{self._wires//3}.png", dpi=300, bbox_inches="tight")
                plt.show()

            case "ket":
                state = self.computeQuantumState(modulo=modulo)
                print(toKet(state, self._wires))

            case "matrix":
                matrix = self.computeMatrix(modulo=modulo)

                with open(f"results/matrix_{self._wires//3}_(modulo={modulo}).txt", "w") as f:
                    for row in matrix:
                        f.write(", ".join(map(str, row)) + "\n")
                #print(matrix) # Console output requires a lot of time for large circuits

            case "state":
                sv = self.computeQuantumState(modulo=modulo)
                print(sv)
    
    def computeMatrix(self, modulo : bool = False) -> np.ndarray:

        matrix = qml.matrix(op = self.circuit(), wire_order = list(range(self._wires)))()

        if modulo == True:
            matrix = [np.sqrt(np.real(x)**2 + np.imag(x)**2) for x in matrix]

        return matrix
    
    def computeQuantumState(self, modulo : bool = False):

        state, _ = self.circuit()()
        
        if modulo == True:
            return np.array([np.sqrt(np.real(x)**2 + np.imag(x)**2) for x in state])
        
        return state
    
    def getCircuitInfos(self):
        specs = qml.specs(self.circuit())()
        return specs
    
    def printMetrics(self, vec1 : np.ndarray, vec2 : np.ndarray):
            print(f"\nFidelity: {qmath.fidelity_statevector(vec1, vec2)}")
            print(f"\nTrace distance: {qmath.trace_distance(qmath.reduce_statevector(vec1, [0]), qmath.reduce_statevector(vec2, [0]))}")            

    def printMatrixComparison(self, mat1 : np.ndarray, mat2 : np.ndarray):
        lambda_diff = [(x - y)**2 for x, y in zip([np.angle(i) for i in mat1], [np.angle(j) for j in mat2])]
        print(f"\nMean squared error: {np.mean(lambda_diff)}")