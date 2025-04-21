# Sun's Algorithm Implementation

## Description
This repository contains all the files related to the first known implementation of Sun's Algorithm for the problem of quantum state preparation of an n-qubit quantum state, using $2n$ ancillary qubits, which is a particular subcase of the range presented in the first theorem of Sun et al. original paper. All modules have been developed in Python with the support of the PennyLane library, which enables quantum computing simulations. The repository includes various files containing classes and functions essential for the complete implementation:
- **utils.py:** Contains all the functions used to operate with vectors, print results, and compute parameters using traditional computing algorithms.
- **circuit\_classes.py:** Contains all the classes used to generate the desired quantum circuit and to simulate quantum states.
- **lambda\_n.py:** Contains a scalable implementation of the $\Lambda_n$ circuit described in the original paper. It also contains a function to test the implementation by comparing the matrix associated with the circuit, with the theoretical matrix produced using algebric calculations.
- **quantum\_state\_preparation.py:** Contains the initialization of the whole quantum state preparation circuit.
- **phases\_linear\_system.py:** Contains all the functions needed to generate the coefficients related to the phases of the QSP.

## Requirements
The project has been implemented using Python 3.12 (Python 3.10 or later should also work) and the latest version (0.38.0) of the *PennyLane* library. The implementation also makes use of the *math*, *operator*, *numpy*, and *matplotlib* modules, which should already be present in the latest version of Python.

## Functions & Methods
- **qspParameters(unit_vector : np.ndarray, n : int):** Takes the parameter vector correspnding to the desire state and the number of qubits of the input register, and returns a tuple containing all the $\alpha$ angles, a list of global phases and the diagonal of the theoretical matrix associated to each $\Lambda$.
- **qspCircuit(N : int, M : int, alphas : list, ucg_0_angles : list):** Takes the number of input and ancillary qubits and all the angles found following Sun's algorithm, and returns the circuit object of quantum state preparation.
- **lambdaCircuit(N : int, M : int, alphas : list):** Takes the number of input and ancillary qubits for the desire $\Lambda$ and a list of $\alpha$ angles, and returns the circuit object of $\Lambda_n$.
- **lambdaTest(N : int, M : int, circuit : 'QuantumCircuit', lambda_diagonals : list):** Takes the number of input and ancillary qubits for the desire $\Lambda$, the circuit of $\Lambda_n$ and the diagonal of the theoretical matrix associated to $\Lambda_n$, and print the comperison metrics between the theoretical matrix and the actual matrix associated to the circuit.
#### Visualize the QSP circuit
- **QuantumCircuit.printCircuit(mode : str, modulo : bool = False):** Accept different modes to visualize aspects of the circuit:
  - *"console"* : Output the circuit diagram on console
  - *"figure"* : Output an interactive plot of the circuit
  - *"ket"* : Print the output state in ket notation without the global phase
  - *"matrix"* : Save of a file the matrix representation of the circuit
  - *"state"* : Print the output quantum state
- **QuantumCircuit.computeMatrix(modulo : bool = False):** Returns the matrix associated to the circuit
- **QuantumCircuit.computeQuantumState(modulo : bool = False):** Run the simulation and returns the output state without global phase
- **QuantumCircuit.getCircuitInfos():** Returns the informations relative to the circuit such as number of gates, depht, ecc.
- **QuantumCircuit.printMetrics(self, vec1 : np.ndarray, vec2 : np.ndarray):** Calculate fidelity and trace distance
- **QuantumCircuit.printMatrixComparison(self, mat1 : np.ndarray, mat2 : np.ndarray):** Calculate the mean squared error between two matrix
  
## Usage
To test the implementation, it is sufficient to run *"main.py"* and follow the instructions printed on terminal console. It is possible to prepare a random , real or complex, N-qubit quantum state or a specific one by selecting a known quantum state or inserting a $2^n$-element vector of $l_2$-norm = 1.
#### Random state exaple
```bash
Specific quantum state [s] or random vector [r]? > r
Number of qubits? > 4
Dense [d] or sparse [s] random state? > d
Real [r] or complex [c] coefficients? > c
```
#### Specific custom state exaple
```bash
Specific quantum state [s] or random vector [r]? > s
Select a known quantum state or insert a custom one:
1) Bell (n=2)
2) GHZ (n=3)
3) GHZ (n=4)
4) W (n=3)
5) W (n=4)
6) Dicke (n=3)
7) Dicke (n=4)
8) Custom
 > 8
Number of qubits? > 2
Insert real coefficients (floats) separated by spaces > 0.5 0.5 0.5 0.5
```
It is also possible to generate and test a single $\Lambda$ circuit by running *"lambda_test.py"* and inserting the nuber of desired input register qubits; plase note that, because the implementetion of the metrics relies on matrix operations, using a value of n greater than 3 could result in long waiting due to the necessity of large memory and high computational power.

#### Quantum Circuit Implementation
Every implemented circuit makes use of some derived classes from the abstract class *Stage*, to define the position of the gates and the values of some useful parameters, and to instantiate a _QuantumCircuit_ object that is basically a collection of stages and allows displaying the circuit and the quantum state in various ways using the **QuantumCircuit.printCircuit(_params_)** method.
It is also possible to directly obtain the output state of the circuit by calling the **QuantumCircuit.computeQuantumState(_params_)** method or to display the circuit metadata using the **QuantumCircuit.getCircuitInfos()** method.

#### Calculation of Parameters
Multiple functions are used to calculate the parameters needed to make all the necessary phase shifts in the computational basis; this process is done using traditional computation and can therefore be executed separately from the quantum computation part. A user can simply call the **qspParameters(_coeff\_vector_, _num\_of\_qubits_)** function to generate $3*(2^n-1)$ $\alpha$ angles from the desired coefficient vector, these angles are then split and divided between the 3 $\Lambda_n$ components used to implement the $UCG_n$.

#### $\Lambda$ Circuit
With the current implementation, it is possible to build the circuit for $\Lambda_n$, with $n$ as a variable, thanks to the adopted stage structure. 
The actual script works by initializing the variable $n$ and then randomly generating the complex coefficient vector, whose values correspond to the probability amplitudes of the desired final quantum state. Then, the necessary angles are calculated from the initial vector and associated with each computational basis string. Finally, using the *Stage* classes, the quantum circuit is constructed and the final state - in this case, the values of the diagonal of the matrix associated with each multiplexor - can be retrieved.
