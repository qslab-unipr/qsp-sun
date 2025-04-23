# Implementation of Sun et al.'s Algorithm (v1.0.0)

## Description
This repository documents the first known implementation of the algorithm proposed by Sun et al. (*arXiv 2108.06150*) concerning the problem of quantum state preparation (QSP) of an $n$-qubit quantum state using $m=2n$ ancillary qubits, which is a particular range of the whole parametric domain presented in the first theorem of Sun et al.'s original paper. All modules have been developed in Python with the support of the *PennyLane* library, which enables quantum computing simulations. The repository is organized into various folders and files:
- the folder called *"Circuits"* contains examples of quantum circuits associated with matrices $\Lambda_n$ for low $n$ and the complete QSP circuit for $n=2$.
- the folder called *"qsp_package"* includes all the classes and functions needed for the complete implementation of the QSP algorithm in the chosen range $m=2n$. Conceived as an executive package, it can be used to prepare a quantum state (either random or specific) or to run tests. The Python modules it contains will be described in more detail below.
- the folders called *"qsp_test_results_date"* contains a backup of the test results conducted on a specific date.
- the folders called *"qsp_test_scripts_date"* contains a backup of the test scripts conducted on a specific date.
- the folder called *"single_case_results_date"* contains a backup of the results of individual preparations (particular states and lambda matrix verification) on a specific date.
- **lambda_test.py** is the script to test the implementation of lambda matrices.
- **main.py** is the main script to easily test the execution of the QSP algorithm.

More details about the folder *"qsp_package"*:
- **utils.py** contains all the functions used to operate with vectors, print results, and compute parameters using traditional computing algorithms.
- **circuit\_classes.py** contains all the classes used to generate the desired quantum circuit and to simulate quantum states.
- **lambda\_n.py** contains a scalable implementation of the $\Lambda_n$ circuit described in the original paper. It also contains a function to test the implementation by comparing the matrix associated to the circuit with the theoretical matrix produced using algebric calculations.
- **quantum\_state\_preparation.py** contains the initialization of the whole QSP circuit.
- **phases\_linear\_system.py** contains all the functions needed to generate the classical parameters related to the parallelized phases of the QSP circuit.

## Requirements
The project has been implemented using Python 3.12 (Python 3.10 or later should also work) and the latest version (0.38.0) of the *PennyLane* library. The implementation also makes use of the *math*, *operator*, *numpy*, and *matplotlib* modules, which should already be present in the latest version of Python.

## Functions & Methods
- **qspParameters(unit_vector : np.ndarray, n : int):** takes the parameter vector corresponding to the desire quantum state (normalized) and the number of qubits of the input register, and returns a tuple containing all the $\alpha$ angles (phase shifts parameters), a list of global phases (derived from the construction of the matrices $\Lambda_n$) and the diagonal elements of the matrix associated to each matrix $\Lambda_n$.
- **qspCircuit(N : int, M : int, alphas : list, ucg_0_angles : list):** takes the number of input and ancillary qubits and all the phase shifts found following the stages of Sun et al.'s algorithm, and returns the QSP circuit object.
- **lambdaCircuit(N : int, M : int, alphas : list):** takes the number of input and ancillary qubits for the desire operator $\Lambda_n$ and a list of the corresponding phase shifts, and returns the $\Lambda_n$ circuit object.
- **lambdaTest(N : int, M : int, circuit : 'QuantumCircuit', lambda_diagonals : list):** takes the number of input and ancillary qubits for the desire operator $\Lambda_n$, the circuit associated to $\Lambda_n$ and the diagonal of the theoretical matrix associated to $\Lambda_n$, and print the comparison metrics between the theoretical matrix and the actual matrix associated to the circuit.
#### Visualize the QSP circuit
- **QuantumCircuit.printCircuit(mode : str, modulo : bool = False):** accepts different modes to visualize aspects of the circuit:
  - *"console"* : outputs the circuit diagram on console
  - *"figure"* : outputs an interactive plot of the circuit
  - *"ket"* : prints the output state in ket notation without the global phase
  - *"matrix"* : saves on a file the matrix representation of the circuit
  - *"state"* : prints the output quantum state
- **QuantumCircuit.computeMatrix(modulo : bool = False):** returns the matrix associated to the circuit
- **QuantumCircuit.computeQuantumState(modulo : bool = False):** runs the simulation and returns the output state without global phase
- **QuantumCircuit.getCircuitInfos():** returns the informations relative to the circuit such as the number of gates, depht, etc.
- **QuantumCircuit.printMetrics(self, vec1 : np.ndarray, vec2 : np.ndarray):** calculates fidelity and trace distance
- **QuantumCircuit.printMatrixComparison(self, mat1 : np.ndarray, mat2 : np.ndarray):** calculates the mean squared error between two matrix
  
## Usage
To test the implementation, it is sufficient to run *"main.py"* and follow the instructions printed on terminal console. It is possible to prepare a random , real or complex, N-qubit quantum state or a specific state by selecting a known quantum state or inserting a $2^n$-element vector of $l_2$-norm = 1.
#### Random state example
```bash
Specific quantum state [s] or random vector [r]? > r
Number of qubits? > 4
Dense [d] or sparse [s] random state? > d
Real [r] or complex [c] coefficients? > c
```
#### Specific custom state example
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
It is also possible to generate and test a single $\Lambda_n$ circuit by running *"lambda_test.py"* and inserting the number of desired input register qubits; plase note that, since the metric evaluation is based on matrix operations, using a value of $n$ greater than 3 could result in long waiting due to the necessity of large memory and high computational power.

#### Quantum Circuit Implementation
Every implemented QSP circuit makes use of some derived classes from the abstract class *Stage*, to define the position of the gates and the values of some useful parameters, and to instantiate a _QuantumCircuit_ object that is basically a collection of stages and allows displaying the circuit and the quantum state in various ways using the **QuantumCircuit.printCircuit(_params_)** method.
It is also possible to directly obtain the output state of the circuit by calling the **QuantumCircuit.computeQuantumState(_params_)** method or to display the circuit metadata using the **QuantumCircuit.getCircuitInfos()** method.

#### Calculation of Parameters
Multiple functions are used to calculate the parameters needed to make all the necessary (parallelized) phase shifts in the computational basis; this process is done using traditional computation and can therefore be executed separately from the quantum computation part. A user can simply call the **qspParameters(_coeff\_vector_, _num\_of\_qubits_)** function to generate $3*(2^n-1)$ $\alpha$ angles (phase shifts parameters) from the desired coefficient vector; these latter angles are then split between the 3 matrices $\Lambda_n$ used to implement (according to the decomposition explained in the original paper by Sun et al.) each $UCG_n$ of the whole traditional QSP ladder structure.

#### $\Lambda$ Circuit
With the current implementation, it is possible to build the circuit associated to the operator $\Lambda_n$, with $n$ as a variable, thanks to the adopted stage structure. 
The actual script works by initializing the variable $n$ and then randomly generating the complex coefficient vector, whose values correspond to the probability amplitudes of the desired final quantum state. Then, the necessary angles are calculated from the initial vector and associated with each computational basis string. Finally, using the *Stage* classes, the quantum circuit is constructed and the final state - in this case, the values of the diagonal of the matrix associated with each multiplexor - can be retrieved.
