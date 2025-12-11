# Implementation of the OSUN Algorithm (v1.0.0), the optimized version of SUN

## Description
This repository documents the first optimization of the Sun et al.'s QSP algorithm in the ancillary range $m=2n$, that is, the lower bound of the first parametric range defined by Theorem 1 of [1]. All modules have been developed in Python with the support of the *PennyLane* library, which enables quantum computing simulations. The repository contains the folder named *qsp_package*, which includes all the classes and functions needed for the complete implementation of the OSUN algorithm. Conceived as an executive package, it can be used to prepare a quantum state (either random or specific) or to run tests. The Python modules it contains will be described in more detail below:
- **main\_opt.py** is the main script to easily test the execution of the OSUN algorithm.
- **utils\_mod.py** contains all the functions used to operate with vectors, print results, and compute parameters using traditional computing algorithms.
- **circuit\_classes.py** contains all the classes used to generate the desired quantum circuit and to simulate quantum states.
- **lambda\_final.py** computes and analyzes the parameters for the last $\Lambda$-type constructor, the one dedicated to the preparation of complex phases for each entry of the desired state.
- **quantum\_state\_preparation\_mod.py** contains the initialization of the whole QSP circuit in the reduced and expanded versions.

[1] Sun, Xiaoming, et al. *Asymptotically Optimal Circuit Depth for Quantum State Preparation and General Unitary Synthesis*, arXiv e-prints, 2021, arXiv:2108.06150.

## Requirements
The project has been implemented using Python 3.12 (Python 3.10 or later should also work) and the version 0.38.0 of the *PennyLane* library. The implementation also makes use of the *math*, *operator*, *numpy*, and *matplotlib* modules, which should already be present in the latest version of Python.

## Functions & Methods
- **qspParameters\_mod(unit_vector : np.ndarray, n : int):** takes the parameter vector corresponding to the desire quantum state (normalized) and the number of qubits of the input register, and returns a tuple containing all the $\alpha$ angles (phase shifts parameters), a list of global phases (derived from the construction of the matrices $\Lambda_n$) and the diagonal elements of the matrix associated to each matrix $\Lambda_n$.
- **qspCircuitReduced(N : int, M : int, alphas : list, ucg_0_angles : list):** takes the number of input and ancillary qubits and all the phase shifts found following the OSUN scheme in the real case, and returns the optimized QSP circuit object. 
- **qspCircuitReduced\_mod(N : int, M : int, alphas_mag : list, alphas_pha : list, ucg_0_angle : float):** takes the number of input and ancillary qubits and all the phase shifts found following the OSUN scheme in the complex case, and returns the optimized QSP circuit object.

#### Visualize the QSP circuit
Refer to the repository **SUN_v1.0.0**.
  
## Usage
To test the implementation, it is sufficient to run *"main_opt.py"* and follow the instructions printed on the terminal console. It is possible to prepare a random, real positive, real negative or complex, $n$-qubit quantum state or a specific state by selecting a known quantum state or inserting a $2^n$-element vector of $l_2$-norm = 1.
The script will automatically select the simplest version of the circuit (expanded or reduced) based on the input coefficients vector.

#### Quantum Circuit Implementation
Refer to the repository **SUN_v1.0.0**, keeping in mind that the preparation of the desired state is now split into its real part and its complex part. This allows us to consider only one $\Lambda$-type constructor for each UCG level, adding a highest-order $\Lambda$-operator at the end of the circuit for the complex component only.

#### Calculation of Parameters
The calculation of the parameters is done similarly to the original version **SUN_v1.0.0**.
