# ANNOTAZIONI DEL CONFRONTO CON IL PROGETTO ORIGINALE
# [MODIFICATO] = codice esistente riscritto o con comportamento diverso.
# [AGGIUNTO] = nuovo codice; [RIMOSSO/SOSTITUITO] = descrizione del codice tolto.
# [RIORGANIZZATO] = logica esistente spostata o resa riutilizzabile.
# Le motivazioni descrivono lo scopo tecnico deducibile dai due sorgenti:
# non sono dichiarazioni dell'autore sulle sue intenzioni storiche.
# Le annotazioni sono solo commenti: il codice della versione aggiornata
# rimane invariato. I riferimenti di riga indicano i file ricevuti, prima
# dell'inserimento di questi commenti.
# [NOTA] Le firme con N e M sono mantenute; QuantumCircuit usa ancora N+M fili. Il passaggio
# effettivo a N soli fili avviene nei chiamanti aggiornati, che impostano M=0. I blocchi UnitPre,
# UnitPost e la prima RotY non sono stati riscritti qui.
# The optimized version of the Sun et al.'s QSP algorithm is based on a new decomposition
# which prepares the real part of the desired state separately from the complex one.
# For state with real positive coefficients there is an even more reduced version.
# This module contains the definition of the optimized QSP circuit
# for the general complex case and for its reductions.

from circuit_classes import *

def qspCircuitReduced(N : int, M : int, alphas : list, ucg_0_angle : float) -> 'QuantumCircuit':

    # QUANTUM STATE PREPARATION CIRCUIT FOR REAL CASE
    circuit = QuantumCircuit(N + M)

    for n in range(1, N+1): 

        match n:
            case 1:

                ry = RotY(N, 1, ucg_0_angle, 0)
                circuit.addStage(ry)

            case _:

                basis_vectors = compBasis(n)[1:]
                alpha_vector = {basis_vectors[idx]: coeff for idx, coeff in enumerate(alphas[n-2])}
                    
                # UNITATRY OPERATORS (S+, H)
                pre = UnitPre(n, n, n-1)
                circuit.addStage(pre)

                # LAMBDA WITHOUT ANCILLA
                # [MODIFICATO / SOSTITUITO] Sostituiti PrefixCopyStage, GrayInitialStage,
                # SuffixCopyStage, GrayPathStage e InverseStage, con le rispettive cinque addStage,
                # da un unico LambdaNoAncilla. La sintesi delle parita avviene sui qubit di
                # ingresso, eliminando la necessita dei registri di copia e fase. Qui si realizza il
                # blocco Lambda per la preparazione dei moduli.
                lamb = LambdaNoAncilla(n, alpha_vector)
                circuit.addStage(lamb)

                # UNITATRY OPERATORS (H, S)
                post = UnitPost(n, n, n-1)
                circuit.addStage(post)

    return circuit


def qspCircuitReduced_mod(N : int, M : int, alphas_mag : list, alphas_pha : list, ucg_0_angle : float) -> 'QuantumCircuit':

    # QUANTUM STATE PREPARATION CIRCUIT WITH LAMBDA2 CONSTRUCTORS
    circuit = QuantumCircuit(N + M)

    for n in range(1, N+1): 

        match n:
            case 1:

                ry = RotY(N, 1, ucg_0_angle, 0)
                circuit.addStage(ry)

            case _:

                basis_vectors = compBasis(n)[1:]
                alpha_vector1 = {basis_vectors[idx]: coeff for idx, coeff in enumerate(alphas_mag[n-2])}
                    
                # UNITARY OPERATORS (S+, H)
                pre = UnitPre(n, n, n-1)
                circuit.addStage(pre)

                # LAMBDA WITHOUT ANCILLA
                # [MODIFICATO / SOSTITUITO] Sostituiti PrefixCopyStage, GrayInitialStage,
                # SuffixCopyStage, GrayPathStage e InverseStage, con le rispettive cinque addStage,
                # da un unico LambdaNoAncilla. La sintesi delle parita avviene sui qubit di
                # ingresso, eliminando la necessita dei registri di copia e fase. Qui si realizza il
                # blocco Lambda per la preparazione dei moduli.
                lamb = LambdaNoAncilla(n, alpha_vector1)
                circuit.addStage(lamb)

                # UNITARY OPERATORS (H, S)
                post = UnitPost(n, n, n-1)
                circuit.addStage(post)
    
    basis_vectors = compBasis(N)[1:]
    alpha_vector2 = {basis_vectors[idx]: coeff for idx, coeff in enumerate(alphas_pha)}

    # FINAL LAMBDA FOR COMPLEX PHASES
    # [MODIFICATO / SOSTITUITO] Sostituiti PrefixCopyStage, GrayInitialStage, SuffixCopyStage,
    # GrayPathStage e InverseStage, con le rispettive cinque addStage, da un unico LambdaNoAncilla.
    # La sintesi delle parita avviene sui qubit di ingresso, eliminando la necessita dei registri di
    # copia e fase. Qui si applicano le fasi finali dello stato.
    lamb = LambdaNoAncilla(N, alpha_vector2)
    circuit.addStage(lamb)

    return circuit