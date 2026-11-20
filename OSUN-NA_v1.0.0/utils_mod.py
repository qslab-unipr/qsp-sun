# ANNOTAZIONI DEL CONFRONTO CON IL PROGETTO ORIGINALE
# [MODIFICATO] = codice esistente riscritto o con comportamento diverso.
# [AGGIUNTO] = nuovo codice; [RIMOSSO/SOSTITUITO] = descrizione del codice tolto.
# [RIORGANIZZATO] = logica esistente spostata o resa riutilizzabile.
# Le motivazioni descrivono lo scopo tecnico deducibile dai due sorgenti:
# non sono dichiarazioni dell'autore sulle sue intenzioni storiche.
# Le annotazioni sono solo commenti: il codice della versione aggiornata
# rimane invariato. I riferimenti di riga indicano i file ricevuti, prima
# dell'inserimento di questi commenti.
# [MODIFICATO / RIMOSSO] Riscritto il preprocessing. Eliminato l import math.acos, sostituito da
# np.arccos vettoriale. Rimosse la classe Node e generateBinTree: generateThetaVect lavora
# direttamente sulle probabilita. Tolti anche commenti, stampe disabilitate e resti della vecchia
# decomposizione; questa pulizia da sola non cambia l esecuzione.
# This module contains the vector utilities and the classical preprocessing for QSP.

import numpy as np


# [MODIFICATO] Sostituita la lista temporanea e bin(sum(v))[-1] con somma modulo 2. Le annotazioni
# x,y passano da int a str, coerentemente con l accesso ai singoli bit. Lo scopo matematico resta il
# prodotto interno binario.
def binaryInnProd(x : str, y : str, n : int) -> str:
    return str(sum(int(x[index]) * int(y[index]) for index in range(n)) % 2)


# [MODIFICATO] Lista costruita con comprehension e format invece di append e bin/zfill. Corretto il
# tipo restituito in list[str]; per n positivo conserva l ordine computazionale crescente.
def compBasis(n : int) -> list[str]:
    return [format(value, f"0{n}b") for value in range(pow(2, n))]


# [MODIFICATO] Converte esplicitamente l input in array NumPy e restituisce direttamente modulo e
# argomento; elimina la variabile temporanea angles. La funzione e stata anche spostata nel file.
def cartesian_to_polar(vector):
    vector = np.asarray(vector)
    return np.abs(vector), np.angle(vector)


# [MODIFICATO] Aggiunti include_zeros e tolerance: permette l espansione completa per i report e
# omette normalmente ampiezze entro 1e-12. Prima si escludevano soltanto gli zeri esatti. Usa join,
# quindi rimuove anche gli spazi finali della vecchia rappresentazione.
def toKet(vector : list, n : int, include_zeros : bool = False, tolerance : float = 1e-12) -> str:
    basis = compBasis(n)
    terms = []

    for index, value in enumerate(vector):
        if include_zeros or abs(value) > tolerance:
            terms.append(str(value) + '|' + basis[index] + '>')

    return '  '.join(terms)


# [MODIFICATO] Sostituita la costruzione per concatenazioni con g=i XOR (i>>1) e operazioni
# vettoriali. Aggiunti il caso n<=0 e ValueError per n_gray non valido. ATTENZIONE: cambia anche l
# ordine dei bit; il vecchio tipo 1 costruiva i bit nell ordine opposto, mentre il nuovo usa prima
# il piu significativo. Il tipo 2 ora e una rotazione di colonne. Le classi con ancilla conservate
# possono quindi ricevere sequenze diverse; LambdaNoAncilla usa invece _ruler_index.
def grayCode(n_gray : int, n : int):
    if n <= 0:
        return np.zeros((1, 0), dtype=int)

    values = np.arange(pow(2, n), dtype=np.uint64)
    gray = values ^ (values >> 1)
    bits = ((gray[:, None] >> np.arange(n - 1, -1, -1, dtype=np.uint64)) & 1).astype(int)

    if n_gray == 1:
        return bits
    if n_gray == 2:
        return np.roll(bits, shift=1, axis=1)
    raise ValueError("n_gray must be 1 or 2.")


# [MODIFICATO] Il prodotto di Kronecker parte da np.array([1.0]) anziche dallo scalare 1. Con zero
# ancilla restituisce quindi un vettore di lunghezza 1, mantenendo il fattore neutro del prodotto.
def stateToVector(state : np.ndarray):
    state_vector = np.array([1.0])
    for value in state:
        state_vector = np.kron(state_vector, np.array([1 - value, value]))
    return state_vector


# [MODIFICATO: formato di ritorno] Accetta una matrice oppure la sua diagonale e seleziona un
# elemento ogni 2**m. Prima costruiva una matrice diagonale ridotta 2D; ora restituisce un vettore
# 1D copiato. Evita la matrice densa; il parametro n rimane in firma ma non viene usato.
def circuitMatrixReduction(diagonal, n : int, m : int):
    diagonal = np.asarray(diagonal)
    if diagonal.ndim == 2:
        diagonal = np.diag(diagonal)
    return diagonal[::pow(2, m)].copy()


# MULTIPLEXOR ANGLES
# [MODIFICATO / SOSTITUITO] Sostituisce l albero Node/generateBinTree e le visite ricorsive con
# somme delle probabilita |v|^2 per blocchi. Mantiene la raccolta per livelli, controlla la
# dimensione e applica arccos(sqrt(p_sinistra/p_padre)). Lo scopo e ridurre oggetti Python e accessi
# ricorsivi.
def generateThetaVect(unit_vector : np.ndarray, n : int, tolerance : float = 1e-12) -> np.ndarray:
    vector = np.asarray(unit_vector, dtype=complex)
    expected_size = pow(2, n)

    if vector.size != expected_size:
        raise ValueError(f"Wrong vector size: {vector.size} != {expected_size}.")

    probabilities = np.abs(vector) ** 2
    theta_levels = []
    zero_threshold = tolerance ** 2

    for depth in range(1, n + 1):
        block_size = pow(2, n - depth + 1)
        half_size = block_size // 2
        blocks = probabilities.reshape(-1, block_size)
        parent_norms = blocks.sum(axis=1)
        left_norms = blocks[:, :half_size].sum(axis=1)

        # [MODIFICATO: rami a probabilita nulla] Divide solo se il peso del padre supera
        # tolerance**2. Altrimenti assume rapporto 1 e quindi theta=0; prima, per padre esattamente
        # zero, usava acos(0)=pi/2. Sono rami senza ampiezza, sui quali si evitano rotazioni
        # inutili. Pesi molto piccoli sono ora trattati numericamente come nulli.
        ratios = np.divide(
            left_norms,
            parent_norms,
            out=np.ones_like(left_norms),
            where=parent_norms > zero_threshold
        )
        # [AGGIUNTO] Limita i rapporti a [0,1] per evitare errori di dominio; porta vicino a 0/1 i
        # valori entro tolleranza, stabilizzando gli angoli limite.
        ratios = np.clip(ratios, 0.0, 1.0)
        ratios[ratios <= tolerance] = 0.0
        ratios[1.0 - ratios <= tolerance] = 1.0
        theta_levels.append(np.arccos(np.sqrt(ratios)))

    return np.concatenate(theta_levels) if theta_levels else np.array([], dtype=float)


# [MODIFICATO] Costruisce A usando AND bit a bit e parita di bit_count, evitando stringhe e ripetuti
# np.append. La matrice densa resta disponibile per diagnostica/validazione; non viene piu usata dal
# calcolo veloce degli alpha.
def generateInnProdMatrix(n : int) -> np.ndarray:
    basis = np.arange(1, pow(2, n), dtype=np.uint64)
    return np.array([
        [int((int(left) & int(right)).bit_count() % 2) for right in basis]
        for left in basis
    ], dtype=float)


# ALPHA COEFFICIENTS
# [MODIFICATO: algoritmo e firma] Prima: generateAlphaVect(A, theta), con np.linalg.solve. Ora:
# generateAlphaVect(theta, n), con controllo di 2**n-1 elementi e trasformata Walsh-Hadamard a
# farfalla sul vettore [0, theta]. Il fattore -2**(1-n) restituisce gli alpha delle parita non
# nulle. Con d=2**n, il calcolo richiede O(d log d) operazioni e O(d) memoria invece di costruire e
# risolvere un sistema denso. Questa e una proprieta del calcolo degli alpha, non una misura delle
# prestazioni dell intero QSP.
def generateAlphaVect(theta_vector : np.ndarray, n : int) -> np.ndarray:
    theta_vector = np.asarray(theta_vector, dtype=float)
    expected_size = pow(2, n) - 1

    if theta_vector.size != expected_size:
        raise ValueError(f"Wrong theta vector size: {theta_vector.size} != {expected_size}.")

    transformed = np.zeros(pow(2, n), dtype=float)
    transformed[1:] = theta_vector
    block_size = 1

    while block_size < transformed.size:
        blocks = transformed.reshape(-1, 2 * block_size)
        left = blocks[:, :block_size].copy()
        right = blocks[:, block_size:].copy()
        blocks[:, :block_size] = left + right
        blocks[:, block_size:] = left - right
        block_size *= 2

    return -pow(2, 1 - n) * transformed[1:]


# [MODIFICATO] Sostituisce la ricorsione e np.union1d con generazione diretta, livello per livello,
# delle stringhe che terminano in zero. binary_vector resta nella firma ma non viene usato: input
# iniziali diversi da quello standard non mantengono necessariamente la vecchia semantica.
# generateThetaVect non chiama piu questa funzione.
def encodeMultiplexor(binary_vector : np.ndarray, n : int) -> np.ndarray:
    strings = []
    for depth in range(1, n + 1):
        strings.extend(format(parent << 1, f"0{depth}b") for parent in range(pow(2, depth - 1)))
    return np.asarray(strings)


# [AGGIUNTO] Porta gli angoli in [-pi, pi) e azzera i residui numerici. Serve a trattare
# uniformemente fasi periodiche e identita approssimate.
def _wrap_angles(angles : np.ndarray, tolerance : float = 1e-12) -> np.ndarray:
    wrapped = (np.asarray(angles, dtype=float) + np.pi) % (2 * np.pi) - np.pi
    wrapped[np.abs(wrapped) <= tolerance] = 0.0
    return wrapped


# [AGGIUNTO] Decide se includere Lambda finale ignorando ampiezze sotto soglia e angoli quasi nulli.
# Riconosce anche coefficienti reali negativi (fase pi). Non elimina in generale una fase globale
# comune: puo selezionare Lambda anche se tutte le ampiezze non nulle hanno la stessa fase.
def has_nontrivial_phases(state_vector, tolerance : float = 1e-12) -> bool:
    vector = np.asarray(state_vector, dtype=complex)
    mask = np.abs(vector) > tolerance
    if not np.any(mask):
        return False
    return bool(np.any(np.abs(_wrap_angles(np.angle(vector[mask]), tolerance)) > tolerance))


# FINAL PHASES
# [MODIFICATO] Vettorializza l estrazione delle fasi con una maschera e rende configurabile la
# soglia (prima 1e-10, ora 1e-12 di default). Conserva la divisione per phase0 e usa _wrap_angles.
# Rimuove la chiave "Diagonal matrix": il risultato contiene il vettore delle fasi, evitando np.diag
# e la sua memoria quadratica. Eventuali chiamanti di quella chiave devono essere adattati.
def matrixLambdaFinal(state_vector, tolerance : float = 1e-12):
    vector = np.asarray(state_vector, dtype=complex)
    magnitudes = np.abs(vector)
    phases = np.ones(vector.size, dtype=complex)
    mask = magnitudes > tolerance
    phases[mask] = vector[mask] / magnitudes[mask]

    phase0 = phases[0]
    new_phases = phases / phase0

    return {
        "Magnitudes": magnitudes,
        "New complex phases": new_phases,
        "New phase angles": _wrap_angles(np.angle(new_phases), tolerance),
        "Global phase": phase0
    }


# QSP PARAMETERS
# [MODIFICATO] Aggiunti include_phases e tolerance, validazioni attive di dimensione/norma e tipi
# restituiti piu precisi. Riorganizza i livelli senza i vecchi rami commentati; usa gli alpha veloci
# e conserva le diagonali ideali per la validazione. Restituisce sempre cinque valori.
def qspParameters_mod(
    unit_vector : np.ndarray,
    n : int,
    include_phases : bool = True,
    tolerance : float = 1e-12
) -> tuple[list, np.ndarray, float, list, list]:
    vector = np.asarray(unit_vector, dtype=complex)
    expected_size = pow(2, n)

    if vector.size != expected_size:
        raise ValueError(f"Wrong vector size: {vector.size} != {expected_size}.")
    if not np.isclose(np.linalg.norm(vector), 1.0, atol=1e-8):
        raise ValueError("Coefficient vector norm is not 1.")

    theta_vector = generateThetaVect(vector, n, tolerance)
    gamma0 = float(theta_vector[0]) if theta_vector.size else 0.0
    middle_alphas = []
    global_phases = []
    # [MODIFICATO] Inizializza una sola volta l elenco. Nell originale era azzerato a ogni
    # iterazione del ciclo k: si perdevano le diagonali dei livelli precedenti. Qui vengono
    # conservate tutte per i controlli dei singoli Lambda.
    lambda_diagonals = []

    for level in range(1, n):
        start = pow(2, level) - 1
        stop = pow(2, level + 1) - 1
        angles = theta_vector[start:stop]
        reference = float(angles[0])

        # [MODIFICATO] Calcola direttamente reference-angle e reference+angle, poi applica il
        # wrapping. Sostituisce la costruzione di esponenziali, la divisione per il primo e la
        # successiva estrazione degli argomenti; conserva la convenzione di fase globale
        # exp(-i*reference).
        combinations = np.empty(pow(2, level + 1), dtype=float)
        combinations[0::2] = reference - angles
        combinations[1::2] = reference + angles
        combinations = _wrap_angles(combinations, tolerance)

        middle_alphas.append(generateAlphaVect(combinations[1:], level + 1))
        global_phases.append(np.exp(-1j * reference))
        lambda_diagonals.append(np.exp(1j * combinations))

    # [MODIFICATO] La Lambda delle fasi era sempre precomputata; ora il chiamante puo saltarla. Nel
    # ramo False phases_alphas e vuoto e non si aggiungono ne la relativa fase globale ne la
    # diagonale finale.
    if include_phases:
        final_parameters = matrixLambdaFinal(vector, tolerance)
        phases_alphas = generateAlphaVect(final_parameters["New phase angles"][1:], n)
        global_phases.append(final_parameters["Global phase"])
        lambda_diagonals.append(
            np.asarray(final_parameters["New complex phases"], dtype=complex)
        )
    else:
        phases_alphas = np.array([], dtype=float)

    return middle_alphas, phases_alphas, gamma0, global_phases, lambda_diagonals
