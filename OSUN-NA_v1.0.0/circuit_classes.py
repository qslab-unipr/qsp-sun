# ANNOTAZIONI DEL CONFRONTO CON IL PROGETTO ORIGINALE
# [MODIFICATO] = codice esistente riscritto o con comportamento diverso.
# [AGGIUNTO] = nuovo codice; [RIMOSSO/SOSTITUITO] = descrizione del codice tolto.
# [RIORGANIZZATO] = logica esistente spostata o resa riutilizzabile.
# Le motivazioni descrivono lo scopo tecnico deducibile dai due sorgenti:
# non sono dichiarazioni dell'autore sulle sue intenzioni storiche.
# Le annotazioni sono solo commenti: il codice della versione aggiornata
# rimane invariato. I riferimenti di riga indicano i file ricevuti, prima
# dell'inserimento di questi commenti.
# [MODIFICATO: intestazione] Il modulo ora descrive entrambe le costruzioni di Lambda. Le classi
# originali con ancilla restano presenti; i costruttori QSP aggiornati usano LambdaNoAncilla.
# This module contains the classes used to generate the quantum circuit, including the ancillary
# and no-ancilla constructions of the $\Lambda$ operator from arXiv:2108.06150.

from utils_mod import *
from operator import add
# [MODIFICATO: import] Aggiunto ceil per dimensionare i bit di sindrome nella costruzione degli
# insiemi di suffissi.
from math import ceil, floor, log2
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


# Compute the rank of a Boolean matrix over F2.
# [AGGIUNTO] Calcola il rango su F2 con eliminazione tramite XOR. Serve a verificare indipendenza
# binaria, che non va valutata con il rango sui reali.
def _gf2_rank(matrix : np.ndarray) -> int:
    reduced = np.asarray(matrix, dtype=np.uint8).copy() % 2

    if reduced.size == 0:
        return 0

    rows, columns = reduced.shape
    rank = 0

    for column in range(columns):
        pivot = next((row for row in range(rank, rows) if reduced[row, column]), None)

        if pivot is None:
            continue

        if pivot != rank:
            reduced[[rank, pivot]] = reduced[[pivot, rank]]

        for row in range(rows):
            if row != rank and reduced[row, column]:
                reduced[row] ^= reduced[rank]

        rank += 1
        if rank == rows:
            break

    return rank


# Compute the inverse of an invertible Boolean matrix over F2.
# [AGGIUNTO] Inverte una matrice binaria con Gauss-Jordan su [A|I]. Serve a passare tra basi di
# parita e a ripristinare la base iniziale; segnala una matrice singolare.
def _gf2_inverse(matrix : np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.uint8) % 2
    size = matrix.shape[0]
    augmented = np.concatenate((matrix.copy(), np.eye(size, dtype=np.uint8)), axis=1)

    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row, column]), None)

        if pivot is None:
            raise ValueError("The Boolean matrix is not invertible.")

        if pivot != column:
            augmented[[column, pivot]] = augmented[[pivot, column]]

        for row in range(size):
            if row != column and augmented[row, column]:
                augmented[row] ^= augmented[column]

    return augmented[:, size:]


# Check whether a collection of Boolean vectors is linearly independent.
# [AGGIUNTO] Converte gli interi in vettori di bit e confronta rango e numero di vettori. Evita di
# scegliere suffissi linearmente dipendenti.
def _is_independent(vectors : list[int], width : int) -> bool:
    if not vectors:
        return True

    matrix = np.array([
        [(vector >> (width - 1 - bit)) & 1 for bit in range(width)]
        for vector in vectors
    ], dtype=np.uint8)

    return _gf2_rank(matrix) == len(vectors)


# Extend an independent collection to a basis of F2^n.
# [AGGIUNTO] Scarta zero, duplicati e dipendenze, poi completa la base di F2^width. Occorre una base
# completa per ottenere trasformazioni reversibili realizzabili con CNOT.
def _extend_basis(vectors : list[int], width : int) -> tuple[int, ...]:
    basis = []

    for vector in vectors:
        if vector != 0 and vector not in basis and _is_independent(basis + [vector], width):
            basis.append(vector)

    for vector in range(1, pow(2, width)):
        if len(basis) == width:
            break

        if _is_independent(basis + [vector], width):
            basis.append(vector)

    if len(basis) != width:
        raise ValueError("Unable to complete a basis over F2.")

    return tuple(basis)


# Build the linearly independent suffix sets used by the no-ancilla construction.
# [AGGIUNTO] Costruisce famiglie di basi usando sindromi e vicinati binari; verifica che coprano
# tutti i suffissi non nulli. Le basi possono sovrapporsi: LambdaNoAncilla evita di applicare due
# volte le stesse fasi. Il caso width=1 e trattato separatamente.
def _independent_suffix_sets(width : int) -> tuple[tuple[int, ...], ...]:
    if width == 1:
        return ((1,),)

    syndrome_bits = ceil(log2(width + 1))
    all_ones = pow(2, syndrome_bits) - 1
    standard_basis = [1 << (width - 1 - bit) for bit in range(width)]
    l_zero = []
    l_one = []

    for vector in range(pow(2, width)):
        syndrome = 0

        for bit in range(width):
            if (vector >> (width - 1 - bit)) & 1:
                syndrome ^= bit + 1

        if syndrome == 0:
            l_zero.append(vector)
        if syndrome == all_ones:
            l_one.append(vector)

    set_l = list(dict.fromkeys(l_zero + l_one))
    suffix_sets = [tuple(standard_basis)]

    for vector in set_l:
        if vector == 0:
            continue

        neighborhood = [vector ^ basis_vector for basis_vector in standard_basis]
        first_part = []

        for candidate in neighborhood:
            if candidate != 0 and _is_independent(first_part + [candidate], width):
                first_part.append(candidate)

            if len(first_part) == width - 1:
                break

        second_part = [candidate for candidate in neighborhood if candidate not in first_part and candidate != 0]

        if vector not in second_part:
            second_part.append(vector)

        suffix_sets.append(_extend_basis(first_part, width))
        suffix_sets.append(_extend_basis(second_part, width))

    covered = {vector for suffix_set in suffix_sets for vector in suffix_set}
    expected = set(range(1, pow(2, width)))

    if covered != expected:
        raise ValueError("The independent suffix sets do not cover F2^n - {0}.")

    return tuple(suffix_sets)


# Convert an integer into a Boolean row vector.
# [AGGIUNTO] Espande un intero in bit dal piu significativo al meno significativo, mantenendo la
# convenzione usata per stringhe computazionali e matrici di base.
def _binary_vector(value : int, width : int) -> list[int]:
    return [(value >> (width - 1 - bit)) & 1 for bit in range(width)]


# Return the CNOT sequence implementing an invertible linear map over F2.
# [AGGIUNTO] Traduce una mappa lineare invertibile su F2 in CNOT. Ogni XOR di righe corrisponde a
# una CNOT; uno scambio usa tre CNOT. La sequenza di riduzione a identita viene invertita per
# realizzare la mappa richiesta.
def _linear_cnot_sequence(matrix : np.ndarray, wires : list[int]) -> list[tuple[int, int]]:
    reduced = np.asarray(matrix, dtype=np.uint8).copy() % 2
    size = len(wires)
    sequence = []

    for column in range(size):
        pivot = next((row for row in range(column, size) if reduced[row, column]), None)

        if pivot is None:
            raise ValueError("The Boolean matrix is not invertible.")

        if pivot != column:
            reduced[column] ^= reduced[pivot]
            sequence.append((wires[pivot], wires[column]))
            reduced[pivot] ^= reduced[column]
            sequence.append((wires[column], wires[pivot]))
            reduced[column] ^= reduced[pivot]
            sequence.append((wires[pivot], wires[column]))

        for row in range(size):
            if row != column and reduced[row, column]:
                reduced[row] ^= reduced[column]
                sequence.append((wires[column], wires[row]))

    return sequence[::-1]


# Return the zero-based ruler-function value used by the shifted Gray codes.
# [AGGIUNTO] Individua il bit meno significativo a 1. Determina quale controllo cambiare a ogni
# passo dei percorsi Gray traslati, senza costruire la tabella Gray completa.
def _ruler_index(value : int) -> int:
    return (value & -value).bit_length() - 1


# Check whether a rotation angle is numerically zero.
# [AGGIUNTO] Introduce la soglia assoluta 1e-12 per omettere rotazioni numericamente trascurabili. E
# una semplificazione numerica: angoli piccoli ma non nulli vengono approssimati a zero.
def _is_zero_angle(angle : float, tolerance : float = 1e-12) -> bool:
    return abs(float(angle)) <= tolerance


# Normalize a phase angle and remove values equivalent to the identity.
# [AGGIUNTO] Riduce una fase in [-pi, pi) e azzera i residui entro tolleranza. Per PhaseShift,
# angoli che differiscono di 2*pi rappresentano la stessa porta.
def _normalize_phase(angle : float, tolerance : float = 1e-12) -> float:
    normalized = (float(angle) + np.pi) % (2 * np.pi) - np.pi
    return 0.0 if abs(normalized) <= tolerance else normalized


# Merge phase rotations and cancel adjacent identical CNOT gates.
# [AGGIUNTO] Accumula PhaseShift consecutive per filo, riduce gli angoli modulo 2*pi e cancella
# coppie di CNOT identiche adiacenti. Ripete fino a stabilizzazione, poiche una cancellazione puo
# renderne possibili altre. Le fasi vengono scaricate prima di ogni CNOT: non si assume che
# commutino con essa.
def _optimize_lambda_operations(operations : list[tuple], tolerance : float = 1e-12) -> list[tuple]:
    def simplify_once(source : list[tuple]) -> list[tuple]:
        simplified = []
        pending_phases = {}
        pending_order = []

        def flush_phases() -> None:
            for wire in pending_order:
                angle = _normalize_phase(pending_phases[wire], tolerance)
                if angle != 0.0:
                    simplified.append(("PhaseShift", wire, angle))
            pending_phases.clear()
            pending_order.clear()

        for operation in source:
            if operation[0] == "PhaseShift":
                _, wire, angle = operation
                if _is_zero_angle(angle, tolerance):
                    continue
                if wire not in pending_phases:
                    pending_phases[wire] = 0.0
                    pending_order.append(wire)
                pending_phases[wire] += float(angle)
                continue

            flush_phases()
            if simplified and simplified[-1] == operation and operation[0] == "CNOT":
                simplified.pop()
            else:
                simplified.append(operation)

        flush_phases()
        return simplified

    optimized = list(operations)
    while True:
        simplified = simplify_once(optimized)
        if simplified == optimized:
            return simplified
        optimized = simplified


# [AGGIUNTO: intera classe] Nuova sintesi di Lambda mediante sole CNOT e PhaseShift sui qubit di
# ingresso. Sostituisce nei costruttori QSP i cinque stadi con registri ausiliari. La divisione
# controllo/target riguarda lo stesso registro, senza aggiungere qubit. Non eredita Stage: espone
# comunque circuit(), richiesto da QuantumCircuit.
class LambdaNoAncilla:
    def __init__(self, n_lambda : int, alphas : dict[str, float], wires : list[int] | None = None):
        self._n_lambda = n_lambda
        self._alphas = alphas
        # [AGGIUNTO] Consente di scegliere i fili; in assenza di una mappa usa i primi n_lambda
        # qubit.
        self._wires = list(range(n_lambda)) if wires is None else wires

    def circuit(self) -> None:
        # [AGGIUNTO] Emette in PennyLane la lista di operazioni gia sintetizzata e semplificata.
        for operation in self.operations():
            if operation[0] == "CNOT":
                _, control, target = operation
                qml.CNOT(wires=[control, target])
            else:
                _, wire, angle = operation
                qml.PhaseShift(angle, wires=wire)

    def operations(self) -> list[tuple]:
        operations = []
        # [AGGIUNTO] Separa la sintesi classica dalla creazione delle porte per poter ottimizzare l
        # intera sequenza prima di emetterla.
        self._build(self._n_lambda, self._alphas, self._wires, operations)
        return _optimize_lambda_operations(operations)

    def _append_linear_map(self, matrix : np.ndarray, wires : list[int], operations : list[tuple]) -> None:
        for control, target in _linear_cnot_sequence(matrix, wires):
            operations.append(("CNOT", control, target))

    def _build(self, n_lambda : int, alphas : dict[str, float], wires : list[int], operations : list[tuple]) -> None:
        # [AGGIUNTO] Un blocco con tutti gli alpha trascurabili e omesso; per un solo qubit basta la
        # PhaseShift associata a 1.
        if not any(not _is_zero_angle(angle) for angle in alphas.values()):
            return

        if n_lambda == 1:
            angle = alphas.get('1', 0.0)
            if not _is_zero_angle(angle):
                operations.append(("PhaseShift", wires[0], angle))
            return

        # [AGGIUNTO] Divide i qubit in ceil(n/2) controlli e floor(n/2) target. Sui target si
        # rappresentano le parita dei suffissi; sui controlli si percorrono i prefissi.
        target_size = n_lambda // 2
        control_size = n_lambda - target_size
        control_wires = wires[:control_size]
        target_wires = wires[control_size:]
        suffix_sets = _independent_suffix_sets(target_size)
        previous_basis = np.eye(target_size, dtype=np.uint8)
        seen_suffixes = set()

        # [AGGIUNTO] Seleziona i suffissi associati ad almeno una fase significativa per saltare
        # intere basi prive di rotazioni utili.
        nonzero_suffixes = {
            suffix
            for suffix in range(1, pow(2, target_size))
            if any(
                not _is_zero_angle(alphas.get(
                    format(prefix, f'0{control_size}b') + format(suffix, f'0{target_size}b'),
                    0.0
                ))
                for prefix in range(pow(2, control_size))
            )
        }

        for suffix_set in suffix_sets:
            # [AGGIUNTO] seen_suffixes impedisce duplicazioni tra basi sovrapposte; una base senza
            # suffissi nuovi attivi viene saltata.
            active_suffixes = [
                suffix in nonzero_suffixes and suffix not in seen_suffixes
                for suffix in suffix_set
            ]

            if not any(active_suffixes):
                continue

            current_basis = np.array([
                _binary_vector(vector, target_size) for vector in suffix_set
            ], dtype=np.uint8)
            # [AGGIUNTO] Passa direttamente dalla base precedente alla successiva tramite B_nuova *
            # B_precedente^-1, evitando un ripristino intermedio per ogni base.
            transition = (current_basis @ _gf2_inverse(previous_basis)) % 2
            self._append_linear_map(transition, target_wires, operations)

            for index, suffix in enumerate(suffix_set):
                if active_suffixes[index]:
                    # [AGGIUNTO] Applica prima le fasi con prefisso nullo; gli altri prefissi
                    # saranno visitati nel percorso Gray.
                    key = ('0' * control_size) + format(suffix, f'0{target_size}b')
                    angle = alphas.get(key, 0.0)
                    if not _is_zero_angle(angle):
                        operations.append(("PhaseShift", target_wires[index], angle))

            current_prefixes = [0] * target_size

            # [AGGIUNTO] I Gray traslati aggiornano una parita per target con una CNOT per passo e
            # applicano gli alpha corrispondenti, senza copiarli su ancilla.
            for step in range(1, pow(2, control_size)):
                for index in range(target_size):
                    changed_bit = (_ruler_index(step) + index) % control_size
                    operations.append(("CNOT", control_wires[changed_bit], target_wires[index]))
                    current_prefixes[index] ^= 1 << (control_size - 1 - changed_bit)

                for index, suffix in enumerate(suffix_set):
                    if active_suffixes[index]:
                        prefix = format(current_prefixes[index], f'0{control_size}b')
                        key = prefix + format(suffix, f'0{target_size}b')
                        angle = alphas.get(key, 0.0)
                        if not _is_zero_angle(angle):
                            operations.append(("PhaseShift", target_wires[index], angle))

            for index in range(target_size):
                # [AGGIUNTO] Chiude il ciclo Gray per annullare il contributo dei controlli prima
                # del successivo cambio di base.
                changed_bit = (index - 1) % control_size
                operations.append(("CNOT", control_wires[changed_bit], target_wires[index]))

            seen_suffixes.update(
                suffix for index, suffix in enumerate(suffix_set) if active_suffixes[index]
            )
            previous_basis = current_basis

        # [AGGIUNTO] Ripristina la base dei target: Lambda deve lasciare invariati i bit e
        # modificarne solo le fasi.
        self._append_linear_map(_gf2_inverse(previous_basis), target_wires, operations)

        # [AGGIUNTO] Le stringhe con suffisso nullo dipendono solo dai controlli. La ricorsione
        # sintetizza queste fasi sul sottoregistro rimanente.
        recursive_alphas = {
            format(prefix, f'0{control_size}b'): alphas.get(
                format(prefix, f'0{control_size}b') + ('0' * target_size),
                0.0
            )
            for prefix in range(1, pow(2, control_size))
        }
        self._build(control_size, recursive_alphas, control_wires, operations)

class RotY(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int, angle : float, wire : int):
        Stage.__init__(self, qsp_qubits, n_lambda)
        self._angle = angle
        self._wire = wire

    def circuit(self) -> list:
        # [MODIFICATO] Prima la rotazione veniva sempre emessa; ora si omette quando l angolo
        # effettivo 2*angle e entro 1e-12, riducendo le porte quasi identita.
        if not _is_zero_angle(2 * self._angle):
            qml.RY(2 * self._angle, self._wire)

class RotZ(Stage):
    def __init__(self, qsp_qubits : int, n_lambda : int, angle : float, wire : int):
        Stage.__init__(self, qsp_qubits, n_lambda)
        self._angle = angle
        self._wire = wire

    def circuit(self) -> list:
        # [MODIFICATO] Prima la rotazione veniva sempre emessa; ora si omette quando l angolo
        # effettivo 2*angle e entro 1e-12, riducendo le porte quasi identita.
        if not _is_zero_angle(2 * self._angle):
            qml.RZ(2 * self._angle, self._wire)

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
                # [MODIFICATO] Il nome usa il numero reale di fili anziche wires//3: la precedente
                # divisione presupponeva N qubit piu 2N ancilla, ipotesi non piu valida con M=0.
                plt.savefig(f"results/circuit_{self._wires}.png", dpi=300, bbox_inches="tight")
                #plt.show()

            case "ket":
                state = self.computeQuantumState(modulo=modulo)
                print(toKet(state, self._wires))

            case "matrix":
                matrix = self.computeMatrix(modulo=modulo)

                # [MODIFICATO] Anche il nome della matrice usa wires anziche wires//3, coerentemente
                # con il circuito senza ancilla.
                with open(f"results/matrix_{self._wires}_(modulo={modulo}).txt", "w") as f:
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