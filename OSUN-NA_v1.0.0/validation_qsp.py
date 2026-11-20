# ANNOTAZIONI DEL CONFRONTO CON IL PROGETTO ORIGINALE
# [MODIFICATO] = codice esistente riscritto o con comportamento diverso.
# [AGGIUNTO] = nuovo codice; [RIMOSSO/SOSTITUITO] = descrizione del codice tolto.
# [RIORGANIZZATO] = logica esistente spostata o resa riutilizzabile.
# Le motivazioni descrivono lo scopo tecnico deducibile dai due sorgenti:
# non sono dichiarazioni dell'autore sulle sue intenzioni storiche.
# Le annotazioni sono solo commenti: il codice della versione aggiornata
# rimane invariato. I riferimenti di riga indicano i file ricevuti, prima
# dell'inserimento di questi commenti.
# [FILE INTERAMENTE AGGIUNTO] Introduce controlli sui singoli Lambda: confronto con la diagonale
# ideale e prova di evoluzione su uno stato casuale. I test documentano il caso selezionato, non una
# dimostrazione generale per ogni N; le diagonali ideali provengono dal preprocessing corrente.
"""Validation of the Lambda blocks used by the selected QSP circuit."""

# [AGGIUNTO: import] Supporta annotazioni di tipo posticipate; re/Path/shutil gestiscono output,
# NumPy le matrici, PennyLane i circuiti e Matplotlib le immagini.
from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pennylane as qml

from circuit_classes import LambdaNoAncilla
from utils_mod import compBasis, generateInnProdMatrix


# [AGGIUNTO] Fissa atol=1e-9 e limita le stampe estese a N<=4. Il limite di stampa non evita la
# costruzione delle matrici dense nei test.
DEFAULT_ATOL = 1e-9
PRINT_FULL_MATRICES_MAX_QUBITS = 4


# [AGGIUNTO] Produce un nome di cartella per il percorso predefinito del validatore; main_opt passa
# invece una directory esplicita.
def _slugify(value: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9_-]+", "_", value.strip()).strip("_").lower()
    return slug or "qsp_state"


# [AGGIUNTO] Verifica i 2**N-1 alpha e li associa alle stringhe non nulle nell ordine
# computazionale, adattando il vettore alla firma di LambdaNoAncilla.
def alpha_vector_to_dict(
    alpha_vector: Sequence[float],
    n_qubits: int,
) -> dict[str, float]:
    values = np.asarray(alpha_vector, dtype=float).reshape(-1)
    expected = 2**n_qubits - 1

    if values.size != expected:
        raise ValueError(
            f"Wrong alpha count for Lambda n={n_qubits}: "
            f"{values.size} != {expected}."
        )

    basis = compBasis(n_qubits)[1:]
    return {basis_state: float(alpha) for basis_state, alpha in zip(basis, values)}


# [AGGIUNTO] Isola un singolo Lambda su N fili per estrarne matrice, disegno e risorse.
def _make_lambda_qnode(n_qubits: int, alphas: Mapping[str, float]):
    device = qml.device("default.qubit", wires=n_qubits)
    stage = LambdaNoAncilla(n_qubits, dict(alphas))

    @qml.qnode(device)
    def lambda_qnode():
        stage.circuit()
        return qml.state()

    return lambda_qnode


# [AGGIUNTO] Prepara uno stato arbitrario con StatePrep e applica Lambda. Serve a controllarne l
# azione oltre alla sola estrazione della matrice; StatePrep e parte del test, non della sintesi
# Lambda misurata separatamente.
def _make_lambda_evolution_qnode(n_qubits: int, alphas: Mapping[str, float]):
    device = qml.device("default.qubit", wires=n_qubits)
    stage = LambdaNoAncilla(n_qubits, dict(alphas))
    wires = list(range(n_qubits))

    @qml.qnode(device)
    def evolve(initial_state):
        qml.StatePrep(initial_state, wires=wires)
        stage.circuit()
        return qml.state()

    return evolve


# [AGGIUNTO] Genera uno stato complesso normalizzato con seed esplicito per riprodurre il controllo
# di evoluzione.
def _random_state(n_qubits: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    state = rng.normal(size=2**n_qubits) + 1j * rng.normal(size=2**n_qubits)
    state /= np.linalg.norm(state)
    return state.astype(complex)


# [AGGIUNTO] Stima la fase relativa con l overlap dei vettori appiattiti, anche per matrici, per
# distinguere uguaglianza esatta ed equivalenza a fase globale.
def _global_phase(reference: np.ndarray, candidate: np.ndarray) -> complex:
    overlap = np.vdot(reference.reshape(-1), candidate.reshape(-1))
    if abs(overlap) <= 1e-15:
        return 1.0 + 0.0j
    return overlap / abs(overlap)


# [AGGIUNTO] Usa l errore massimo elemento per elemento, evidenziando anche una discrepanza
# localizzata.
def _max_abs_error(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(first) - np.asarray(second))))


# [AGGIUNTO] Estrae le risorse in un dizionario e gestisce campi non disponibili, per includerle nei
# risultati della validazione.
def _resource_dict(qnode) -> dict:
    specs = qml.specs(qnode)()
    resources = specs.get("resources")

    if resources is None:
        return {}

    result = {
        "num_wires": getattr(resources, "num_wires", None),
        "num_gates": getattr(resources, "num_gates", None),
        "depth": getattr(resources, "depth", None),
        "gate_types": dict(getattr(resources, "gate_types", {}) or {}),
        "gate_sizes": dict(getattr(resources, "gate_sizes", {}) or {}),
    }

    return {key: value for key, value in result.items() if value is not None}


# [AGGIUNTO] Salva un immagine per blocco e chiude la figura per non accumulare finestre/oggetti
# grafici.
def _draw_lambda_png(qnode, image_path: Path) -> None:
    fig, _ = qml.draw_mpl(qnode, decimals=4, show_all_wires=True)()
    fig.savefig(image_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# [AGGIUNTO] Cancella ricorsivamente la directory di validazione scelta e la ricrea. Scopo: non
# lasciare immagini di blocchi assenti nella nuova esecuzione; i precedenti contenuti di quella
# directory vengono rimossi.
def _prepare_output_directory(output_dir: Path) -> None:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


# [AGGIUNTO] Stampa forma e, per N piccolo, tutti gli elementi: limita la verbosita del report.
def _print_matrix(name: str, matrix: np.ndarray, n_qubits: int) -> None:
    print(f"{name}: shape={matrix.shape}")
    if n_qubits <= PRINT_FULL_MATRICES_MAX_QUBITS:
        print(matrix)
    else:
        print(f"  Matrix omitted for n>{PRINT_FULL_MATRICES_MAX_QUBITS}.")


# [AGGIUNTO] Mostra il confronto delle componenti con etichette computazionali per N piccolo e lo
# omette per dimensioni maggiori.
def _print_component_comparison(reference: np.ndarray, candidate: np.ndarray, n_qubits: int) -> None:
    print("Component comparison:")
    if n_qubits <= PRINT_FULL_MATRICES_MAX_QUBITS:
        for index, (expected, actual) in enumerate(zip(reference, candidate)):
            bit_string = format(index, f"0{n_qubits}b")
            error = abs(expected - actual)
            print(
                f"  |{bit_string}>  expected={expected!s:>24}  "
                f"circuit={actual!s:>24}  |diff|={error:.3e}"
            )
    else:
        print(f"  Component table omitted for n>{PRINT_FULL_MATRICES_MAX_QUBITS}.")


# [AGGIUNTO] Verifica dimensioni, costruisce matrice ideale e matrice del circuito, confronta
# entrambe e applica i due operatori allo stesso stato casuale. Restituisce errori, esiti esatti e a
# fase globale, fidelita, risorse e percorso immagine.
def validate_lambda(
    *,
    n_qubits: int,
    alphas: Sequence[float] | Mapping[str, float],
    ideal_diagonal: Sequence[complex],
    label: str,
    seed: int,
    atol: float,
    output_dir: Path,
) -> dict:
    if isinstance(alphas, Mapping):
        alpha_dict = {key: float(value) for key, value in alphas.items()}
    else:
        alpha_dict = alpha_vector_to_dict(alphas, n_qubits)

    diagonal = np.asarray(ideal_diagonal, dtype=complex).reshape(-1)
    expected_dimension = 2**n_qubits

    if diagonal.size != expected_dimension:
        raise ValueError(
            f"Wrong ideal diagonal size for Lambda n={n_qubits}: "
            f"{diagonal.size} != {expected_dimension}."
        )

    print("\n" + "-" * 78)
    print(f"LAMBDA CHECK: {label}")
    print("-" * 78)

    target_size = n_qubits // 2
    control_size = n_qubits - target_size
    control_wires = list(range(control_size))
    target_wires = list(range(control_size, n_qubits))
    print(f"Lambda size: n={n_qubits}")
    print(f"Register split: r_c={control_size}, r_t={target_size}")
    print(f"Control wires: {control_wires}")
    print(f"Target wires:  {target_wires}")

    print("\nPhase rotations associated with computational strings:")
    for bit_string in compBasis(n_qubits)[1:]:
        print(f"  {bit_string} -> alpha = {alpha_dict.get(bit_string, 0.0): .12f} rad")

    # [AGGIUNTO] Ricostruisce A solo per diagnostica: il preprocessing veloce non la usa. Questo
    # test conserva un costo di memoria quadratico in 2**N.
    inner_product_matrix = generateInnProdMatrix(n_qubits)
    print("\nBinary inner-product matrix A:")
    if n_qubits <= PRINT_FULL_MATRICES_MAX_QUBITS:
        print(inner_product_matrix)
    else:
        print(f"  A shape={inner_product_matrix.shape}")

    # [AGGIUNTO] Materializza la matrice ideale diagonale per confrontarla con qml.matrix. Anche
    # quando le stampe sono omesse, entrambe le matrici dense vengono calcolate.
    ideal_matrix = np.diag(diagonal)
    print("\nIdeal Lambda diagonal:")
    print(diagonal)
    _print_matrix("Ideal Lambda matrix", ideal_matrix, n_qubits)

    lambda_qnode = _make_lambda_qnode(n_qubits, alpha_dict)
    circuit_matrix = np.asarray(
        qml.matrix(lambda_qnode, wire_order=list(range(n_qubits)))(),
        dtype=complex,
    )
    _print_matrix("Lambda matrix extracted from circuit", circuit_matrix, n_qubits)

    # [AGGIUNTO] Registra sia il confronto diretto sia quello dopo rimozione della fase globale; gli
    # allclose usano rtol=0 per imporre la sola tolleranza assoluta.
    matrix_phase = _global_phase(ideal_matrix, circuit_matrix)
    aligned_circuit_matrix = circuit_matrix / matrix_phase
    matrix_error = _max_abs_error(ideal_matrix, circuit_matrix)
    matrix_phase_error = _max_abs_error(ideal_matrix, aligned_circuit_matrix)
    matrix_exact_pass = bool(np.allclose(ideal_matrix, circuit_matrix, atol=atol, rtol=0.0))
    matrix_physical_pass = bool(
        np.allclose(ideal_matrix, aligned_circuit_matrix, atol=atol, rtol=0.0)
    )

    print("\nMatrix implementation test:")
    print(f"  Result (atol={atol:g}): {'PASS' if matrix_exact_pass else 'FAIL'}")
    print(f"  Maximum matrix error: {matrix_error:.3e}")

    # [AGGIUNTO] Controlla l azione sull intera sovrapposizione casuale: Lambda su |0...0> da solo
    # non verificherebbe tutte le fasi.
    initial_state = _random_state(n_qubits, seed)
    print("\nRandom-state evolution test:")
    print(f"  Seed: {seed}")
    print(f"  Input-state norm: {np.linalg.norm(initial_state):.16f}")

    matrix_output = ideal_matrix @ initial_state
    evolution_qnode = _make_lambda_evolution_qnode(n_qubits, alpha_dict)
    circuit_output = np.asarray(evolution_qnode(initial_state), dtype=complex)

    state_phase = _global_phase(matrix_output, circuit_output)
    aligned_circuit_output = circuit_output / state_phase
    state_error = _max_abs_error(matrix_output, circuit_output)
    state_phase_error = _max_abs_error(matrix_output, aligned_circuit_output)
    state_exact_pass = bool(np.allclose(matrix_output, circuit_output, atol=atol, rtol=0.0))
    state_physical_pass = bool(
        np.allclose(matrix_output, aligned_circuit_output, atol=atol, rtol=0.0)
    )
    fidelity = float(abs(np.vdot(matrix_output, circuit_output)) ** 2)

    _print_component_comparison(matrix_output, circuit_output, n_qubits)
    print(f"  State evolution test: {'PASS' if state_exact_pass else 'FAIL'}")
    print(f"  Fidelity: {fidelity:.16f}")
    print(f"  Maximum state error: {state_error:.3e}")

    image_path = output_dir / f"{label}.png"
    _draw_lambda_png(lambda_qnode, image_path)
    resources = _resource_dict(lambda_qnode)
    print(f"  Lambda circuit image: {image_path}")
    if resources:
        print(
            "  Lambda resources: "
            f"wires={resources.get('num_wires', '?')}, "
            f"gates={resources.get('num_gates', '?')}, "
            f"depth={resources.get('depth', '?')}, "
            f"gate_types={resources.get('gate_types', {})}"
        )

    return {
        "label": label,
        "n_qubits": n_qubits,
        "image": str(image_path),
        "register_split": {
            "r_c": control_size,
            "r_t": target_size,
            "control_wires": control_wires,
            "target_wires": target_wires,
        },
        "matrix_test": {
            "exact_pass": matrix_exact_pass,
            "pass_up_to_global_phase": matrix_physical_pass,
            "max_abs_error": matrix_error,
            "max_abs_error_up_to_global_phase": matrix_phase_error,
        },
        "random_state_test": {
            "seed": seed,
            "exact_pass": state_exact_pass,
            "pass_up_to_global_phase": state_physical_pass,
            "fidelity": fidelity,
            "max_abs_error": state_error,
            "max_abs_error_up_to_global_phase": state_phase_error,
        },
        "resources": resources,
    }


# [AGGIUNTO] Valida numero/dimensione dei Lambda, prepara la cartella ed esegue i controlli di ogni
# livello dei moduli piu l eventuale Lambda finale. Restituisce anche le risorse del QSP completo;
# la verifica finale dello stato QSP viene fatta in main_opt.
def validate_current_qsp_run(
    *,
    n_qubits: int,
    state_label: str,
    qsp_circuit,
    middle_alphas: Sequence[Sequence[float]],
    phases_alphas: Sequence[float] | np.ndarray | None,
    lambda_diagonals: Sequence[Sequence[complex]],
    has_final_phase_lambda: bool,
    output_root: str | Path = "results/validation_current",
    output_dir: str | Path | None = None,
    seed: int = 2026,
    atol: float = DEFAULT_ATOL,
) -> dict:
    if n_qubits < 1:
        raise ValueError("n_qubits must be positive.")

    expected_middle_count = max(0, n_qubits - 1)
    if len(middle_alphas) != expected_middle_count:
        raise ValueError(
            "Unexpected number of magnitude Lambdas: "
            f"{len(middle_alphas)} != {expected_middle_count}."
        )

    expected_diagonal_count = expected_middle_count + int(has_final_phase_lambda)
    if len(lambda_diagonals) != expected_diagonal_count:
        raise ValueError(
            "The ideal Lambda diagonals do not match the current circuit: "
            f"{len(lambda_diagonals)} != {expected_diagonal_count}."
        )

    if output_dir is None:
        slug = _slugify(state_label)
        output_dir = Path(output_root) / slug
    else:
        output_dir = Path(output_dir)
    _prepare_output_directory(output_dir)

    print("\n")
    print("=" * 78)
    print("LAMBDA VALIDATION")
    print("=" * 78)
    print(f"Selected state: {state_label}")
    print(f"Current QSP size: N={n_qubits}")

    lambda_reports = []

    for index, alpha_vector in enumerate(middle_alphas):
        lambda_n = index + 2
        lambda_reports.append(
            validate_lambda(
                n_qubits=lambda_n,
                alphas=alpha_vector,
                ideal_diagonal=lambda_diagonals[index],
                label=f"lambda_magnitude_n{lambda_n}",
                seed=seed + index,
                atol=atol,
                output_dir=output_dir,
            )
        )

    if has_final_phase_lambda:
        phase_vector = np.asarray(phases_alphas, dtype=float).reshape(-1)
        if phase_vector.size != 2**n_qubits - 1:
            raise ValueError(
                "Wrong alpha count for the final phase Lambda: "
                f"{phase_vector.size} != {2**n_qubits - 1}."
            )

        lambda_reports.append(
            validate_lambda(
                n_qubits=n_qubits,
                alphas=phase_vector,
                ideal_diagonal=lambda_diagonals[-1],
                label=f"lambda_final_phase_n{n_qubits}",
                seed=seed + len(lambda_reports),
                atol=atol,
                output_dir=output_dir,
            )
        )

    # [AGGIUNTO] Il risultato complessivo richiede i test ESATTI di matrice ed evoluzione. Gli esiti
    # a fase globale sono diagnostici e non sostituiscono exact_pass. Per N=1 senza Lambda
    # intermedi, la lista vuota viene accettata: la RotY finale e verificata dal controllo QSP in
    # main_opt.
    all_pass = all(
        report["matrix_test"]["exact_pass"]
        and report["random_state_test"]["exact_pass"]
        for report in lambda_reports
    ) if lambda_reports else n_qubits == 1

    qsp_qnode = qsp_circuit.circuit()
    full_resources = _resource_dict(qsp_qnode)

    print("\n" + "=" * 78)
    print(f"LAMBDA VALIDATION RESULT: {'PASS' if all_pass else 'FAIL'}")
    print("=" * 78)

    return {
        "state_label": state_label,
        "n_qubits": n_qubits,
        "lambda_output_directory": str(output_dir),
        "full_circuit_resources": full_resources,
        "lambda_tests": lambda_reports,
        "all_pass": all_pass,
    }