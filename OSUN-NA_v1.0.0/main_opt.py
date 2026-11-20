# ANNOTAZIONI DEL CONFRONTO CON IL PROGETTO ORIGINALE
# [MODIFICATO] = codice esistente riscritto o con comportamento diverso.
# [AGGIUNTO] = nuovo codice; [RIMOSSO/SOSTITUITO] = descrizione del codice tolto.
# [RIORGANIZZATO] = logica esistente spostata o resa riutilizzabile.
# Le motivazioni descrivono lo scopo tecnico deducibile dai due sorgenti:
# non sono dichiarazioni dell'autore sulle sue intenzioni storiche.
# Le annotazioni sono solo commenti: il codice della versione aggiornata
# rimane invariato. I riferimenti di riga indicano i file ricevuti, prima
# dell'inserimento di questi commenti.
# [MODIFICATO: struttura del file] Il vecchio script eseguiva input, simulazione e report al livello
# del modulo. Ora li organizza in funzioni e main(), con esecuzione protetta da __name__. Parte
# della logica e solo spostata (ad esempio stati 1-10 e generazione casuale); non va considerata un
# nuovo algoritmo.
# [MODIFICATO: import] Aggiunti re, shutil, Path e redirect_stdout per nomi/cartelle/report;
# aggiunti has_nontrivial_phases e il nuovo validatore. os e rimasto importato anche se non usato.
import os
import re
import sys
import math
import time
import shutil
from pathlib import Path
from contextlib import redirect_stdout

import numpy as np

from quantum_state_preparation_mod import qspCircuitReduced, qspCircuitReduced_mod
from utils_mod import (
    toKet,
    stateToVector,
    qspParameters_mod,
    has_nontrivial_phases,
)
from validation_qsp import validate_current_qsp_run


# [MODIFICATO] Precisione di stampa da 3 a 6 cifre e larghezza da 100 a 140, per rendere piu
# leggibili i dettagli numerici; non modifica i calcoli.
np.set_printoptions(precision=6, suppress=True, linewidth=140)

# [MODIFICATO] Path e mkdir sostituiscono stringa results/ e controllo os.path.exists. ATOL=1e-9
# centralizza la tolleranza assoluta dei nuovi confronti finali.
RESULT_ROOT = Path("results")
ATOL = 1e-9
RESULT_ROOT.mkdir(parents=True, exist_ok=True)


# [AGGIUNTO] Centralizza le etichette e aggiunge due stati periodici (11-12). La selezione Custom
# passa da 11 a 13.
PREDEFINED_STATE_NAMES = {
    1: "Bell Phi+",
    2: "Bell Phi-",
    3: "Bell Psi+",
    4: "Bell Psi-",
    5: "GHZ n=3",
    6: "GHZ n=4",
    7: "W n=3",
    8: "W n=4",
    9: "Dicke n=3",
    10: "Dicke n=4",
    11: "Periodic n=2",
    12: "Periodic n=4",
}


# [AGGIUNTO / SOSTITUITO] Sostituisce map(float,...): accetta coefficienti complessi, i o j e
# virgole come separatori. Permette input manuali con fasi; la normalizzazione viene verificata
# separatamente.
def parse_custom_coefficients(text: str) -> list[complex]:
    values = []
    for token in text.replace(",", " ").split():
        values.append(complex(token.replace("i", "j")))
    return values


# [AGGIUNTO] Controlla che gli elementi siano 2**N e che la norma sia circa 1 prima della
# costruzione. Converte in complesso; non normalizza automaticamente. np.isclose conserva qui la
# tolleranza relativa predefinita.
def check_input_state(coefficients, n_qubits: int) -> np.ndarray:
    vector = np.asarray(coefficients, dtype=complex)
    expected_size = 2**n_qubits

    if vector.size != expected_size:
        raise ValueError(
            f"The state must contain 2**N={expected_size} coefficients; "
            f"received {vector.size}."
        )

    norm = np.linalg.norm(vector)
    if not np.isclose(norm, 1.0, atol=1e-8):
        raise ValueError(f"The input state is not normalized: norm={norm}.")

    return vector


# [AGGIUNTO] Genera lo stato normalizzato a segni alternati +1,-1. Introduce un caso strutturato con
# fasi, oltre ai dieci stati gia presenti.
def periodic_state(n_qubits: int) -> list[float]:
    size = 2**n_qubits
    vector = np.array([1.0 if index % 2 == 0 else -1.0 for index in range(size)])
    vector /= np.linalg.norm(vector)
    return vector.tolist()


# [AGGIUNTO] Trasforma le etichette in nomi di cartella e distingue esplicitamente plus/minus,
# evitando di confondere Bell Phi+ e Bell Phi-.
def state_slug(state_label: str) -> str:
    text = state_label.strip().lower()
    text = text.replace("+", " plus ").replace("-", " minus ").replace("=", " ")
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_") or "state"


# [RIORGANIZZATO / AGGIUNTO] Estrae dal vecchio menu i dieci stati noti, conservandone coefficienti
# e dimensioni. Restituisce anche il nome e aggiunge i casi periodici 11 e 12, cosi la stessa
# selezione serve sia al menu sia alla modalita batch.
def get_predefined_state(selection: int):
    match selection:
        case 1:
            return "Bell Phi+", 2, [math.sqrt(1 / 2), 0, 0, math.sqrt(1 / 2)]
        case 2:
            return "Bell Phi-", 2, [math.sqrt(1 / 2), 0, 0, -math.sqrt(1 / 2)]
        case 3:
            return "Bell Psi+", 2, [0, math.sqrt(1 / 2), math.sqrt(1 / 2), 0]
        case 4:
            return "Bell Psi-", 2, [0, math.sqrt(1 / 2), -math.sqrt(1 / 2), 0]
        case 5:
            return "GHZ n=3", 3, [math.sqrt(1 / 2), 0, 0, 0, 0, 0, 0, math.sqrt(1 / 2)]
        case 6:
            return "GHZ n=4", 4, [
                math.sqrt(1 / 2), 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, math.sqrt(1 / 2),
            ]
        case 7:
            return "W n=3", 3, [
                0, math.sqrt(1 / 3), math.sqrt(1 / 3), 0,
                math.sqrt(1 / 3), 0, 0, 0,
            ]
        case 8:
            return "W n=4", 4, [
                0, math.sqrt(1 / 4), math.sqrt(1 / 4), 0,
                math.sqrt(1 / 4), 0, 0, 0,
                math.sqrt(1 / 4), 0, 0, 0, 0, 0, 0, 0,
            ]
        case 9:
            return "Dicke n=3", 3, [
                0, 0, 0, math.sqrt(1 / 3),
                0, math.sqrt(1 / 3), math.sqrt(1 / 3), 0,
            ]
        case 10:
            return "Dicke n=4", 4, [
                0, 0, 0, math.sqrt(1 / 6),
                0, math.sqrt(1 / 6), math.sqrt(1 / 6), 0,
                0, math.sqrt(1 / 6), math.sqrt(1 / 6), 0,
                math.sqrt(1 / 6), 0, 0, 0,
            ]
        case 11:
            return "Periodic n=2", 2, periodic_state(2)
        case 12:
            return "Periodic n=4", 4, periodic_state(4)
        case _:
            raise ValueError("Invalid predefined-state selection.")


# [AGGIUNTO] Stima la fase globale dall overlap e la rimuove prima del confronto delle componenti;
# evita di segnalare come diversi due vettori fisicamente equivalenti. Se l overlap e quasi zero non
# tenta l allineamento.
def align_global_phase(reference: np.ndarray, candidate: np.ndarray):
    overlap = np.vdot(reference, candidate)
    if abs(overlap) <= 1e-15:
        return candidate.copy(), 1.0 + 0.0j
    phase = overlap / abs(overlap)
    return candidate / phase, phase


# [AGGIUNTO / SOSTITUITO] Amplia la precedente stampa di metriche: confronta ogni componente dopo
# allineamento, calcola fidelita e distanza di traccia dell intero stato puro. Prima printMetrics
# calcolava la distanza tra stati ridotti al solo qubit 0, che puo nascondere differenze globali. Il
# test usa atol con rtol=0.
def print_component_test(target: np.ndarray, output: np.ndarray, n_qubits: int, atol: float = ATOL):
    aligned_output, phase = align_global_phase(target, output)
    differences = np.abs(target - aligned_output)
    overlap = np.vdot(target, output)
    fidelity = float(abs(overlap) ** 2)
    trace_distance = float(np.sqrt(max(0.0, 1.0 - min(1.0, fidelity))))
    component_pass = bool(np.allclose(target, aligned_output, atol=atol, rtol=0.0))

    print("\n" + "=" * 78)
    print("FINAL STATE CHECK")
    print("=" * 78)
    print(f"Target-state norm: {np.linalg.norm(target):.16f}")
    print(f"Output-state norm: {np.linalg.norm(output):.16f}")
    print(f"Raw overlap <target|output>: {overlap}")
    print(f"Global phase found from overlap: {phase}")
    print(f"Fidelity: {fidelity:.16f}")
    print(f"Trace distance: {trace_distance:.16e}")

    print("\nComponent-by-component comparison after global-phase alignment:")
    for index, (expected, actual, difference) in enumerate(zip(target, aligned_output, differences)):
        bit_string = format(index, f"0{n_qubits}b")
        print(
            f"  |{bit_string}>  target={expected!s:>24}  "
            f"output={actual!s:>24}  |diff|={difference:.3e}"
        )

    print(f"\nMaximum component error: {np.max(differences):.3e}")
    print(f"Component test: {'PASS' if component_pass else 'FAIL'}")

    return {
        "fidelity": fidelity,
        "trace_distance": trace_distance,
        "max_component_error": float(np.max(differences)),
        "component_pass": component_pass,
        "global_phase": phase,
        "aligned_output": aligned_output,
    }


# [AGGIUNTO / RIORGANIZZATO] Porta in console un riepilogo leggibile; i dettagli rimangono nel
# report. Prima stdout veniva rediretto globalmente e in console rimaneva soltanto il messaggio
# finale.
def print_terminal_summary(
    coefficients,
    desired_state,
    state_check,
    parameter_time,
    execution_time,
    resources,
    n_wires,
):
    print(f"> Coefficients: {coefficients.tolist()}")
    print(f"> Modulo: {np.abs(coefficients).tolist()}")
    print("_" * 50)
    print("\nDesired quantum state after global-phase correction:")
    print(toKet(desired_state, n_wires))
    print(f"\nFidelity: {state_check['fidelity']}")
    print(f"\nTrace distance: {state_check['trace_distance']}")
    print(f"\nTime taken to generate parameters: {parameter_time}s.")
    print(f"Time taken to execute the circuit: {execution_time}s.")
    print()
    print(f"num_wires: {getattr(resources, 'num_wires', '?')}")
    print(f"num_gates: {getattr(resources, 'num_gates', '?')}")
    print(f"depth: {getattr(resources, 'depth', '?')}")


# [AGGIUNTO] Sposta la figura dalla cartella generale a quella dello stato, verifica che esista e
# sostituisce la precedente destinazione. Tiene insieme immagine e risultati della stessa
# esecuzione.
def move_circuit_image(n_wires: int, output_dir: Path) -> Path:
    source = RESULT_ROOT / f"circuit_{n_wires}.png"
    destination = output_dir / f"circuit_{n_wires}.png"

    if not source.exists():
        raise FileNotFoundError(f"Circuit image not found: {source}")

    if destination.exists():
        destination.unlink()
    shutil.move(str(source), str(destination))
    return destination


# [RIORGANIZZATO / MODIFICATO] Raccoglie in una funzione riutilizzabile l intera esecuzione:
# controllo input, parametri, validazione Lambda, simulazione, confronto finale e report
# strutturato. I nuovi controlli interni sono commentati nei rispettivi blocchi.
def run_qsp_state(state_label: str, n_qubits: int, coefficients) -> dict:
    coefficients = check_input_state(coefficients, n_qubits)
    # [AGGIUNTO] Organizza risultati per stato con sottocartelle validation/output. Una nuova
    # esecuzione della stessa etichetta riscrive il report e le immagini associati.
    run_slug = state_slug(state_label)
    state_dir = RESULT_ROOT / run_slug
    validation_dir = state_dir / "validation"
    output_dir = state_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "qsp_output.txt"

    # [MODIFICATO CENTRALE] Prima M=2*N: il dispositivo aveva 3N qubit. Ora M=0, coerentemente con
    # LambdaNoAncilla, e il dispositivo ne ha N. La simulazione del vettore di stato richiede quindi
    # 2**N ampiezze anziche 2**(3N).
    M = 0

    param_start_time = time.perf_counter()
    # [MODIFICATO] Sostituisce il confronto esatto np.angle(coefficients)==0 con la selezione a
    # tolleranza; include_phases evita anche il preprocessing delle fasi quando non necessario.
    has_phases = has_nontrivial_phases(coefficients)
    middle_alphas, phases_alphas, gamma0, global_phases, lambda_diagonals = qspParameters_mod(
        coefficients,
        n_qubits,
        include_phases=has_phases,
    )

    if has_phases:
        circuit = qspCircuitReduced_mod(n_qubits, M, middle_alphas, phases_alphas, gamma0)
    else:
        circuit = qspCircuitReduced(n_qubits, M, middle_alphas, gamma0)

    # [MODIFICATO: significato della misura] Il tempo include ora anche la costruzione dell oggetto
    # circuito, che nell originale avveniva dopo param_end_time. I due tempi non misurano
    # esattamente lo stesso intervallo.
    parameter_time = time.perf_counter() - param_start_time

    validation_start_time = time.perf_counter()

    # [MODIFICATO] Usa redirect_stdout dentro contesti che chiudono il file e ripristinano la
    # console automaticamente. Sostituisce l assegnazione globale sys.stdout e aggiunge un report
    # completo di input e parametri.
    with open(output_path, "w", encoding="utf-8") as output_file:
        with redirect_stdout(output_file):
            print("=" * 78)
            print("INPUT STATE")
            print("=" * 78)
            print(f"Selected state: {state_label}")
            print(f"Number of qubits: N={n_qubits}")
            print(f"Coefficient count: {len(coefficients)}")
            print(f"Input-state norm: {np.linalg.norm(coefficients):.16f}")
            print(f"Coefficients: {coefficients.tolist()}")
            print(f"Moduli: {np.abs(coefficients).tolist()}")
            print("\nComplete target expansion:")
            print(toKet(coefficients, n_qubits, include_zeros=True))

            print("\n" + "=" * 78)
            print("PARAMETER GENERATION")
            print("=" * 78)
            print(f"Non-trivial phases present: {has_phases}")
            print(f"First RY parameter gamma0: {gamma0:.12f} rad")
            for index, alpha_vector in enumerate(middle_alphas, start=2):
                print(f"Magnitude Lambda n={index}: {np.asarray(alpha_vector)}")
            if has_phases:
                print(f"Final phase Lambda n={n_qubits}: {np.asarray(phases_alphas)}")
            print(
                "Global-phase product C (correction factor): "
                f"{math.prod(global_phases) if global_phases else 1.0}"
            )
            print(f"Parameter-generation time: {parameter_time:.9f} s")

            try:
                # [AGGIUNTO] Valida i Lambda effettivamente presenti prima della simulazione QSP.
                # Gli errori sono registrati e tradotti in un risultato di fallimento.
                validation_summary = validate_current_qsp_run(
                    n_qubits=n_qubits,
                    state_label=state_label,
                    qsp_circuit=circuit,
                    middle_alphas=middle_alphas,
                    phases_alphas=phases_alphas,
                    lambda_diagonals=lambda_diagonals,
                    has_final_phase_lambda=has_phases,
                    output_dir=validation_dir,
                )
            except Exception as error:
                print("\nVALIDATION ERROR")
                print(f"{type(error).__name__}: {error}")
                validation_summary = None

    validation_time = time.perf_counter() - validation_start_time

    # [AGGIUNTO] Interrompe questa esecuzione se la validazione non e completabile o se un Lambda
    # non supera i test, evitando di presentare il resto come riuscito.
    if validation_summary is None:
        return {
            "success": False,
            "state_label": state_label,
            "state_dir": state_dir,
            "message": f"Validation could not be completed. See {output_path} for details.",
        }

    if not validation_summary["all_pass"]:
        return {
            "success": False,
            "state_label": state_label,
            "state_dir": state_dir,
            "message": f"Validation failed. See {output_path} and {validation_dir}.",
        }

    circ_start_time = time.perf_counter()
    # [AGGIUNTO] Elimina la vecchia immagine attesa prima del disegno: se il nuovo disegno fallisce
    # non si riutilizza per errore una figura precedente.
    default_circuit_image = RESULT_ROOT / f"circuit_{n_qubits + M}.png"
    if default_circuit_image.exists():
        default_circuit_image.unlink()

    with open(output_path, "a", encoding="utf-8") as output_file:
        with redirect_stdout(output_file):
            print("\n" + "=" * 78)
            print("QSP EXECUTION")
            print("=" * 78)
            circuit.printCircuit("figure", modulo=False)
            output_state = np.asarray(circuit.computeQuantumState(modulo=False), dtype=complex)

    execution_time = time.perf_counter() - circ_start_time
    circuit_image_path = move_circuit_image(n_qubits + M, output_dir)

    # [RIORGANIZZATO] Conserva la costruzione del target con eventuali ancilla; con M=0 il fattore
    # aggiunto e [1] e il target resta quello di ingresso.
    new_qubits = [0] * M
    state_vector = stateToVector(new_qubits)
    target_coefficients = np.kron(coefficients, state_vector).astype(complex)
    # [MODIFICATO / AGGIUNTO] La correzione globale era gia moltiplicata all output per la stampa.
    # Ora viene conservata e controllata esplicitamente: si verifica C, si predice l output grezzo
    # come target/C e si controlla output*C.
    global_phase_correction = complex(math.prod(global_phases) if global_phases else 1.0 + 0.0j)
    if abs(global_phase_correction) <= 1e-15:
        raise ValueError("The collected global-phase correction cannot be zero.")

    predicted_raw_phase = 1.0 / global_phase_correction
    predicted_raw_output = target_coefficients * predicted_raw_phase
    corrected_output_state = output_state * global_phase_correction

    with open(output_path, "a", encoding="utf-8") as output_file:
        with redirect_stdout(output_file):
            print("\nComplete target-state expansion:")
            print(toKet(target_coefficients, n_qubits + M, include_zeros=True))

            print("\nComplete raw output expansion:")
            print(toKet(output_state, n_qubits + M, include_zeros=True))

            print("\n" + "=" * 78)
            print("GLOBAL PHASE CHECK")
            print("=" * 78)
            print(f"Collected correction factor C: {global_phase_correction}")
            print(f"|C|: {abs(global_phase_correction):.16f}")
            # [AGGIUNTO] Controlla che C abbia modulo 1: una fase globale non deve modificare la
            # norma.
            phase_modulus_pass = bool(
                np.isclose(abs(global_phase_correction), 1.0, atol=ATOL, rtol=0.0)
            )
            print(f"Unit-modulus phase test: {'PASS' if phase_modulus_pass else 'FAIL'}")
            print(f"Predicted raw-output phase C^(-1): {predicted_raw_phase}")

            # [AGGIUNTO] Verifica l output grezzo contro target/C; controlla anche la convenzione
            # del fattore di fase raccolto, non solo l equivalenza fisica.
            raw_prediction_error = float(np.max(np.abs(output_state - predicted_raw_output)))
            raw_prediction_pass = bool(
                np.allclose(output_state, predicted_raw_output, atol=ATOL, rtol=0.0)
            )
            print(f"Maximum raw-output prediction error: {raw_prediction_error:.3e}")
            print(f"Raw-output phase prediction test: {'PASS' if raw_prediction_pass else 'FAIL'}")

            print("\nOutput after global-phase correction:")
            print(toKet(corrected_output_state, n_qubits + M, include_zeros=True))
            # [AGGIUNTO] Confronta il vettore corretto con il target senza un ulteriore
            # allineamento: testa direttamente la correzione calcolata dal preprocessing.
            correction_error = float(np.max(np.abs(target_coefficients - corrected_output_state)))
            correction_pass = bool(
                np.allclose(target_coefficients, corrected_output_state, atol=ATOL, rtol=0.0)
            )
            print(f"Corrected-output norm: {np.linalg.norm(corrected_output_state):.16f}")
            print(f"Maximum corrected-output error: {correction_error:.3e}")
            print(f"Collected global-phase correction test: {'PASS' if correction_pass else 'FAIL'}")

            state_check = print_component_test(
                target_coefficients,
                output_state,
                n_qubits + M,
            )
            # [AGGIUNTO] Confronta la fase prevista con quella misurata dall overlap, come controllo
            # indipendente della coerenza di C.
            overlap_phase = state_check["global_phase"]
            normalized_predicted_phase = predicted_raw_phase / abs(predicted_raw_phase)
            phase_consistency_error = abs(overlap_phase - normalized_predicted_phase)
            phase_consistency_pass = bool(phase_consistency_error <= ATOL)
            print("\nGlobal-phase consistency test:")
            print(f"  Predicted phase: {normalized_predicted_phase}")
            print(f"  Measured phase:  {overlap_phase}")
            print(f"  Error: {phase_consistency_error:.3e}")
            print(f"  Result: {'PASS' if phase_consistency_pass else 'FAIL'}")

            if not state_check["component_pass"]:
                print("\nWARNING: final component test failed.")

            print("\n" + "=" * 78)
            print("TIMING AND CIRCUIT METRICS")
            print("=" * 78)
            print(f"Time taken to generate parameters: {parameter_time:.9f} s")
            print(f"Time taken to validate the current circuit: {validation_time:.9f} s")
            print(f"Time taken to execute the circuit: {execution_time:.9f} s")

            # [MODIFICATO] Esplicita fili, porte, profondita e distribuzioni, invece di stampare
            # soltanto l oggetto resources. Il tempo di validazione e riportato separatamente; il
            # tempo di esecuzione comprende anche disegno e output nel contesto cronometrato.
            specs = circuit.getCircuitInfos()
            resources = specs["resources"]
            print(f"num_wires: {getattr(resources, 'num_wires', '?')}")
            print(f"num_gates: {getattr(resources, 'num_gates', '?')}")
            print(f"depth: {getattr(resources, 'depth', '?')}")
            print(f"gate_types: {dict(getattr(resources, 'gate_types', {}) or {})}")
            print(f"gate_sizes: {dict(getattr(resources, 'gate_sizes', {}) or {})}")
            print(f"Circuit image: {circuit_image_path}")
            print(f"Validation directory: {validation_dir}")


    print_terminal_summary(
        coefficients,
        corrected_output_state,
        state_check,
        parameter_time,
        execution_time,
        resources,
        n_qubits + M,
    )
    print(f"\nResults saved in {state_dir} directory.")

    # [AGGIUNTO] Il successo richiede tutti i controlli Lambda, componenti e fase globale. Il
    # dizionario restituito serve al riepilogo batch e al codice di uscita.
    all_checks_pass = bool(
        validation_summary["all_pass"]
        and state_check["component_pass"]
        and phase_modulus_pass
        and raw_prediction_pass
        and correction_pass
        and phase_consistency_pass
    )

    return {
        "success": all_checks_pass,
        "state_label": state_label,
        "state_dir": state_dir,
        "output_path": output_path,
        "validation_dir": validation_dir,
        "circuit_image": circuit_image_path,
        "message": "PASS" if all_checks_pass else "One or more final checks failed.",
    }


# [AGGIUNTO / RIORGANIZZATO] Genera il menu dalla tabella dei nomi, includendo i nuovi casi e
# Custom=13, invece della lunga stringa incorporata nell input.
def predefined_menu() -> str:
    lines = ["Select a known quantum state or insert a custom one:", ""]
    for index in range(1, 13):
        lines.append(f"{index}. {PREDEFINED_STATE_NAMES[index]}")
    lines.append("13. Custom")
    lines.append("")
    return "\n".join(lines)


# [RIORGANIZZATO / MODIFICATO] Isola la selezione dello stato dal calcolo; utilizza il nuovo parser
# complesso per Custom. Restituisce nome, N e coefficienti al chiamante.
def read_single_state():
    selection = int(input(predefined_menu() + "> ").strip())
    if 1 <= selection <= 12:
        return get_predefined_state(selection)
    if selection == 13:
        n_qubits = int(input("Number of qubits? > ").strip())
        coefficients = parse_custom_coefficients(
            input(
                "Insert 2**N coefficients separated by spaces "
                "(real or complex, e.g. 0.5+0.5j) > "
            ).strip()
        )
        return f"Custom n={n_qubits}", n_qubits, coefficients
    raise ValueError("Invalid input.")


# [RIORGANIZZATO] La generazione densa/sparsa e la normalizzazione erano gia presenti. Ora
# restituisce anche un etichetta e usa ValueError. "Real negative" diventa "real signed": la normale
# genera valori sia positivi sia negativi, come gia avveniva prima.
def read_random_state():
    n_qubits = int(input("Number of qubits? > ").strip())

    match input("Dense [d] or sparse [s] random state? > ").strip().lower():
        case "d":
            density_label = "dense"
            match input(
                "Real positive [p], real signed [n] or complex [c] coefficients? > "
            ).strip().lower():
                case "p":
                    coefficient_label = "real positive"
                    numbers = np.abs(np.random.normal(size=2**n_qubits))
                case "n":
                    coefficient_label = "real signed"
                    numbers = np.random.normal(size=2**n_qubits)
                case "c":
                    coefficient_label = "complex"
                    numbers = np.random.normal(size=2**n_qubits) + 1j * np.random.normal(size=2**n_qubits)
                case _:
                    raise ValueError("Invalid input.")

            numbers /= np.linalg.norm(numbers, ord=2)
            coefficients = numbers.tolist()

        case "s":
            density_label = "sparse"
            match input(
                "Real positive [p], real signed [n] or complex [c] coefficients? > "
            ).strip().lower():
                case "p":
                    coefficient_label = "real positive"
                    numbers = np.abs(np.random.normal(size=2**n_qubits))
                case "n":
                    coefficient_label = "real signed"
                    numbers = np.random.normal(size=2**n_qubits)
                case "c":
                    coefficient_label = "complex"
                    numbers = np.random.normal(size=2**n_qubits) + 1j * np.random.normal(size=2**n_qubits)
                case _:
                    raise ValueError("Invalid input.")

            zero_indices = np.random.choice(
                2**n_qubits,
                np.random.randint(1, 2**n_qubits),
                replace=False,
            )
            numbers[zero_indices] = 0
            numbers /= np.linalg.norm(numbers, ord=2)
            coefficients = numbers.tolist()

        case _:
            raise ValueError("Invalid input.")

    return f"Random {density_label}, {coefficient_label}, n={n_qubits}", n_qubits, coefficients


# [AGGIUNTO] Esegue i dodici casi senza ripetere il menu, raccoglie PASS/FAIL e prosegue anche se un
# caso solleva un eccezione. Consente un controllo ripetibile sui casi noti.
def run_all_predefined_states():
    results = []

    print("\nRunning all predefined states (1-12).")
    for selection in range(1, 13):
        state_label, n_qubits, coefficients = get_predefined_state(selection)
        print("\n" + "#" * 78)
        print(f"[{selection}/12] {state_label}")
        print("#" * 78)

        try:
            result = run_qsp_state(state_label, n_qubits, coefficients)
        except Exception as error:
            result = {
                "success": False,
                "state_label": state_label,
                "state_dir": RESULT_ROOT / state_slug(state_label),
                "message": f"{type(error).__name__}: {error}",
            }
            print(f"Run failed: {result['message']}")

        results.append(result)

    print("\n" + "=" * 78)
    print("ALL PREDEFINED STATES - SUMMARY")
    print("=" * 78)
    for result in results:
        status = "PASS" if result["success"] else "FAIL"
        print(f"{status:4}  {result['state_label']:<20}  {result['state_dir']}")
        if not result["success"]:
            print(f"      {result['message']}")

    passed = sum(1 for result in results if result["success"])
    print(f"\nPassed: {passed}/{len(results)}")
    print(f"Failed: {len(results) - passed}/{len(results)}")

    return all(result["success"] for result in results)


# [RIORGANIZZATO / AGGIUNTO] Coordina le modalita specifica, casuale e la nuova modalita all. Usa
# exit 2 per verifiche fallite ed exit 1 per eccezioni: rende l esito interpretabile anche da
# automazioni.
def main():
    try:
        mode = input(
            "Specific quantum state [s], random vector [r] or ALL predefined states [a]? > "
        ).strip().lower()

        if mode in {"a", "all"}:
            success = run_all_predefined_states()
            if not success:
                sys.exit(2)
            return

        if mode == "s":
            state_label, n_qubits, coefficients = read_single_state()
        elif mode == "r":
            state_label, n_qubits, coefficients = read_random_state()
        else:
            raise ValueError("Invalid input.")

        result = run_qsp_state(state_label, n_qubits, coefficients)
        if not result["success"]:
            print(result["message"])
            sys.exit(2)

    except Exception as error:
        print(error)
        sys.exit(1)


# [AGGIUNTO] Il menu si avvia solo eseguendo il file come programma. Importarlo permette di riusare
# le funzioni senza avviare input e simulazione; la creazione di results resta comunque al livello
# del modulo.
if __name__ == "__main__":
    main()