# ANNOTAZIONI DEL CONFRONTO CON IL PROGETTO ORIGINALE
# [MODIFICATO] = codice esistente riscritto o con comportamento diverso.
# [AGGIUNTO] = nuovo codice; [RIMOSSO/SOSTITUITO] = descrizione del codice tolto.
# [RIORGANIZZATO] = logica esistente spostata o resa riutilizzabile.
# Le motivazioni descrivono lo scopo tecnico deducibile dai due sorgenti:
# non sono dichiarazioni dell'autore sulle sue intenzioni storiche.
# Le annotazioni sono solo commenti: il codice della versione aggiornata
# rimane invariato. I riferimenti di riga indicano i file ricevuti, prima
# dell'inserimento di questi commenti.
# [FILE INTERAMENTE AGGIUNTO] Automatizza misure al crescere di N per cinque famiglie di stati, con
# CSV e grafici. Misura la versione aggiornata: non esegue in parallelo il progetto originale e non
# dimostra da solo un miglioramento rispetto ad esso.
# [AGGIUNTO: import] argparse configura le prove da terminale; csv salva dati riutilizzabili; time
# misura intervalli; NumPy genera stati e metriche; Matplotlib produce i grafici.
import argparse
import csv
import math
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from quantum_state_preparation_mod import qspCircuitReduced, qspCircuitReduced_mod
from utils_mod import qspParameters_mod, has_nontrivial_phases


# [AGGIUNTO] Separa stati sparsi complessi, densi complessi, positivi, a segni misti e periodici,
# per osservare come struttura e fasi influenzino le risorse.
CATEGORIES = (
    "sparse_complex",
    "dense_complex",
    "dense_positive",
    "dense_signed",
    "periodic",
)


# [AGGIUNTO] Normalizza ogni vettore generato e rifiuta il vettore nullo, affinche le metriche per
# stati puri siano applicabili.
def normalize(vector):
    vector = np.asarray(vector, dtype=complex)
    norm = np.linalg.norm(vector)
    if norm == 0:
        raise ValueError("Cannot normalize the zero vector.")
    return vector / norm


# [AGGIUNTO] Costruisce gli input delle cinque categorie. Gli sparsi hanno supporto max(2,2**N//3);
# i periodici alternano i segni. Il generatore passato come argomento permette di riprodurre gli
# input casuali con il seed.
def make_state(category: str, n_qubits: int, rng: np.random.Generator):
    size = 2**n_qubits

    if category == "sparse_complex":
        vector = np.zeros(size, dtype=complex)
        support_size = max(2, size // 3)
        support = rng.choice(size, support_size, replace=False)
        vector[support] = rng.normal(size=support_size) + 1j * rng.normal(size=support_size)
        return normalize(vector)

    if category == "dense_complex":
        return normalize(rng.normal(size=size) + 1j * rng.normal(size=size))

    if category == "dense_positive":
        return normalize(np.abs(rng.normal(size=size)))

    if category == "dense_signed":
        return normalize(rng.normal(size=size))

    if category == "periodic":
        pattern = np.array([1.0, -1.0], dtype=complex)
        repetitions = math.ceil(size / len(pattern))
        return normalize(np.tile(pattern, repetitions)[:size])

    raise ValueError(f"Unknown category: {category}")


# [AGGIUNTO] Rimuove dall output la fase ricavata dall overlap prima di misurare l errore per
# componente.
def align_global_phase(reference, candidate):
    overlap = np.vdot(reference, candidate)
    if abs(overlap) <= 1e-15:
        return candidate.copy()
    return candidate / (overlap / abs(overlap))


# [AGGIUNTO] Converte le risorse PennyLane in colonne numeriche: fili, porte, profondita e conteggi
# per tipo. I valori predefiniti zero gestiscono campi assenti ma non costituiscono una misura in
# tal caso.
def resource_values(circuit):
    resources = circuit.getCircuitInfos()["resources"]
    gate_types = dict(getattr(resources, "gate_types", {}) or {})
    return {
        "num_wires": getattr(resources, "num_wires", 0),
        "num_gates": getattr(resources, "num_gates", 0),
        "depth": getattr(resources, "depth", 0),
        "cnot": gate_types.get("CNOT", 0),
        "phase_shift": gate_types.get("PhaseShift", 0),
        "ry": gate_types.get("RY", 0),
        "rz": gate_types.get("RZ", 0),
        "qubit_unitary": gate_types.get("QubitUnitary", 0),
    }


# [AGGIUNTO] Esegue una singola prova con seed noto, M=0 e Lambda finale solo quando selezionata.
# Separa tempo di preprocessing/costruzione circuito e tempo di simulazione; registra fidelita,
# distanza di traccia e massimo errore allineato. Non include i costosi test matriciali di
# validation_qsp e non impone qui una soglia PASS/FAIL.
def run_case(category, n_qubits, seed):
    rng = np.random.default_rng(seed)
    coefficients = make_state(category, n_qubits, rng)

    param_start = time.perf_counter()
    has_phases = has_nontrivial_phases(coefficients)
    middle_alphas, phases_alphas, gamma0, _, _ = qspParameters_mod(
        coefficients,
        n_qubits,
        include_phases=has_phases,
    )
    if has_phases:
        circuit = qspCircuitReduced_mod(
            n_qubits, 0, middle_alphas, phases_alphas, gamma0
        )
    else:
        circuit = qspCircuitReduced(n_qubits, 0, middle_alphas, gamma0)
    param_time = time.perf_counter() - param_start

    execution_start = time.perf_counter()
    output = np.asarray(circuit.computeQuantumState(modulo=False), dtype=complex)
    execution_time = time.perf_counter() - execution_start

    overlap = np.vdot(coefficients, output)
    fidelity = float(abs(overlap) ** 2)
    trace_distance = float(np.sqrt(max(0.0, 1.0 - min(1.0, fidelity))))
    aligned_output = align_global_phase(coefficients, output)
    max_component_error = float(np.max(np.abs(coefficients - aligned_output)))

    row = {
        "category": category,
        "n": n_qubits,
        "seed": seed,
        "parameter_time_s": param_time,
        "execution_time_s": execution_time,
        "fidelity": fidelity,
        "trace_distance": trace_distance,
        "max_component_error": max_component_error,
    }
    row.update(resource_values(circuit))
    return row


# [AGGIUNTO] Salva una riga per ogni prova, incluse dimensione, categoria e seed, per poter
# rianalizzare i dati oltre ai grafici.
def write_csv(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# [AGGIUNTO] Raggruppa per categoria e N e calcola la media della metrica tra prove; non calcola
# intervalli di confidenza.
def aggregate(rows, metric):
    result = {}
    for category in CATEGORIES:
        points = []
        n_values = sorted({row["n"] for row in rows if row["category"] == category})
        for n_qubits in n_values:
            values = [
                float(row[metric])
                for row in rows
                if row["category"] == category and row["n"] == n_qubits
            ]
            points.append((n_qubits, float(np.mean(values))))
        result[category] = points
    return result


# [AGGIUNTO] Disegna le medie al crescere di N. Nei grafici logaritmici limita i valori visualizzati
# a 1e-18 per rendere rappresentabili gli zeri; il CSV conserva i valori originali.
def plot_metric(rows, metric, ylabel, title, output_path: Path, log_y=False):
    data = aggregate(rows, metric)
    fig, ax = plt.subplots(figsize=(8, 5))

    for category, points in data.items():
        if not points:
            continue
        x = [point[0] for point in points]
        y = [point[1] for point in points]
        if log_y:
            y = [max(value, 1e-18) for value in y]
        ax.plot(x, y, marker="o", label=category)

    ax.set_xlabel("Number of qubits n")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()
    if log_y:
        ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# [AGGIUNTO] Configura N minimo/massimo, ripetizioni, seed e cartella; verifica l intervallo, esegue
# tutte le combinazioni e salva CSV e sette grafici. I default sono N=2..6, tre prove, seed iniziale
# 2026.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=2)
    parser.add_argument("--max-n", type=int, default=6)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output-dir", default="results/scaling")
    args = parser.parse_args()

    if args.min_n < 1 or args.max_n < args.min_n or args.trials < 1:
        raise ValueError("Invalid benchmark range.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    counter = 0
    for n_qubits in range(args.min_n, args.max_n + 1):
        for category in CATEGORIES:
            for trial in range(args.trials):
                # [AGGIUNTO] Assegna un seed distinto a ogni prova. La riproducibilita vale a parita
                # di configurazione/ordine; lo stato periodico resta identico fra le prove, che ne
                # ripetono comunque le misure temporali.
                seed = args.seed + counter
                row = run_case(category, n_qubits, seed)
                rows.append(row)
                counter += 1
                print(
                    f"n={n_qubits} category={category} trial={trial + 1}/{args.trials} "
                    f"F={row['fidelity']:.12f} depth={row['depth']} gates={row['num_gates']}"
                )

    csv_path = output_dir / "scaling_metrics.csv"
    write_csv(rows, csv_path)

    plot_metric(
        rows,
        "parameter_time_s",
        "Parameter generation time [s]",
        "Classical parameter-generation cost",
        output_dir / "scaling_parameter_time.png",
        log_y=True,
    )
    plot_metric(
        rows,
        "execution_time_s",
        "Circuit simulation time [s]",
        "Circuit execution cost",
        output_dir / "scaling_execution_time.png",
        log_y=True,
    )
    plot_metric(
        rows,
        "num_gates",
        "Total gates",
        "Total gate count versus n",
        output_dir / "scaling_total_gates.png",
    )
    plot_metric(
        rows,
        "depth",
        "Circuit depth",
        "Circuit depth versus n",
        output_dir / "scaling_depth.png",
    )
    plot_metric(
        rows,
        "cnot",
        "CNOT count",
        "CNOT count versus n",
        output_dir / "scaling_cnot.png",
    )
    plot_metric(
        rows,
        "trace_distance",
        "Trace distance",
        "Preparation error versus n",
        output_dir / "scaling_trace_distance.png",
        log_y=True,
    )
    plot_metric(
        rows,
        "max_component_error",
        "Maximum component error",
        "Component-wise preparation error versus n",
        output_dir / "scaling_component_error.png",
        log_y=True,
    )

    print(f"\nCSV saved in: {csv_path}")
    print(f"Plots saved in: {output_dir}")


# [AGGIUNTO] Rende il benchmark eseguibile da terminale e importabile senza avviare automaticamente
# le prove.
if __name__ == "__main__":
    main()
