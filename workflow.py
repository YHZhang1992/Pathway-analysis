#!/usr/bin/env python3
"""Dependency-free reference runner for the curated analysis workflow packages.

The portable profile proves the complete I/O and audit path on synthetic data.
For production inference, use the package's advanced implementation and approved
statistical environment documented in docs/ENVIRONMENT.md.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import random
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


PALETTE = {"reference": "#1F77B4", "active": "#D62728", "positive": "#009E73", "warning": "#D55E00", "neutral": "#9AA5AA"}


def run_time_utc() -> str:
    """Return an auditable run time, honoring SOURCE_DATE_EPOCH for fixtures."""
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    instant = datetime.fromtimestamp(int(epoch), timezone.utc) if epoch else datetime.now(timezone.utc)
    return instant.isoformat()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Input not found: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("Input must contain at least one data row")
    return rows


def require_columns(rows: list[dict[str, str]], required: Iterable[str]) -> None:
    missing = sorted(set(required) - set(rows[0]))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def number(value: str, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric; received {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("Refusing to write an empty result table")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else math.nan


def sd(values: list[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


def normal_two_sided(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def welch(a: list[float], b: list[float]) -> tuple[float, float]:
    estimate = mean(b) - mean(a)
    se2 = (sd(a) ** 2 / len(a) if a else 0) + (sd(b) ** 2 / len(b) if b else 0)
    p = normal_two_sided(estimate / math.sqrt(se2)) if se2 > 0 else (1.0 if estimate == 0 else 0.0)
    return estimate, p


def bh_adjust(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [1.0] * len(p_values)
    running = 1.0
    size = len(p_values)
    for rank_index in range(size - 1, -1, -1):
        index = order[rank_index]
        running = min(running, p_values[index] * size / (rank_index + 1))
        adjusted[index] = min(1.0, running)
    return adjusted


def grouped_difference(rows: list[dict[str, str]], feature: str, group: str, value: str, reference: str, active: str, transform=lambda x: x) -> list[dict[str, object]]:
    require_columns(rows, (feature, group, value))
    buckets: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        buckets[row[feature]][row[group]].append(transform(number(row[value], value)))
    results = []
    for name, groups in sorted(buckets.items()):
        a, b = groups.get(reference, []), groups.get(active, [])
        if not a or not b:
            continue
        estimate, p = welch(a, b)
        results.append({feature: name, "reference": reference, "active": active, "n_reference": len(a), "n_active": len(b), "mean_reference": round(mean(a), 6), "mean_active": round(mean(b), 6), "estimate": round(estimate, 6), "p_value": p})
    if not results:
        raise ValueError(f"No {feature} has observations in both {reference!r} and {active!r}")
    adjusted = bh_adjust([float(row["p_value"]) for row in results])
    for row, value_adjusted in zip(results, adjusted):
        row["adjusted_p_value"] = value_adjusted
    return results


def correlation(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or len(x) < 3:
        return math.nan
    mx, my = mean(x), mean(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return numerator / denominator if denominator else math.nan


def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> float:
    n = a + b + c + d
    row1, col1 = a + b, a + c
    low, high = max(0, row1 - (n - col1)), min(row1, col1)
    def probability(x: int) -> float:
        return math.comb(col1, x) * math.comb(n - col1, row1 - x) / math.comb(n, row1)
    observed = probability(a)
    return min(1.0, sum(probability(x) for x in range(low, high + 1) if probability(x) <= observed + 1e-12))


def logistic_fit(x: list[list[float]], y: list[int], iterations: int = 2500, rate: float = 0.05) -> list[float]:
    if not x or len(set(y)) < 2:
        raise ValueError("Logistic model requires rows and both outcome classes")
    beta = [0.0] * (len(x[0]) + 1)
    scale = max(1, len(y))
    for _ in range(iterations):
        gradient = [0.0] * len(beta)
        for features, outcome in zip(x, y):
            values = [1.0] + features
            linear = max(-30.0, min(30.0, sum(coef * value for coef, value in zip(beta, values))))
            probability = 1.0 / (1.0 + math.exp(-linear))
            for index, value in enumerate(values):
                gradient[index] += (probability - outcome) * value
        for index in range(len(beta)):
            beta[index] -= rate * gradient[index] / scale
    return beta


def logistic_predict(x: list[list[float]], beta: list[float]) -> list[float]:
    output = []
    for features in x:
        linear = max(-30.0, min(30.0, beta[0] + sum(coef * value for coef, value in zip(beta[1:], features))))
        output.append(1.0 / (1.0 + math.exp(-linear)))
    return output


def auc(y: list[int], probability: list[float]) -> float:
    positives = [p for p, outcome in zip(probability, y) if outcome == 1]
    negatives = [p for p, outcome in zip(probability, y) if outcome == 0]
    if not positives or not negatives:
        return math.nan
    wins = sum(1 if pos > neg else 0.5 if pos == neg else 0 for pos in positives for neg in negatives)
    return wins / (len(positives) * len(negatives))


def svg_bar(path: Path, rows: list[dict[str, object]], label: str, value: str, title: str) -> None:
    values = [float(row[value]) for row in rows if math.isfinite(float(row[value]))]
    if not values:
        return
    width, height, margin = 900, max(300, 90 + len(rows) * 34), 75
    bound = max(abs(v) for v in values) or 1
    center, usable = width / 2, width / 2 - margin - 20
    items = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', f'<text x="{margin}" y="34" font-family="Liberation Sans,Arial" font-size="18" font-weight="bold">{title}</text>', f'<line x1="{center}" x2="{center}" y1="55" y2="{height-30}" stroke="#666666" stroke-width="1"/>']
    for index, row in enumerate(rows):
        y = 72 + index * 34
        current = float(row[value])
        bar_width = abs(current) / bound * usable
        x = center if current >= 0 else center - bar_width
        color = PALETTE["active"] if current >= 0 else PALETTE["reference"]
        items.extend([f'<text x="{margin}" y="{y+14}" font-family="Liberation Sans,Arial" font-size="11">{str(row[label])[:32]}</text>', f'<rect x="{x}" y="{y}" width="{bar_width}" height="18" fill="{color}" opacity="0.88"/>', f'<text x="{x+bar_width+5 if current>=0 else x-45}" y="{y+14}" font-family="Liberation Sans,Arial" font-size="10">{current:.3g}</text>'])
    items.append('</svg>')
    path.write_text("\n".join(items), encoding="utf-8")


def run_profile(profile: str, rows: list[dict[str, str]], config: dict, output: Path) -> list[dict[str, object]]:
    reference, active = config.get("reference", "Control"), config.get("active", "Treatment")
    if profile == "oligo":
        require_columns(rows, ("candidate_id", "sequence", "type"))
        results = []
        for row in rows:
            sequence = row["sequence"].upper().replace("U", "T")
            if set(sequence) - set("ACGT"):
                raise ValueError(f"Invalid nucleotide in {row['candidate_id']}")
            gc = 100 * sum(base in "GC" for base in sequence) / len(sequence)
            run_penalty = 2 if any(base * 4 in sequence for base in "ACGT") else 0
            score = max(0.0, 10 - abs(gc - 45) / 5 - run_penalty)
            results.append({"candidate_id": row["candidate_id"], "type": row["type"], "sequence": sequence, "length": len(sequence), "gc_percent": round(gc, 2), "score": round(score, 3), "status": "review" if score >= 6 else "deprioritize"})
        return sorted(results, key=lambda row: float(row["score"]), reverse=True)
    if profile == "bulk_transcriptomics":
        return grouped_difference(rows, "gene", "group", "count", reference, active, lambda value: math.log2(value + 1))
    if profile == "proteomics":
        return grouped_difference(rows, "analyte", "treatment", "value", reference, active)
    if profile == "single_cell":
        require_columns(rows, ("cell_id", "sample_id", "condition", "n_umi", "n_gene", "mito_percent"))
        results = []
        for row in rows:
            umi, genes, mito = number(row["n_umi"], "n_umi"), number(row["n_gene"], "n_gene"), number(row["mito_percent"], "mito_percent")
            passed = umi >= config.get("min_umi", 500) and genes >= config.get("min_genes", 200) and mito <= config.get("max_mito", 20)
            cluster = "high_complexity" if genes >= config.get("high_complexity_genes", 1000) else "standard"
            results.append({"cell_id": row["cell_id"], "sample_id": row["sample_id"], "condition": row["condition"], "qc_pass": int(passed), "qc_reason": "pass" if passed else "threshold_failure", "portable_cluster": cluster if passed else "not_assigned"})
        return results
    if profile == "perturb_seq":
        return grouped_difference(rows, "target_gene", "condition", "count", reference, active, lambda value: math.log2(value + 1))
    if profile == "multi_omics":
        require_columns(rows, ("sample_id", "rna_feature", "rna_value", "protein_feature", "protein_value"))
        pairs: dict[tuple[str, str], tuple[list[float], list[float]]] = defaultdict(lambda: ([], []))
        for row in rows:
            key = (row["rna_feature"], row["protein_feature"])
            pairs[key][0].append(number(row["rna_value"], "rna_value")); pairs[key][1].append(number(row["protein_value"], "protein_value"))
        return [{"rna_feature": key[0], "protein_feature": key[1], "n": len(values[0]), "pearson_correlation": round(correlation(*values), 6)} for key, values in sorted(pairs.items())]
    if profile == "network_analysis":
        require_columns(rows, ("source", "target", "weight", "interaction_type"))
        nodes: dict[str, dict[str, object]] = defaultdict(lambda: {"neighbors": set(), "weighted_degree": 0.0, "interaction_types": set()})
        seen_edges = set()
        for row in rows:
            source, target = row["source"].strip(), row["target"].strip()
            if not source or not target or source == target: raise ValueError("Network edges require two different non-empty nodes")
            edge = tuple(sorted((source, target)))
            if edge in seen_edges: raise ValueError(f"Duplicate undirected edge: {edge}")
            seen_edges.add(edge); weight = number(row["weight"], "weight")
            if weight < 0: raise ValueError("Network weights must be non-negative")
            for node, neighbor in ((source, target), (target, source)):
                nodes[node]["neighbors"].add(neighbor); nodes[node]["weighted_degree"] += weight; nodes[node]["interaction_types"].add(row["interaction_type"])
        ranked = sorted(nodes, key=lambda node: (-float(nodes[node]["weighted_degree"]), node))
        rank = {node: index + 1 for index, node in enumerate(ranked)}
        return [{"node": node, "degree": len(nodes[node]["neighbors"]), "weighted_degree": round(float(nodes[node]["weighted_degree"]), 6), "hub_rank": rank[node], "interaction_types": ";".join(sorted(nodes[node]["interaction_types"]))} for node in ranked]
    if profile == "pathway_analysis":
        require_columns(rows, ("pathway", "feature", "statistic"))
        pathways: dict[str, list[tuple[str, float]]] = defaultdict(list)
        seen = set()
        for row in rows:
            key = (row["pathway"], row["feature"])
            if key in seen: raise ValueError(f"Duplicate pathway-feature pair: {key}")
            seen.add(key); pathways[row["pathway"]].append((row["feature"], number(row["statistic"], "statistic")))
        output = []
        for pathway, values in sorted(pathways.items()):
            score = mean([value for _, value in values]); top = max(values, key=lambda item: abs(item[1]))
            output.append({"pathway": pathway, "n_features": len(values), "enrichment_score": round(score, 6), "direction": "positive" if score > 0 else "negative" if score < 0 else "neutral", "top_feature": top[0], "top_statistic": top[1]})
        return sorted(output, key=lambda row: abs(float(row["enrichment_score"])), reverse=True)
    if profile == "biomarker":
        require_columns(rows, ("subject_id", "treatment", "visit", "biomarker", "value"))
        buckets: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        for row in rows: buckets[(row["treatment"], row["visit"], row["biomarker"])].append(number(row["value"], "value"))
        return [{"treatment": key[0], "visit": key[1], "biomarker": key[2], "n": len(values), "mean": round(mean(values), 6), "standard_error": round(sd(values) / math.sqrt(len(values)), 6)} for key, values in sorted(buckets.items())]
    if profile == "mmrm_screen":
        require_columns(rows, ("subject_id", "treatment", "visit", "baseline", "outcome"))
        changes = [dict(row, change=str(number(row["outcome"], "outcome") - number(row["baseline"], "baseline"))) for row in rows]
        result = grouped_difference(changes, "visit", "treatment", "change", reference, active)
        for row in result: row["method"] = "portable change-score screen; use advanced/GxP MMRM for inference"
        return result
    if profile == "glmm_screen":
        require_columns(rows, ("subject_id", "treatment", "visit", "baseline", "response"))
        treatments = sorted({row["treatment"] for row in rows}); visits = sorted({row["visit"] for row in rows})
        x = [[number(row["baseline"], "baseline"), float(row["treatment"] == active), float(row["visit"] == visits[-1])] for row in rows]
        y = [int(number(row["response"], "response")) for row in rows]
        beta = logistic_fit(x, y)
        names = ["(Intercept)", "baseline", f"treatment[{active}]", f"visit[{visits[-1]}]"]
        return [{"term": name, "coefficient": round(value, 6), "odds_ratio": round(math.exp(max(-20, min(20, value))), 6), "method": "portable logistic screen; use advanced/GxP GLMM for repeated inference"} for name, value in zip(names, beta)]
    if profile == "responder":
        require_columns(rows, ("subject_id", "treatment", "responder"))
        counts = {group: [0, 0] for group in (reference, active)}
        for row in rows:
            if row["treatment"] in counts: counts[row["treatment"]][int(number(row["responder"], "responder"))] += 1
        ref, act = counts[reference], counts[active]
        if min(sum(ref), sum(act)) == 0: raise ValueError("Both treatment groups are required")
        ref_rate, act_rate = ref[1] / sum(ref), act[1] / sum(act)
        odds = ((act[1] + .5) * (ref[0] + .5)) / ((act[0] + .5) * (ref[1] + .5))
        return [{"reference": reference, "active": active, "n_reference": sum(ref), "n_active": sum(act), "response_rate_reference": ref_rate, "response_rate_active": act_rate, "risk_difference": act_rate - ref_rate, "odds_ratio_haldane": odds, "p_value_fisher_exact": fisher_exact_2x2(act[1], act[0], ref[1], ref[0])}]
    if profile == "bql":
        require_columns(rows, ("subject_id", "visit", "parameter", "result", "lloq"))
        results = []
        seen = set()
        for row in rows:
            key = (row["subject_id"], row["visit"], row["parameter"])
            if key in seen: raise ValueError(f"Duplicate analysis key: {key}")
            seen.add(key); value, lloq = number(row["result"], "result"), number(row["lloq"], "lloq"); flag = value < lloq
            results.append({"subject_id": key[0], "visit": key[1], "parameter": key[2], "original_value": value, "lloq": lloq, "analysis_value": lloq / 2 if flag else value, "bql_flag": int(flag), "bql_rule": "LLOQ/2"})
        return results
    if profile == "survival":
        require_columns(rows, ("subject_id", "group", "time", "event"))
        output = []
        for group in sorted({row["group"] for row in rows}):
            members = [(number(row["time"], "time"), int(number(row["event"], "event"))) for row in rows if row["group"] == group]
            survival = 1.0
            for current in sorted({time for time, event in members if event == 1}):
                at_risk = sum(time >= current for time, _ in members); events = sum(time == current and event == 1 for time, event in members)
                survival *= 1 - events / at_risk
                output.append({"group": group, "time": current, "at_risk": at_risk, "events": events, "survival_probability": survival})
        return output
    if profile in {"feature_selection", "nested_cv"}:
        require_columns(rows, ("sample_id", "outcome")); features = sorted(name for name in rows[0] if name.startswith("feature_")); y = [int(number(row["outcome"], "outcome")) for row in rows]
        if profile == "feature_selection":
            result = []
            for feature in features:
                zeros = [number(row[feature], feature) for row, outcome in zip(rows, y) if outcome == 0]; ones = [number(row[feature], feature) for row, outcome in zip(rows, y) if outcome == 1]
                estimate, p = welch(zeros, ones); result.append({"feature": feature, "mean_difference": estimate, "p_value": p})
            adjusted = bh_adjust([row["p_value"] for row in result])
            for row, adj in zip(result, adjusted): row["adjusted_p_value"] = adj
            return sorted(result, key=lambda row: abs(float(row["mean_difference"])), reverse=True)
        folds = [index % 5 for index in range(len(rows))]; probabilities = [math.nan] * len(rows)
        for fold in range(5):
            train = [i for i in range(len(rows)) if folds[i] != fold]; test = [i for i in range(len(rows)) if folds[i] == fold]
            beta = logistic_fit([[number(rows[i][f], f) for f in features] for i in train], [y[i] for i in train])
            predicted = logistic_predict([[number(rows[i][f], f) for f in features] for i in test], beta)
            for index, value in zip(test, predicted): probabilities[index] = value
        return [{"model": "logistic", "folds": 5, "n": len(y), "accuracy": sum(int(p >= .5) == outcome for p, outcome in zip(probabilities, y)) / len(y), "auc": auc(y, probabilities), "log_loss": -mean([outcome * math.log(max(p, 1e-12)) + (1-outcome)*math.log(max(1-p, 1e-12)) for p, outcome in zip(probabilities, y)])}]
    if profile == "frozen_transfer":
        require_columns(rows, ("sample_id", "feature_1", "feature_2")); model = config["model"]; beta = [model["intercept"], model["coefficients"]["feature_1"], model["coefficients"]["feature_2"]]
        probabilities = logistic_predict([[number(row["feature_1"], "feature_1"), number(row["feature_2"], "feature_2")] for row in rows], beta)
        return [{"sample_id": row["sample_id"], "predicted_probability": probability, "predicted_class": int(probability >= model.get("threshold", .5)), "model_version": model["version"]} for row, probability in zip(rows, probabilities)]
    if profile == "external_validation":
        require_columns(rows, ("patient_id", "pasi_base", "pasi_fu", "predicted_score")); results = []
        for row in rows:
            base, follow = number(row["pasi_base"], "pasi_base"), number(row["pasi_fu"], "pasi_fu")
            observed = 100 * (base - follow) / base
            results.append({"patient_id": row["patient_id"], "observed_percent_change": observed, "predicted_score": number(row["predicted_score"], "predicted_score"), "absolute_error": abs(observed - number(row["predicted_score"], "predicted_score"))})
        return results
    if profile == "decision_support":
        require_columns(rows, ("case_id", "question", "context")); results = []
        for row in rows:
            text = f"{row['question']} {row['context']}".lower()
            recommendation = "Use a prespecified GLMM" if "binary" in text and "repeat" in text else "Route to statistical review"
            results.append({"case_id": row["case_id"], "recommendation": recommendation, "review_status": "human review required", "model": "deterministic portable rule set"})
        return results
    raise ValueError(f"Unknown profile: {profile}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("examples/input.csv"))
    parser.add_argument("--config", type=Path, default=Path("config/workflow.json"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    rows = read_csv(args.input)
    args.output.mkdir(parents=True, exist_ok=True)
    results = run_profile(config["profile"], rows, config, args.output)
    result_path = args.output / "results.csv"
    write_csv(result_path, results)
    numeric = next((name for name in ("estimate", "mean_difference", "risk_difference", "score", "pearson_correlation", "weighted_degree", "enrichment_score") if name in results[0]), None)
    label = next((name for name in ("gene", "analyte", "target_gene", "feature", "candidate_id", "rna_feature", "node", "pathway") if name in results[0]), None)
    if numeric and label: svg_bar(args.output / "summary.svg", results[:20], label, numeric, config["title"])
    manifest = {"workflow": config["profile"], "title": config["title"], "run_time_utc": run_time_utc(), "python": sys.version, "platform": platform.platform(), "input": str(args.input), "input_sha256": sha256(args.input), "config": str(args.config), "config_sha256": sha256(args.config), "result": str(result_path), "result_sha256": sha256(result_path), "row_count": len(results), "palette": PALETTE}
    (args.output / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"status": "ok", "rows": len(results), "result": str(result_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
