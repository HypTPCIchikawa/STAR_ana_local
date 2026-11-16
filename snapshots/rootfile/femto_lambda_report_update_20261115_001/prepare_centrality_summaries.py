#!/usr/bin/env python3
"""Export existing centrality metrics for ROOT plots without recalculating CFs."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import yaml


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text(value):
    return "NA" if value is None else format(value, ".17g")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-json", type=Path, required=True)
    parser.add_argument("--companion-root", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--plot-dir", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.analysis_json.read_text())
    config = data["config"]
    args.work.mkdir(parents=True, exist_ok=True)
    centralities = []
    for cent in config["centralities"]:
        if cent["min"] != cent["max"]:
            continue
        low, high = [float(x) for x in cent["label"].rstrip("%").split("-")]
        centralities.append(dict(cent, percent_low=low, percent_high=high))
    centralities.sort(key=lambda c: c["percent_low"])
    assert len(centralities) == len({c["min"] for c in centralities})
    pair = config["pair_rapidity"]
    merge = pair["adjacent_merge_groups"][1]
    width = (pair["maximum"] - pair["minimum"]) / pair["n_bins"]
    merge_low = pair["minimum"] + (merge[0] - 1) * width
    merge_high = pair["minimum"] + merge[-1] * width
    display = {
        "purpose": "Additional display only; existing CF values/errors are unchanged",
        "datasets": ["EP1", "EP8"],
        "colors": {key: config["display"][key + "_color"] for key in ["EP1", "EP8"]},
        "single_cent9_bins": centralities,
        "map": {
            "tail_color_range": [0.8, 1.2],
            "shape_color_range": [0.85, 1.15],
            "color_steps": 96,
            "palette": "blue below unity, white at unity, orange above unity",
            "missing": "gray cell labelled NA; no numerical replacement",
            "outside_color_range": "saturated endpoint color, hatch, and unchanged numerical value/error",
            "large_error": "! if error >= the full color-range width; display flag only, no exclusion",
            "numerical_label": "central value and one-sigma statistical error on separate lines",
        },
        "metrics": {
            "tail": {"source": "windows.k080_100", "window": config["high_windows"]["k080_100"]},
            "shape": {"source": "shape_ratio", "numerator": config["high_windows"]["k080_100"],
                      "denominator": config["high_windows"]["k050_070"]},
        },
        "normalization": config["normalization"],
        "pair_rapidity": pair,
        "pair_merge02": {"source": "pair_merge02", "bins": merge, "range": [merge_low, merge_high]},
        "he4_rapidity_threshold": config["he4_rapidity_threshold"],
        "uncertainties": data["statistical_assumption"],
        "plot_directory": str(args.plot_dir.resolve()),
    }
    yaml_path = args.work / "summary_display.yaml"
    with yaml_path.open("x") as stream:
        yaml.safe_dump(display, stream, sort_keys=False, allow_unicode=True)
    selections = ["inclusive", "pair_merge02", "he4_below", "he4_above"]
    selections += ["pair_bin%02d" % i for i in range(1, pair["n_bins"] + 1)]
    output = args.work / "summary_metrics.tsv"
    rows = []
    with output.open("x") as stream:
        for dataset in display["datasets"]:
            for cent in centralities:
                for selection in selections:
                    key = "__".join([dataset, cent["id"], selection])
                    curve = data["curves"][key]
                    for metric in ["tail", "shape"]:
                        source = curve["windows"]["k080_100"] if metric == "tail" else curve["shape_ratio"]
                        value, error = source["value"], source["error"]
                        assert (value is None) == (error is None), key
                        if value is not None:
                            assert math.isfinite(value) and math.isfinite(error) and error >= 0, key
                        fields = [dataset, selection, cent["id"], cent["label"], str(cent["min"]),
                                  text(cent["percent_low"]), text(cent["percent_high"]), metric,
                                  text(value), text(error), text(curve["normSE"]), text(curve["normME"]),
                                  text(curve["totalSE"]), text(curve["totalME"])]
                        stream.write("\t".join(fields) + "\n")
                        rows.append(dict(key=key, metric=metric, value=value, error=error))
    parameters = {
        "color.EP1": display["colors"]["EP1"], "color.EP8": display["colors"]["EP8"],
        "map.tail.min": display["map"]["tail_color_range"][0], "map.tail.max": display["map"]["tail_color_range"][1],
        "map.shape.min": display["map"]["shape_color_range"][0], "map.shape.max": display["map"]["shape_color_range"][1],
        "map.color_steps": display["map"]["color_steps"], "cent.count": len(centralities),
        "pair.min": pair["minimum"], "pair.max": pair["maximum"], "pair.bins": pair["n_bins"],
        "merge02.min": merge_low, "merge02.max": merge_high,
        "he4.threshold": config["he4_rapidity_threshold"],
        "tail.low": config["high_windows"]["k080_100"][0], "tail.high": config["high_windows"]["k080_100"][1],
        "reference.low": config["high_windows"]["k050_070"][0], "reference.high": config["high_windows"]["k050_070"][1],
        "norm.low": config["normalization"][0], "norm.high": config["normalization"][1],
        "header": "STAR " + config["production"] + "; #Lambda-^{4}He raw signal CF",
    }
    parameter_path = args.work / "summary_parameters.tsv"
    with parameter_path.open("x") as stream:
        for key, value in parameters.items():
            stream.write(key + "\t" + (value if isinstance(value, str) else text(value)) + "\n")
    provenance = {
        "analysis_json": {"path": str(args.analysis_json.resolve()), "sha256": digest(args.analysis_json)},
        "companion_cf_root": {"path": str(args.companion_root.resolve()), "sha256": digest(args.companion_root),
                              "role": "companion source; summary metrics are read from analysis_json"},
        "display_yaml": {"path": str(yaml_path.resolve()), "sha256": digest(yaml_path)},
        "metrics_tsv": {"path": str(output.resolve()), "sha256": digest(output)},
        "parameters_tsv": {"path": str(parameter_path.resolve()), "sha256": digest(parameter_path)},
        "rows": len(rows), "existing_metrics_only": True, "cf_recalculated": False,
        "single_cent9_bins": len(centralities),
        "statistical_assumption": data["statistical_assumption"],
    }
    with (args.work / "plot_input_provenance.json").open("x") as stream:
        json.dump(provenance, stream, indent=2); stream.write("\n")
    print("PREPARE_CENTRALITY_SUMMARIES_OK rows=%d single_cent9=%d" % (len(rows), len(centralities)))


if __name__ == "__main__":
    main()
