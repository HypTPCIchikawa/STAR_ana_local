"""Compare explicit CF normalization intervals without refitting mass spectra.

The saved mass-fit extraction is the sole numerical input. All fit validity
gates and the shared Method6 purity errors are retained by analyze_group.
The raw-ratio plateau checks use count-only diagonal errors and are diagnostic;
pair/event correlations and mass-model systematic errors are not available.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import yaml

from method6_statistics import analyze_group


METHODS = ("method3", "method5", "method6")
SIDES = ("SE", "ME")


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def tsv(path, header, rows):
    def cell(value):
        if value is None:
            return "nan"
        return str(value).replace("\t", " ").replace("\n", " ")
    path.write_text("\t".join(header) + "\n" + "".join(
        "\t".join(cell(v) for v in row) + "\n" for row in rows))


def estimate_fields(estimate):
    return (estimate["value"], estimate["error"]) if estimate else (None, None)


def aligned_indices(edges, limits):
    if not all(any(abs(x-edge) < 1.e-9 for edge in edges) for x in limits):
        raise ValueError("Normalization/diagnostic limits must coincide with CF edges")
    return [i for i, (a, b) in enumerate(zip(edges, edges[1:]))
            if a >= limits[0]-1.e-9 and b <= limits[1]+1.e-9]


def raw_cells(group, settings):
    out = {side: [] for side in SIDES}
    for side in SIDES:
        source = group["method5Base"][side]
        for first in range(0, len(source), settings["rebin"]):
            rows = source[first:first+settings["rebin"]]
            out[side].append({"value": sum(sum(row["counts"]) for row in rows),
                              "variance": sum(sum(row["error2"]) for row in rows)})
    return out


def raw_ratio(group, settings, norm):
    edges = group["kstarEdges"]
    cells = raw_cells(group, settings)
    norm_indices = aligned_indices(edges, norm)
    total = {side: sum(cells[side][i]["value"] for i in norm_indices) for side in SIDES}
    out = []
    for i, (low, high) in enumerate(zip(edges, edges[1:])):
        same, mixed = (cells[side][i] for side in SIDES)
        s, m = same["value"], mixed["value"]
        ratio = s/m if m > 0 else None
        variance = (same["variance"]/m**2 + s*s*mixed["variance"]/m**4) if m > 0 else None
        cf = ratio*total["ME"]/total["SE"] if ratio is not None and total["SE"] > 0 else None
        cf_variance = None
        if cf is not None and s > 0 and m > 0:
            cf_variance = 0.0
            for j in range(len(edges)-1):
                selected = 1.0 if j in norm_indices else 0.0
                ds = cf*((1.0/s if i == j else 0.0)-selected/total["SE"])
                dm = cf*(selected/total["ME"]-(1.0/m if i == j else 0.0))
                cf_variance += ds*ds*cells["SE"][j]["variance"] + dm*dm*cells["ME"][j]["variance"]
        out.append(dict(bin=i+1, low=low, high=high, SE=s, ME=m,
                        ratio=ratio, error=math.sqrt(variance) if variance is not None else None,
                        cf=cf, cfError=math.sqrt(cf_variance) if cf_variance is not None else None))
    return out


def plateau(points, limits):
    rows = [p for p in points if p["low"] >= limits[0]-1.e-9 and
            p["high"] <= limits[1]+1.e-9 and p["ratio"] is not None and p["error"] > 0]
    if len(rows) < 3:
        return dict(range=limits, valid=False, reason="fewer than three valid ratio points")
    center = sum(limits)/2
    s = sx = sxx = sy = sxy = 0.0
    for row in rows:
        x = (row["low"]+row["high"])/2-center
        w = 1/row["error"]**2
        s += w; sx += w*x; sxx += w*x*x
        sy += w*row["ratio"]; sxy += w*x*row["ratio"]
    mean = sy/s
    det = s*sxx-sx*sx
    intercept = (sxx*sy-sx*sxy)/det
    slope = (s*sxy-sx*sy)/det
    const_chi2 = sum((row["ratio"]-mean)**2/row["error"]**2 for row in rows)
    linear_chi2 = sum((row["ratio"]-intercept-slope*((row["low"]+row["high"])/2-center))**2/
                      row["error"]**2 for row in rows)
    return dict(range=limits, valid=True, points=len(rows), mean=mean, meanError=math.sqrt(1/s),
                constantChi2=const_chi2, constantNdf=len(rows)-1,
                linearReference=center, linearIntercept=intercept, slope=slope,
                slopeError=math.sqrt(s/det), linearChi2=linear_chi2, linearNdf=len(rows)-2,
                relativeChangeAcrossRange=slope*(limits[1]-limits[0])/intercept,
                errorModel="count-only diagonal delta errors; no pair/event covariance")


def failed_fits(data):
    edges = data["massEdges"]
    fit_low, fit_high = data["settings"]["fitRange"]
    bins = [i for i in range(len(edges)-1) if fit_low <= (edges[i]+edges[i+1])/2 <= fit_high]
    result = []
    for group in data["groups"]:
        for side in SIDES:
            for index, row in enumerate(group["methods"]["method5"][side]):
                if row["full"]["valid"]:
                    continue
                result.append(dict(group=group["id"], side=side, purityBin=index+1,
                    purityRange=group["method5PurityEdges"][index:index+2],
                    slidingWindow=[group["originalKstarEdges"][row["windowOriginalFirstBin"]-1],
                                   group["originalKstarEdges"][row["windowOriginalLastBin"]]],
                    storedMassRange=[edges[0], edges[-1]],
                    storedMassObservations=row["fit"]["observationsInStoredMassRange"],
                    fitMassRange=[fit_low, fit_high], fitMassBinRange=[bins[0]+1, bins[-1]+1],
                    observationsInFitRange=sum(row["massCounts"][i] for i in bins),
                    observationsInCentralFitRange=sum(row["centralMassCounts"][i] for i in bins),
                    fit=row["fit"], full=row["full"], regionReason=row.get("regionReason")))
    return result


def write_selected(out, results):
    dump(out/"methods_comparison_selected.json", {"schemaVersion": 1, "groups": results})
    cf_rows, purity_rows = [], []
    for group in results:
        label = group.get("label", group["id"])
        for method, result in group["methods"].items():
            for point in result["points"]:
                cf_rows.append([method, group["id"], label, point["bin"], point["low"], point["high"],
                                *estimate_fields(point["cf"])])
            purity = result["purity"]
            for side in SIDES:
                for i, value in enumerate(purity[side]):
                    purity_rows.append([method, group["id"], label, side, i+1,
                                        purity["edges"][i], purity["edges"][i+1], *estimate_fields(value)])
    tsv(out/"cf_points.tsv", ["method", "group", "label", "bin", "low", "high", "value", "error"], cf_rows)
    tsv(out/"purity_points.tsv", ["method", "group", "label", "side", "bin", "low", "high", "value", "error"], purity_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", required=True)
    args = parser.parse_args()
    config_path = Path(args.settings).resolve()
    config = yaml.safe_load(config_path.read_text())
    required = ("input", "output_dir", "candidates", "selected", "plot_range", "comparison_range",
                "plateau_ranges", "labels")
    for key in required:
        if key not in config:
            raise ValueError("Missing setting: "+key)
    input_path, out = Path(config["input"]).resolve(), Path(config["output_dir"]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out/"normalization_audit.json").exists():
        raise FileExistsError("Refusing to overwrite a completed normalization audit")
    data = json.loads(input_path.read_text())
    for group in data["groups"]:
        group["label"] = config["labels"][group["id"]]
    selected = [c for c in config["candidates"] if c["id"] == config["selected"]]
    if len(selected) != 1:
        raise ValueError("Selected normalization must identify exactly one candidate")
    candidates = selected+[c for c in config["candidates"] if c["id"] != config["selected"]]
    raw_groups = []
    raw_rows = []
    for group in data["groups"]:
        points = raw_ratio(group, data["settings"], selected[0]["range"])
        diagnostics = [plateau(points, limits) for limits in config["plateau_ranges"]]
        raw_groups.append(dict(id=group["id"], label=group["label"], diagnostics=diagnostics, points=points))
        raw_rows.extend([[group["id"], p["bin"], p["low"], p["high"], p["SE"], p["ME"],
                          p["ratio"], p["error"], p["cf"], p["cfError"]] for p in points])
    tsv(out/"raw_ratio_points.tsv", ["group", "bin", "low", "high", "SE", "ME", "ratio", "error", "cf", "cfError"], raw_rows)
    audit = dict(schemaVersion=1, settings=config, input=str(input_path),
                 inputSha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
                 failedFits=failed_fits(data), rawPlateauDiagnostics=raw_groups, candidates=[],
                 recommendation=dict(selected=config["selected"], status="candidate pending inspection",
                    reason="Adjacent disjoint interval avoids the failed purity cell while retaining the fit gates.",
                    limitations=["Changing normalization does not repair failed CF bins or the mass model.",
                                 "Normalization sensitivity and raw-ratio trends must be reported; a plateau is not assumed."]))
    rows = []
    for candidate in candidates:
        print("Starting normalization "+candidate["id"], flush=True)
        results = [analyze_group(group, data["settings"], candidate["range"],
                                 config["plot_range"], config["comparison_range"]) for group in data["groups"]]
        if candidate["id"] == config["selected"]:
            write_selected(out, results)
            print("Selected CF/purity TSVs saved", flush=True)
        groups = []
        for group in results:
            methods = {}
            for method, result in group["methods"].items():
                methods[method] = {key: result[key] for key in ("normalization", "integrated", "valid_plot_bins", "fit_problems")}
                factor, factor_error = estimate_fields(result["normalization"]["factor"])
                cf, cf_error = estimate_fields(result["integrated"]["cf"])
                rows.append([candidate["id"], group["id"], method, *candidate["range"],
                             result["normalization"]["valid"], factor, factor_error, cf, cf_error])
            groups.append(dict(id=group["id"], label=group["label"], methods=methods))
        audit["candidates"].append(dict(id=candidate["id"], normalization=candidate["range"], groups=groups))
        dump(out/"normalization_audit.partial.json", audit)
        tsv(out/"normalization_audit.tsv", ["candidate", "group", "method", "normLow", "normHigh", "valid", "factor", "factorError", "integratedCF", "integratedCFError"], rows)
        print("Finished normalization "+candidate["id"], flush=True)
    reference = {group["id"]: group for group in audit["candidates"][0]["groups"]}
    sensitivity_rows = []
    for candidate in audit["candidates"]:
        for group in candidate["groups"]:
            for method, result in group["methods"].items():
                factor = result["normalization"]["factor"]
                base = reference[group["id"]]["methods"][method]["normalization"]["factor"]
                change = factor["value"]/base["value"]-1 if factor and base else None
                result["relativeScaleChangeFromSelected"] = change
                sensitivity_rows.append([candidate["id"], group["id"], method, change])
    tsv(out/"normalization_sensitivity.tsv", ["candidate", "group", "method", "relativeScaleChangeFromSelected"], sensitivity_rows)
    dump(out/"normalization_audit.json", audit)
    print("Normalization audit complete", flush=True)


if __name__ == "__main__":
    main()
