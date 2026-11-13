#!/usr/bin/env python3
"""Create Lambda diagnostic histogram profiles from explicit mainconfs."""

import argparse
import copy
import math
from pathlib import Path

import yaml


def load(path):
    with path.open() as stream:
        value = yaml.safe_load(stream)
    if not isinstance(value, dict):
        raise ValueError("Expected a YAML mapping: {}".format(path))
    return value


def config_root(mainconf):
    mainconf = mainconf.resolve()
    if mainconf.parent.name != "mainconf" or mainconf.parent.parent.name != "config":
        raise ValueError("mainconf must be below config/mainconf")
    return mainconf.parent.parent


def reference(root, document, key):
    value = document[key]
    if not isinstance(value, str) or not value:
        raise ValueError("Missing explicit config reference: {}".format(key))
    path = Path(value)
    return path if path.is_absolute() else root / path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-mainconf", required=True, type=Path)
    parser.add_argument("--mainconf", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    base_root, target_root = config_root(args.base_mainconf), config_root(args.mainconf)
    base, target = load(args.base_mainconf), load(args.mainconf)
    maker = load(reference(target_root, target, "maker"))
    rapidity = maker["he4RapiditySplit"]
    geometry = maker["mergingGeometry"]
    if rapidity["enabled"] is not True:
        raise ValueError("This profile requires enabled he4RapiditySplit")
    threshold = float(rapidity["yLabThreshold"])
    if not math.isfinite(threshold):
        raise ValueError("Nonfinite He4 lab-rapidity threshold")
    if geometry["method"] != "starHelixFixedGlobalEta":
        raise ValueError("Unsupported mergingGeometry method")
    radius = float(maker["mergingRadiusMeters"])
    extent = float(geometry["tpcHalfLengthCm"])
    tolerance = float(geometry["radiusToleranceCm"])
    if not (radius > 0 and extent > 0 and 0 < tolerance < 100 * radius):
        raise ValueError("Invalid explicit merging geometry")
    if not isinstance(geometry["includePion"], bool):
        raise ValueError("mergingGeometry.includePion must be boolean")
    pair_rapidity = maker.get("pairRapiditySplit")
    if pair_rapidity is not None:
        if not isinstance(pair_rapidity, dict) or set(pair_rapidity) != {
            "enabled", "frame", "nBins", "min", "max"
        }:
            raise ValueError("pairRapiditySplit requires enabled, frame, nBins, min and max")
        if not isinstance(pair_rapidity["enabled"], bool):
            raise ValueError("pairRapiditySplit.enabled must be boolean")
        if pair_rapidity["frame"] != "lab":
            raise ValueError("pairRapiditySplit requires frame: lab")
        pair_bins = pair_rapidity["nBins"]
        pair_min, pair_max = float(pair_rapidity["min"]), float(pair_rapidity["max"])
        if (not isinstance(pair_bins, int) or isinstance(pair_bins, bool) or pair_bins <= 0
                or not math.isfinite(pair_min) or not math.isfinite(pair_max)
                or not pair_min < pair_max or not math.isfinite(pair_max - pair_min)):
            raise ValueError("Invalid pairRapiditySplit binning")
    femto = load(reference(base_root, base, "femtoHist"))
    nuclear = load(reference(base_root, base, "nuclearHist"))
    fh, nh = femto["histograms"], nuclear["histograms"]
    # The EPD provider fills these event-level diagnostics only when enabled.
    epd_axes = {
        "EpdPsi1": {"nBins": 72, "min": -math.pi, "max": math.pi},
        "EpdMixingBin": {"nBins": 8, "min": -0.5, "max": 7.5},
        "EpdStatus": {"nBins": 16, "min": -0.5, "max": 15.5},
    }
    for name, axis in epd_axes.items():
        if name in femto["axes"]:
            raise ValueError("EPD axis already exists: " + name)
        femto["axes"][name] = axis
    cent_axis = fh["hKstarSEVsCent_lambda_he4"]["yAxis"]
    fh["hEpdPsi1VsCent"] = {
        "type": "TH2D", "xAxis": epd_axes["EpdPsi1"], "yAxis": cent_axis,
        "title": "Valid corrected EPD first-order plane;#Psi_{1} [rad];cent9",
    }
    fh["hEpdMixingBinVsCent"] = {
        "type": "TH2D", "xAxis": epd_axes["EpdMixingBin"], "yAxis": cent_axis,
        "title": "Valid Xialei EPD mixing class;EP bin;cent9",
    }
    fh["hEpdStatus"] = {
        "type": "TH1D", "axis": epd_axes["EpdStatus"],
        "title": "EPD calculation status;Status code;Events",
    }
    added = []
    for mode in ("SE", "ME"):
        for mass in ("", "_signal", "_leftSB", "_rightSB"):
            for side, symbol in (("Below", "<"), ("Above", ">")):
                for cent in (False, True):
                    original = "hKstar{}{}_lambda_he4{}".format(mode, "VsCent" if cent else "", mass)
                    name = original + "_YLab" + side
                    if name in fh:
                        raise ValueError("Diagnostic histogram already exists: " + name)
                    definition = copy.deepcopy(fh[original])
                    title = "{} Lambda-He4{}, He4 y_{{lab}} {} {:g};k* [GeV/c];{}".format(
                        mode, mass, symbol, threshold, "cent9" if cent else "Pairs"
                    )
                    definition["title"] = title
                    fh[name] = definition
                    added.append(name)
    pair_count = 0
    if pair_rapidity is not None and pair_rapidity["enabled"]:
        pair_axis = {"nBins": pair_bins, "min": pair_min, "max": pair_max}
        if "PairYLab" in femto["axes"]:
            raise ValueError("PairYLab axis already exists")
        femto["axes"]["PairYLab"] = pair_axis
        for mode in ("SE", "ME"):
            for mass in ("", "_signal", "_leftSB", "_rightSB"):
                original = "hKstar{}VsCent_lambda_he4{}".format(mode, mass)
                name = "hKstar{}VsPairYVsCent_lambda_he4{}".format(mode, mass)
                if name in fh:
                    raise ValueError("Pair-rapidity histogram already exists: " + name)
                fh[name] = {
                    "type": "TH3D",
                    "xAxis": fh[original]["xAxis"],
                    "yAxis": pair_axis,
                    "zAxis": fh[original]["yAxis"],
                    "title": "{} Lambda-He4{}, lab pair rapidity;k* [GeV/c];y_{{pair,lab}};cent9".format(
                        mode, mass
                    ),
                }
                pair_count += 1
    pion_count = 0
    species = []
    for name in list(nh):
        prefix = "hDphiDeta_proton_"
        if name.startswith(prefix) and not name.endswith("_EventField"):
            species.append(name[len(prefix):])
    for species_key in species:
        for suffix in ("", "_EventField"):
            name = "hDphiDeta_proton_" + species_key + suffix
            nh[name]["title"] = (
                "Signal SE proton-{}, STAR helix r={:g} m, global eta;"
                "#Delta#phi* [rad];#Delta#eta"
            ).format(species_key, radius)
            if geometry["includePion"] is True:
                pion_name = "hDphiDeta_pion_" + species_key + suffix
                if pion_name in nh:
                    raise ValueError("Pion merging histogram already exists: " + pion_name)
                definition = copy.deepcopy(nh[name])
                definition["title"] = (
                    "Signal SE pion-{}, STAR helix r={:g} m, global eta;"
                    "#Delta#phi* [rad];#Delta#eta"
                ).format(species_key, radius)
                nh[pion_name] = definition
                pion_count += 1
    for key, document in (("femtoHist", femto), ("nuclearHist", nuclear)):
        output = reference(target_root, target, key)
        if output.resolve() == reference(base_root, base, key).resolve():
            raise ValueError("Refusing to replace the base histogram profile")
        if output.exists() and not args.overwrite:
            raise FileExistsError(str(output))
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w") as stream:
            stream.write("# Diagnostic profile derived from explicit maker settings.\n")
            yaml.safe_dump(document, stream, sort_keys=False, default_flow_style=False)
        print(output)
    print("Added {} physical He4 rapidity views and {} pion merging maps".format(len(added), pion_count))
    if pair_count:
        print("Added {} physical Lambda-He4 pair-rapidity TH3D maps".format(pair_count))


if __name__ == "__main__":
    main()
