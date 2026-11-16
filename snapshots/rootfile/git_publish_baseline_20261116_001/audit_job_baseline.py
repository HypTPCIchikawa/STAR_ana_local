#!/usr/bin/env python3
"""Compare current analysis files with a submitted package; never modify them."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess

import yaml


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def compare_tree(repo, package, group):
    records = []
    for baseline in sorted((package / group).rglob("*")):
        if not baseline.is_file():
            continue
        relative = baseline.relative_to(package)
        current = repo / relative
        entry = {"path": str(relative), "package_sha256": digest(baseline)}
        entry["current_sha256"] = digest(current) if current.is_file() else None
        entry["status"] = "equal" if entry["package_sha256"] == entry["current_sha256"] else (
            "changed" if current.is_file() else "missing_current")
        records.append(entry)
    current_only = []
    for current in sorted((repo / group).rglob("*")):
        relative = current.relative_to(repo)
        if current.is_file() and not (package / relative).exists():
            current_only.append({"path": str(relative), "current_sha256": digest(current)})
    return {"package_files": len(records),
            "equal": sum(r["status"] == "equal" for r in records),
            "changed": sum(r["status"] == "changed" for r in records),
            "missing_current": sum(r["status"] == "missing_current" for r in records),
            "current_only_count": len(current_only), "files": records, "current_only": current_only}


def scalar_strings(node, visited=None):
    if visited is None:
        visited = set()
    if isinstance(node, (dict, list)):
        if id(node) in visited:
            return
        visited.add(id(node))
    if isinstance(node, dict):
        for value in node.values():
            yield from scalar_strings(value, visited)
    elif isinstance(node, list):
        for value in node:
            yield from scalar_strings(value, visited)
    elif isinstance(node, str):
        yield node


def config_closure(tree, repo, main):
    seen = set()
    tree = tree.resolve()
    todo = [(tree / main).resolve()]
    records = []
    while todo:
        path = todo.pop()
        if path in seen:
            continue
        seen.add(path)
        assert len(seen) <= 2000, "Unexpectedly large/cyclic configuration reference closure"
        if not path.is_file():
            records.append({"path": str(path.relative_to(tree)), "missing": True})
            continue
        records.append({"path": str(path.relative_to(tree)), "sha256": digest(path)})
        if path.suffix.lower() not in [".yaml", ".yml"]:
            continue
        document = yaml.safe_load(path.read_text())
        for value in scalar_strings(document):
            if "\n" in value or not value:
                continue
            value_path = Path(value)
            if value_path.suffix.lower() not in {".yaml", ".yml", ".root", ".txt", ".csv", ".dat", ".list", ".xml"}:
                continue
            candidates = []
            if value_path.is_absolute():
                try:
                    candidates.append(tree / value_path.relative_to(repo))
                except ValueError:
                    continue
            else:
                candidates = [tree / value_path, tree / "config" / value_path, path.parent / value_path]
            for candidate in candidates:
                # Configuration closure excludes generated ROOT outputs,
                # input data and operational source files outside config/.
                try:
                    candidate.relative_to(tree / "config")
                except ValueError:
                    continue
                candidate = candidate.resolve()
                try:
                    candidate.relative_to(tree / "config")
                except ValueError:
                    continue
                if candidate.is_file():
                    try:
                        candidate.relative_to(tree)
                    except ValueError:
                        continue
                    todo.append(candidate)
                    break
    return sorted(records, key=lambda r: r["path"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo, package = args.repo.resolve(), args.package.resolve()
    assert package.is_dir()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    groups = {}
    for group in ["StMaker", "include", "analysis", "config"]:
        groups[group] = compare_tree(repo, package, group)
        print("GROUP_OK", group, groups[group]["package_files"], groups[group]["equal"], flush=True)
        if group == "analysis":
            physics_proof = {"created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                             "repo": str(repo), "submitted_package": str(package),
                             "groups": dict(groups), "source_files_modified": False}
            proof_path = args.output_dir / "physics_package_comparison.json"
            if proof_path.exists():
                assert json.loads(proof_path.read_text())["groups"] == physics_proof["groups"], "Source proof changed during audit"
            else:
                with proof_path.open("x") as stream:
                    json.dump(physics_proof, stream, indent=2); stream.write("\n")
    prefix = "auau3p85_anaFemtoLambda_4He_KFParticle_"
    main_paths = ["config/mainconf/main_" + prefix + suffix + ".yaml" for suffix in
                  ["epd8_pairY10_covered25", "epd1_pairY10_covered25", "epd8_kstarSplit_pairY10_covered25"]]
    closures = {}
    for path in main_paths:
        print("CLOSURE_START", path, flush=True)
        closures[path] = {"current": config_closure(repo, repo, Path(path)),
                          "package": config_closure(package, repo, Path(path)) if (package / path).is_file() else None}
        print("CLOSURE_OK", path, len(closures[path]["current"]), flush=True)
    draft_paths = [
        "config/mainconf/main_" + prefix + "epd8_kstarSplit_pairY10_covered25.yaml",
        "config/analysis/analysis_info_" + prefix + "epd8_kstarSplit_pairY10_covered25.yaml",
        "config/hist/hist_auau3p85_femtoLambda_he4_he3_control_epd8KstarPairY10_20261115_kf.yaml",
        "tests/test_femto_epd_kstar_splits.C", "tests/run_femto_epd_kstar_splits_test.C",
    ]
    draft_work = repo / "rootfile/femto_lambda_ep_kstar_validation_20261115_001"
    if draft_work.is_dir():
        draft_paths.extend(str(path.relative_to(repo)) for path in sorted(draft_work.rglob("*")) if path.is_file())
    draft_paths = sorted(set(draft_paths))
    drafts = [{"path": path, "exists": (repo / path).is_file(),
               "sha256": digest(repo / path) if (repo / path).is_file() else None,
               "in_package": (package / path).exists()} for path in draft_paths]
    references = []
    scan = []
    for group in ["config", "tests", "script", "StMaker", "include", "analysis"]:
        scan.extend((repo / group).rglob("*"))
    if draft_work.is_dir():
        scan.extend(draft_work.rglob("*"))
    tokens = {path: {path, Path(path).name} for path in draft_paths}
    for path in sorted(set(scan)):
        if not path.is_file() or path.stat().st_size > 2_000_000:
            continue
        try:
            lines = path.read_text().splitlines()
        except (UnicodeError, OSError):
            continue
        for index, line in enumerate(lines, 1):
            hits = [target for target, words in tokens.items() if any(word in line for word in words)]
            if hits:
                references.append({"path": str(path.relative_to(repo)), "line": index,
                                   "targets": hits, "text": line})
    result = {"created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "repo": str(repo), "submitted_package": str(package),
              "comparison": "Exact SHA256 of file bytes; git HEAD is not the submitted baseline",
              "closure_scope": "YAML and resource paths within config/; generated output/data paths excluded",
              "groups": groups, "mainconf_closures": closures,
              "draft_files": drafts, "draft_references": references,
              "auditor_sha256": digest(Path(__file__)), "source_config_tests_modified": False}
    output = args.output_dir / "baseline_comparison.json"
    with output.open("x") as stream:
        json.dump(result, stream, indent=2); stream.write("\n")
    command = ["git", "status", "--short", "--", "StMaker", "include", "analysis", "config", "tests"]
    status = subprocess.run(command, cwd=str(repo), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
    with (args.output_dir / "git_status_before_archive.txt").open("x") as stream:
        stream.write(status.stdout)
    for group, comparison in groups.items():
        print(group, {key: comparison[key] for key in ["package_files", "equal", "changed", "missing_current", "current_only_count"]})
    print("DRAFT_FILES", len(drafts))
    for entry in drafts:
        print(entry["path"], entry["exists"], entry["sha256"])
    print("BASELINE_AUDIT_OK", output)


if __name__ == "__main__":
    main()
