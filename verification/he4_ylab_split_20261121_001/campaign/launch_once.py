#!/usr/bin/env python3
"""Check or submit this approved 20-period replay once; never retry or hadd.

Use --check for a read-only preflight. Only the parent's tmux invocation should
use --submit. The dedicated He4-rapidity watcher verifies outputs and creates CF plots
before sending one minimal completion email. Shared tooling is not edited.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import fcntl
import getpass
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
import ast


STUDY = Path(__file__).resolve().parent
ROOT = STUDY.parent.parent
SHIPPED = {"analysis", "config", "lib", "include", "StMaker", "StRoot"}
EXPECTED_PERIODS = 20
MINIMUM_FREE_BYTES = 500000000000


def now():
    return datetime.now(timezone.utc).isoformat()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inside(path, directory):
    try:
        Path(path).resolve().relative_to(Path(directory).resolve())
        return True
    except ValueError:
        return False


def load(path):
    with Path(path).open() as source:
        return json.load(source)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save(name, value):
    path = STUDY / name
    require(path.parent == STUDY and path.suffix == ".json" and not path.is_symlink(), "unsafe state destination")
    temporary = STUDY / ("." + name + "." + str(os.getpid()) + ".tmp")
    with temporary.open("x") as output:
        json.dump(value, output, indent=2, ensure_ascii=False)
        output.write("\n")
        output.flush()
        os.fsync(output.fileno())
    temporary.replace(path)


def exclusive_json(path, value):
    with Path(path).open("x") as output:
        json.dump(value, output, indent=2)
        output.write("\n")
        output.flush()
        os.fsync(output.fileno())


def file_record(path):
    path = Path(path)
    stat = path.stat()
    require(path.is_file() and stat.st_size > 0, "missing/empty file: " + str(path))
    return {"path": str(path), "bytes": stat.st_size, "mtimeNs": stat.st_mtime_ns,
            "inode": stat.st_ino, "device": stat.st_dev}


def check_sources(config):
    path = Path(config["inputList"])
    require(inside(path, STUDY), "input list must be frozen inside this study")
    require(sha(path) == config["inputListSha256"], "frozen input list SHA256 changed")
    paths = []
    for line in path.read_text().splitlines():
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        require(len(fields) == 1, "period input list must contain one absolute path per line")
        item = Path(fields[0])
        require(item.is_absolute() and inside(item, ROOT / "rootfile") and item.suffix == ".root"
                and not item.name.endswith(".picoDst.root"), "input is not an approved local period ROOT")
        paths.append(item.resolve())
    require(len(paths) == EXPECTED_PERIODS and len(set(paths)) == EXPECTED_PERIODS
            and len({path.name for path in paths}) == EXPECTED_PERIODS, "need 20 distinct period ROOTs and basenames")
    require(len(config["sources"]) == EXPECTED_PERIODS, "need 20 pinned source-stat records")
    observed = []
    for path, expected in zip(paths, config["sources"]):
        require(path == Path(expected["path"]).resolve(), "source order/path differs from the frozen list")
        actual = file_record(path)
        for key in ("bytes", "mtimeNs", "inode", "device"):
            require(actual[key] == expected[key], "source stat changed: " + str(path) + " " + key)
        observed.append(actual)
    return paths, observed


def flat_yaml(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*:\s*([^#\s]+)\s*(?:#.*)?$", line)
        if match:
            result[match.group(1)] = match.group(2).strip("\"'")
    return result


def check_pins(config):
    pins = {str(Path(path).resolve()): expected for path, expected in config["fingerprints"].items()}
    required = [Path(__file__).resolve(), Path(config["joblist"]).resolve(),
                ROOT / "job/run/submit.sh", ROOT / "script/watch_job_notify.py",
                ROOT / "analysis/FemtoLambdaRunSupport.h", ROOT / "include/FemtoLambdaReducedTree.h",
                ROOT / "include/FemtoLambdaReducedTreeConfig.h", ROOT / "StMaker/StFemtoMaker/StFemtoMaker.cxx",
                ROOT / "StMaker/StFemtoMaker/StFemtoMaker.h"]
    required += [ROOT / "lib" / name for name in ("libStarAnaConfig.so", "libStRefMultCorr.so",
                  "libKFParticle.so", "libStKfParticleCommon.so", "libStCommon.so", "libStFemtoMaker.so")]
    mainconf = ROOT / config["mainconf"]
    require(inside(mainconf, ROOT / "config/mainconf"), "mainconf outside config/mainconf")
    required.append(mainconf)
    refs = flat_yaml(mainconf)
    for key, value in refs.items():
        if value.endswith((".yaml", ".yml")):
            target = Path(value) if value.startswith("/") else ROOT / "config" / value
            require(inside(target, ROOT / "config"), "mainconf dependency outside config: " + key)
            required.append(target.resolve())
    info = flat_yaml(ROOT / "config" / refs["analysis"])
    require(info["anaName"] == config["anaName"] and info["outputFileStem"] == config["outputFileStem"],
            "analysis-info naming mismatch")
    require(info["jobName"] + "_" + config["label"] == config["submitName"], "analysis-info jobName mismatch")
    require(info["maxEvents"] == "-1", "production must replay all stored events")
    required += [ROOT / "analysis" / (info[key] + ".C") for key in ("baseRunMacro", "baseAnaMacro")]
    reduced = flat_yaml(ROOT / "config" / refs["reducedTree"])
    mixing = flat_yaml(ROOT / "config" / refs["mixing"])
    require(reduced["mode"] == "replay" and reduced["requireOriginalMixing"] == "false"
            and reduced["resetPoolAtSegmentBoundary"] == "false", "reduced-tree policy is not false/false replay")
    require(mixing["bufferSize"] == "100" and mixing["maxMixEvents"] == "100",
            "both mixing depth limits must be 100")
    require(mixing["mixingMode"] == "bufferAll" and mixing["mixBothDirections"] == "true",
            "production must retain full bidirectional buffer mixing")
    require(all(str(path.resolve()) in pins for path in required), "fingerprints omit a required code/config/lib/launcher dependency")
    for path, expected in pins.items():
        item = Path(path)
        require(inside(item, ROOT), "pinned file outside project")
        require(item.suffix != ".root" or inside(item, ROOT / "config"), "do not hash production ROOT data")
        require(re.fullmatch(r"[0-9a-f]{64}", expected) and sha(item) == expected, "pinned file changed: " + path)
    pilot = config["pilotComparison"]
    require(Path(pilot["path"]).resolve() == STUDY / "pilot_comparison.json", "unexpected pilot comparison path")
    require(sha(pilot["path"]) == pilot["sha256"] and load(pilot["path"]).get("passed") is True,
            "pinned pilot comparison did not pass")
    return pins, info


def normalized_xml(element):
    return [element.tag, sorted(element.attrib.items()), (element.text or "").strip(),
            [normalized_xml(child) for child in element]]


def check_xml(config, info):
    path = Path(config["joblist"]).resolve()
    require(path.parent == ROOT / "job/joblist" and path.name == "joblist_" + config["submitName"] + ".xml",
            "joblist path/name does not match submitName")
    text = path.read_text()
    require(not re.search(r"__[A-Z_]+__", text), "unfilled joblist placeholders")
    xml = ET.fromstring(text)
    require(xml.get("name") == config["submitName"] and xml.get("minFilesPerProcess") == "1"
            and xml.get("maxFilesPerProcess") == "1" and xml.get("copyInputLocally") == "true"
            and xml.get("fileListSyntax") == "paths", "XML identity/staging/min-max input mismatch")
    require(float(xml.get("filesPerHour", "0")) == 1.0 / 24.0,
            "XML must declare the reviewed 24-hour-per-period scheduler estimate")
    require(len(xml.findall("input")) == 1 and xml.find("input").get("URL") == "filelist:" + config["inputList"],
            "XML input is not the frozen 20-period list")
    expected_resource = ET.fromstring(config["resourceUsageXml"])
    resource = xml.findall("ResourceUsage")
    require(expected_resource.tag == "ResourceUsage" and len(resource) == 1
            and normalized_xml(resource[0]) == normalized_xml(expected_resource), "XML ResourceUsage differs from approved config")
    outputs = xml.findall("output")
    expected_output = info["scratchSubdir"] + "/" + config["outputFileStem"] + "_$JOBID.root"
    require(len(outputs) == 1 and outputs[0].get("fromScratch") == expected_output
            and outputs[0].get("toURL") == "file:" + config["outputDirectory"].rstrip("/") + "/",
            "XML output path/stem mismatch")
    package = xml.find("SandBox/Package")
    require(xml.find("SandBox").get("installer") == "ZIP" and package is not None
            and re.fullmatch(r"[A-Za-z0-9_]+", package.get("name", "")), "XML sandbox is not a named ZIP")
    require({(item.text or "").strip() for item in package.findall("File")}
            == {"file:./" + name + "/" for name in SHIPPED}, "XML sandbox ship list mismatch")
    command = xml.findtext("command") or ""
    call = command.replace("\\", "")
    require("analysis/" + info["baseRunMacro"] + ".C" in call
            and ',-1,"' + config["mainconf"] + '")' in call,
            "XML runner/all-events/mainconf arguments mismatch")
    require(not re.search(r"\bhadd\b|watch_job_merge|merge_root_files", command), "XML contains a merge command")
    return xml, package.get("name")


def notifier_compatible():
    path = ROOT / "script/watch_job_notify.py"
    tree = ast.parse(path.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "output_summary")
    sample = "FemtoLambda tree replay job=test_0 requested=10000 completed=10000 available=118240 status=0 RealTime=3.4 CpuTime=3.0"
    for node in ast.walk(function):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "findall"
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "re" and node.args):
            pattern = getattr(node.args[0], "s", None)
            if isinstance(pattern, str) and re.findall(pattern, sample) == [("10000", "10000", "0")]:
                return True
    raise RuntimeError("shared notifier does not recognize tree replay completion; request a reviewed change, do not submit")


def quota(config):
    expected = ["/usr/lpp/mmfs/bin/mmlsquota", "-u", getpass.getuser(), "-Y", "rcf.bnl.gov:gpfs:star-pwg"]
    require(config["quotaCommand"] == expected, "unexpected quota command")
    result = subprocess.run(expected, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=45)
    require(result.returncode == 0, "quota command failed; no submission")
    lines = result.stdout.splitlines()
    header = next(line for line in lines if line.startswith("mmlsquota:user:HEADER:"))
    row = next(line for line in lines if line.startswith("mmlsquota:user:0:"))
    fields = dict(zip(header.split(":")[3:], row.split(":")[3:]))
    free = (int(fields["blockQuota"]) - int(fields["blockUsage"]) - int(fields["blockInDoubt"])) * 1024
    require(config["minimumSoftHeadroomBytes"] >= MINIMUM_FREE_BYTES
            and free >= config["minimumSoftHeadroomBytes"], "quota below approved >=500GB reserve")
    return {"checkedAtUTC": now(), "rawFieldsKiB": fields, "softHeadroomAfterInDoubtBytes": free}


def fresh(config, package):
    require(not (STUDY / "submit_attempt.json").exists(), "submit_attempt.json exists; automatic repeat is forbidden")
    require(not (STUDY / "notify").exists(), "notification state already exists; inspect before any launch")
    require(not any((ROOT / "job/run").glob(config["submitName"] + "*.report")), "matching SUMS report already exists")
    require(not any((ROOT / "job/run/runmeta").glob("runmeta_" + config["submitName"] + "_*.json")), "matching runmeta already exists")
    require(not (ROOT / "job/run" / (package + ".zip")).exists()
            and not (ROOT / "job/run" / (package + ".package")).exists(), "matching sandbox already exists")
    output = Path(config["outputDirectory"])
    require(inside(output, ROOT / "rootfile") and not any(output.glob("*.root")), "output directory holds ROOTs or is outside project")


def preflight(config):
    required_keys = {"schemaVersion", "projectRoot", "directory", "submitHost", "python", "anaName", "label",
                     "submitName", "outputFileStem", "mainconf", "joblist", "outputDirectory", "inputList",
                     "inputListSha256", "sources", "fingerprints", "pilotComparison", "quotaCommand",
                     "minimumSoftHeadroomBytes", "resourceUsageXml", "nativeResources", "email", "emailContent",
                     "pollSeconds", "settleSeconds"}
    require(required_keys <= set(config), "missing config keys: " + ", ".join(sorted(required_keys - set(config))))
    require(config["schemaVersion"] == 1 and Path(config["projectRoot"]).resolve() == ROOT
            and Path(config["directory"]).resolve() == STUDY, "wrong config/project/study")
    require(STUDY.name == "lambda_he4_ylab_split_20261121_001"
            and socket.getfqdn() == config["submitHost"], "wrong study/submit host")
    require(re.fullmatch(r"[A-Za-z0-9_]+", config["anaName"]) and config["anaName"] == "auau3p85_anaFemtoLambda_4He_tree_he4YLab8_mix100"
            and config["label"] == "he4YLab820261121"
            and config["submitName"] == config["anaName"] + "_" + config["label"]
            and config["outputFileStem"] == config["anaName"], "campaign naming mismatch")
    require(config["emailContent"] == "minimal" and config["pollSeconds"] >= 30 and config["settleSeconds"] >= 0,
            "invalid completion-only notification settings")
    require(re.fullmatch(r"[A-Za-z0-9.!#$%&'*+/=?^_`{|}~-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", config["email"]), "invalid email")
    require(Path(config["python"]).is_file() and os.access(config["python"], os.X_OK), "notification Python unavailable")
    require(set(config["nativeResources"]) >= {"request_memory", "request_disk"}
            and all(isinstance(value, (int, float)) and value > 0 for value in config["nativeResources"].values()),
            "native resource expectations must include positive request_memory/request_disk")
    paths, records = check_sources(config)
    pins, info = check_pins(config)
    xml, package = check_xml(config, info)
    notifier_compatible()
    fresh(config, package)
    current_quota = quota(config)
    return {"passed": True, "checkedAtUTC": now(), "sourceCount": len(paths), "sources": records,
            "inputListSha256": config["inputListSha256"], "pinnedFileCount": len(pins),
            "package": package, "quota": current_quota, "notifierTreeReplayMarkerSupported": True,
            "automaticResubmit": False, "automaticMerge": False}, paths, pins, xml


def verify_artifact(record, path):
    require(isinstance(record, dict) and Path(record["path"]).resolve() == path.resolve(), "missing/wrong submit artifact: " + str(path))
    actual = file_record(path)
    require(actual["bytes"] == record["sizeBytes"] and sha(path) == record["sha256"], "submit artifact size/SHA mismatch: " + str(path))
    return {"path": str(path), "bytes": actual["bytes"], "sha256": record["sha256"]}


def submitted_artifacts(config, jobid, pins, package):
    stem = config["submitName"] + "_" + jobid
    path = ROOT / "job/run/runmeta" / ("runmeta_" + stem + ".json")
    meta = load(path)
    require(meta["jobid"] == jobid and meta["anaName"] == config["submitName"]
            and meta["jobPrefix"] == config["submitName"] + jobid, "submitted runmeta identity mismatch")
    require(Path(meta["projectRoot"]).resolve() == ROOT and meta["submit"]["host"] == config["submitHost"]
            and meta["submit"]["mainconf"] == config["mainconf"], "submitted runmeta context mismatch")
    result = {"runmeta": {"path": str(path), "bytes": file_record(path)["bytes"], "sha256": sha(path)}}
    for name, record, expected in [
        ("configlog", meta["submit"].get("configSnapshot"), ROOT / "job/run/configlog" / ("config_" + stem + ".txt")),
        ("joblistlog", meta["submit"].get("joblistSnapshot"), ROOT / "job/run/joblistlog" / ("joblist_" + stem + ".xml")),
        ("runtimeBundle", meta["runtime"].get("runtimeBundle"), ROOT / "job/run/runmeta" / ("runtime_bundle_" + stem + ".tar.gz")),
        ("sumsArtifacts", meta["runtime"].get("sumsArtifacts"), ROOT / "job/run/runmeta" / ("sums_artifacts_" + stem + ".tar.gz"))]:
        result[name] = verify_artifact(record, expected)
    require(result["joblistlog"]["sha256"] == pins[str(Path(config["joblist"]).resolve())], "submitted XML differs from pinned ResourceUsage joblist")
    archive = ROOT / "job/run" / (package + ".zip")
    result["runtimeZip"] = {"path": str(archive), "bytes": file_record(archive)["bytes"], "sha256": sha(archive)}
    checked = 0
    with zipfile.ZipFile(archive) as frozen:
        for source, expected in pins.items():
            relative = Path(source).relative_to(ROOT)
            if relative.parts[0] in SHIPPED:
                require(hashlib.sha256(frozen.read(relative.as_posix())).hexdigest() == expected,
                        "actual ZIP member differs from approved fingerprint: " + str(relative))
                checked += 1
    result["runtimeZip"]["pinnedShippedMembersVerified"] = checked
    return meta, path, result


def native_coverage(config, meta, paths):
    pattern = re.compile(re.escape(meta["jobPrefix"]) + r"_(\d+)\.list$")
    found = {}
    for path in (ROOT / "job/run").glob(meta["jobPrefix"] + "_*.list"):
        match = pattern.fullmatch(path.name)
        require(match is not None and match.group(1) == str(int(match.group(1))), "noncanonical worker list")
        index = int(match.group(1))
        require(index not in found, "duplicate worker index")
        lines = [line.split() for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
        require(len(lines) == 1 and len(lines[0]) in (1, 2), "every native worker must have exactly one input")
        if len(lines[0]) == 2:
            require(re.fullmatch(r"\d+", lines[0][1]), "invalid SUMS event-count column")
        found[index] = {"subjob": index, "nativeList": str(path), "input": lines[0][0], "basename": Path(lines[0][0]).name}
    require(set(found) == set(range(EXPECTED_PERIODS)), "native worker lists are not contiguous 0..19")
    ordered = [found[index] for index in range(EXPECTED_PERIODS)]
    require([row["basename"] for row in ordered] == [path.name for path in paths], "native 20-input basename order differs from frozen list")
    require([str(Path(row["input"]).resolve()) for row in ordered] == [str(path) for path in paths], "native 20-input full path order changed")
    check_sources(config)
    return {"passed": True, "expectedJobs": EXPECTED_PERIODS, "actualJobs": len(ordered), "filesPerWorker": 1,
            "inputListSha256": config["inputListSha256"], "orderedInputs": ordered}


def condor_records(config, meta):
    files = sorted((ROOT / "job/run").glob(meta["jobPrefix"] + "*.condor"))
    require(files, "no generated Condor description")
    events = set()
    records = []
    for path in files:
        fields = {}
        for line in path.read_text().splitlines():
            match = re.match(r"\s*([+A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$", line)
            if match:
                fields[match.group(1).lower()] = match.group(2)
        for name, expected in config["nativeResources"].items():
            require(name.lower() in fields and float(fields[name.lower()]) == expected, "native Condor resource mismatch: " + name)
        event = fields.get("log", "").strip('"')
        require(event.startswith("/") and meta["jobid"] in event and "$" not in event, "unresolved event-log path")
        events.add(event)
        records.append({"path": str(path), "resources": {name: fields[name.lower()] for name in config["nativeResources"]}})
    require(len(events) == 1, "more than one Condor event log")
    return events.pop(), records


def queue_receipt(config, meta, event):
    argv = ["condor_q", "-name", config["submitHost"], "-constraint", 'regexp("' + meta["jobid"] + '",Out)',
            "-json", "-attributes", "ClusterId,ProcId,JobStatus,Out,RequestMemory,RequestDisk"]
    result = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=45, check=True)
    ads = json.loads(result.stdout)
    require(isinstance(ads, list), "invalid queue response")
    seen = set()
    for ad in ads:
        match = re.search(re.escape(meta["jobid"]) + r"_(\d+)\.", ad.get("Out", ""))
        require(match is not None and 0 <= int(match.group(1)) < EXPECTED_PERIODS, "queue contains an unexpected worker")
        require(int(match.group(1)) not in seen, "duplicate queued worker index")
        seen.add(int(match.group(1)))
        for native, expected in config["nativeResources"].items():
            attr = {"request_memory": "RequestMemory", "request_disk": "RequestDisk"}.get(native.lower())
            if attr:
                require(ad.get(attr) == expected, "live queue resource mismatch: " + attr)
    if not ads:
        # A fast failure must not be called a successful/live submission. Accept
        # zero live jobs only when saved Condor events prove all 20 terminal.
        sys.dont_write_bytecode = True
        spec = importlib.util.spec_from_file_location("period_completion_readonly", str(ROOT / "script/watch_job_notify.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        submitted, terminal = module.parse_events(Path(event).read_text())
        require(len(submitted) == EXPECTED_PERIODS and set(terminal) == submitted, "no live queue jobs and no complete terminal-event proof")
    return {"checkedAtUTC": now(), "queueJobCount": len(ads),
            "states": dict(Counter(ad["JobStatus"] for ad in ads)), "ads": ads,
            "queueEmptyDoesNotProveAnalysisSuccess": True}


def start_monitor(config, meta, runmeta, event, state):
    directory = STUDY / "notify"
    require(not directory.exists(), "notification state unexpectedly exists")
    command = [config["python"], "-u", str(ROOT / "script/watch_job_notify.py"), "--runmeta", str(runmeta),
               "--event-log", event, "--email", config["email"], "--email-content", "minimal",
               "--state-dir", str(directory), "--poll-sec", str(config["pollSeconds"]),
               "--settle-sec", str(config["settleSeconds"]), "--watch"]
    state.update(phase="monitor_starting", monitorCommand=command)
    save("launch_state.json", state)
    with (STUDY / "completion_email_watch.log").open("x") as output:
        monitor = subprocess.Popen(command, cwd=str(ROOT), stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
    state.update(monitorPid=monitor.pid)
    save("launch_state.json", state)
    deadline = time.monotonic() + 45
    while not (directory / "settings.json").is_file():
        require(monitor.poll() is None, "completion watcher exited before recording settings")
        require(time.monotonic() < deadline, "completion watcher settings timeout; do not auto-restart")
        time.sleep(1)
    settings = load(directory / "settings.json")
    require(settings["expectedSubjobs"] == EXPECTED_PERIODS and settings["pid"] == monitor.pid
            and settings["host"] == config["submitHost"] and settings["email"] == config["email"]
            and settings["emailContent"] == "minimal" and settings["runmeta"] == str(runmeta)
            and settings["eventLog"] == event, "completion watcher settings mismatch")
    require(monitor.poll() is None, "completion watcher is not live")
    os.kill(monitor.pid, 0)
    state.update(phase="monitor_started", monitorStartedAtUTC=now(), actualJobs=EXPECTED_PERIODS)
    save("launch_state.json", state)
    return monitor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    require(args.config.resolve().parent == STUDY, "launcher config must reside in this study")
    config = load(args.config)
    if args.check:
        receipt, _, _, _ = preflight(config)
        print(json.dumps(receipt, indent=2), flush=True)
        return 0
    lock = (STUDY / "launcher.lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    require(not (STUDY / "submit_attempt.json").exists(), "prior submission attempt; no automatic repeat")
    state = {"phase": "preflight", "pid": os.getpid(), "host": socket.getfqdn(), "startedAtUTC": now(),
             "anaName": config["anaName"], "submitName": config["submitName"],
             "automaticResubmit": False, "automaticMerge": False}
    try:
        receipt, paths, pins, xml = preflight(config)
        save("launch_preflight.json", receipt)
        # Refresh source stats/quota immediately before the durable one-shot gate.
        check_sources(config)
        receipt["quota"] = quota(config)
        attempt = {"attemptedAtUTC": now(), "config": str(args.config.resolve()), "configSha256": sha(args.config),
                   "joblist": config["joblist"], "joblistSha256": sha(config["joblist"]),
                   "inputListSha256": config["inputListSha256"], "quota": receipt["quota"]}
        exclusive_json(STUDY / "submit_attempt.json", attempt)
        state.update(phase="submitting", submitStartedAtUTC=now())
        save("launch_state.json", state)
        relative_xml = os.path.relpath(config["joblist"], str(ROOT / "job/run"))
        with (STUDY / "submit.log").open("x") as output:
            result = subprocess.run(["/bin/bash", str(ROOT / "job/run/submit.sh"), relative_xml],
                                    cwd=str(ROOT / "job/run"), stdout=output, stderr=subprocess.STDOUT)
        text = (STUDY / "submit.log").read_text(errors="replace")
        ids = set(re.findall(re.escape(config["submitName"]) + r"([0-9A-Fa-f]{32})\.report", text))
        state.update(submitReturnCode=result.returncode, submitReturnedAtUTC=now(), possibleReportIds=sorted(ids))
        save("launch_state.json", state)
        require(result.returncode == 0 and len(ids) == 1, "submission nonzero/ambiguous SUMS ID; no automatic retry")
        jobid = ids.pop()
        state["jobid"] = jobid
        meta, runmeta, artifacts = submitted_artifacts(config, jobid, pins, receipt["package"])
        coverage = native_coverage(config, meta, paths)
        event, descriptions = condor_records(config, meta)
        queue = queue_receipt(config, meta, event)
        save("input_coverage.json", coverage)
        save("submission_artifacts.json", artifacts)
        save("native_resources.json", {"descriptions": descriptions, "queue": queue})
        state.update(phase="submitted_verified", runmeta=str(runmeta), eventLog=event,
                     actualJobs=EXPECTED_PERIODS, sourceCount=EXPECTED_PERIODS)
        save("launch_state.json", state)
        command=[config["python"],"-u",str(STUDY/"watch_he4_y.py"),"--launch-config",str(args.config.resolve()),
                 "--launch-state",str(STUDY/"launch_state.json"),"--settings",str(STUDY/"he4_y_watch_settings.json"),
                 "--state-dir",str(STUDY/"completion_watch_001")]
        with (STUDY/"completion_watch.log").open("x") as output:
            monitor=subprocess.Popen(command,cwd=str(ROOT),stdout=output,stderr=subprocess.STDOUT,start_new_session=True)
        save("watch_started.json",dict(pid=monitor.pid,command=command,startedUTC=now(),jobid=jobid))
        print(json.dumps(dict(phase="submitted_verified",jobid=jobid,actualJobs=20,watchPid=monitor.pid,automaticResubmit=False)),flush=True)
        code=monitor.wait()
        require(code==0 and (STUDY/"completion_watch_001/sent.json").is_file(),"Watcher stopped before confirmed completion/mail; no automatic restart")
        sent=load(STUDY/"completion_watch_001/sent.json")
        require(sent.get("mailStatus")=="accepted_by_local_sendmail" and sent.get("kind")=="success","Watcher did not record successful comparison mail")
        state.update(phase="notification_completed",completedAtUTC=now(),watcherReturnCode=code)
        save("launch_state.json",state)
        return 0
    except Exception as error:
        state.update(phase="attention_required", error=str(error), failedAtUTC=now())
        save("launch_state.json", state)
        print(json.dumps({"phase": "attention_required", "error": str(error),
                          "automaticResubmit": False, "automaticMerge": False}), flush=True)
        return 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(json.dumps({"phase": "refused", "error": str(error)}), file=sys.stderr, flush=True)
        sys.exit(1)
