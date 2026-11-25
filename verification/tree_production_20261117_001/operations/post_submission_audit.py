#!/usr/bin/env python3
"""Observe this single submission and audit its provenance and mail-only watcher.

Only post_submission_audit.json is written. This program never submits jobs,
changes the queue, merges outputs, sends mail, or opens production ROOT data.
"""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import socket
import sys
import time
import xml.etree.ElementTree as ET


CONTROL = Path(__file__).resolve().parent
ROOT = CONTROL.parent.parent
STATE_DIR = CONTROL / "tree"
JOBID = "9037A6F05FAD08FC8E9BD9286C8A38A5"
ANA = "auau3p85_anaFemtoLambda_4He_tree_treeWrite20261117"
MAINCONF = "config/mainconf/main_auau3p85_anaFemtoLambda_4He_tree.yaml"
HOST = "starsub03.sdcc.bnl.gov"
EMAIL = "yudai.ichikawa.d3@tohoku.ac.jp"
JOBS, INPUTS = 4058, 73054
INPUT_SHA256 = "b469b4a926dd4781256d3e289c08655c53dd71f41f343411967df342070fd9c5"
POLL_SECONDS, TIMEOUT_SECONDS = 30, 3 * 60 * 60
REPORT = CONTROL / "post_submission_audit.json"
RUNMETA = ROOT / "job/run/runmeta" / ("runmeta_" + ANA + "_" + JOBID + ".json")


class Pending(Exception):
    """A monitor-startup artifact has not appeared yet."""


def now():
    return datetime.now(timezone.utc).isoformat()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    with Path(path).open() as source:
        value = json.load(source)
    require(isinstance(value, dict), "JSON object required: " + str(path))
    return value


def save_report(value):
    # The destination is fixed; source files and launcher state are never written.
    require(REPORT.parent == CONTROL and REPORT.suffix == ".json", "unsafe audit destination")
    require(not REPORT.is_symlink(), "audit destination is a symlink")
    temporary = CONTROL / (".post_submission_audit." + str(os.getpid()) + ".tmp")
    descriptor = os.open(str(temporary), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        with os.fdopen(descriptor, "w") as output:
            json.dump(value, output, ensure_ascii=False, indent=2)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(str(temporary), str(REPORT))
    finally:
        if temporary.exists():
            temporary.unlink()


def digest_file(path):
    path = Path(path)
    before = path.stat()
    require(path.is_file() and before.st_size > 0, "missing/empty artifact: " + str(path))
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns")
    require(all(getattr(before, f) == getattr(after, f) for f in fields),
            "artifact changed during read: " + str(path))
    return {"path": str(path), "sizeBytes": before.st_size, "sha256": digest.hexdigest()}


def artifact(record, expected):
    require(isinstance(record, dict), "missing artifact metadata: " + str(expected))
    require(Path(record.get("path", "")).resolve() == expected.resolve(),
            "unexpected artifact path: " + str(expected))
    measured = digest_file(expected)
    if record.get("sizeBytes") is not None:
        require(measured["sizeBytes"] == record["sizeBytes"], "artifact size mismatch: " + str(expected))
    recorded_sha = record.get("sha256")
    if recorded_sha is not None:
        require(isinstance(recorded_sha, str) and re.fullmatch(r"[0-9a-f]{64}", recorded_sha),
                "invalid recorded SHA256: " + str(expected))
        require(measured["sha256"] == recorded_sha, "artifact SHA256 mismatch: " + str(expected))
    measured["recordedSha256Matched"] = True if recorded_sha is not None else None
    return measured


def command_option(command, option):
    require(command.count(option) == 1, "missing/duplicate watcher option: " + option)
    index = command.index(option)
    require(index + 1 < len(command), "watcher option has no value: " + option)
    return command[index + 1]


def check_watcher(state, config):
    if state.get("monitorPid") is None:
        raise Pending("completion watcher PID not yet recorded")
    settings_path = STATE_DIR / "notify/settings.json"
    if not settings_path.is_file():
        raise Pending("completion watcher settings not yet available")
    settings = read_json(settings_path)
    require(settings.get("expectedSubjobs") == JOBS, "watcher job count mismatch")
    require(settings.get("host") == HOST and settings.get("email") == EMAIL,
            "watcher host/email mismatch")
    require(settings.get("emailContent") == "minimal", "watcher email must be minimal")
    require(Path(settings.get("runmeta", "")).resolve() == RUNMETA.resolve(), "watcher runmeta mismatch")
    require(Path(settings.get("stateDirectory", "")).resolve() == (STATE_DIR / "notify").resolve(),
            "watcher state directory mismatch")
    require(settings.get("eventLog") == state.get("eventLog") and JOBID in settings.get("eventLog", ""),
            "watcher event log mismatch")
    require(settings.get("pollSeconds") == config["pollSeconds"]
            and settings.get("settleSeconds") == config["settleSeconds"], "watcher timing mismatch")
    pid = settings.get("pid")
    require(isinstance(pid, int) and pid > 1 and pid == state.get("monitorPid"), "watcher PID mismatch")
    try:
        os.kill(pid, 0)
        process = Path("/proc") / str(pid)
        proc_state = (process / "stat").read_text().rsplit(")", 1)[1].split()[0]
        require(proc_state not in ("Z", "X"), "watcher is a zombie/exited process")
        command = [part.decode() for part in (process / "cmdline").read_bytes().split(b"\0") if part]
    except OSError as error:
        raise ValueError("watcher is not live/readable: " + str(error))
    watcher = str(ROOT / "script/watch_job_notify.py")
    require(command.count(watcher) == 1 and command.count("--watch") == 1,
            "PID is not the completion-only watcher")
    require(Path(command_option(command, "--runmeta")).resolve() == RUNMETA.resolve(),
            "live watcher runmeta mismatch")
    require(command_option(command, "--email") == EMAIL
            and command_option(command, "--email-content") == "minimal", "live watcher email mismatch")
    require(command_option(command, "--event-log") == settings["eventLog"], "live watcher event log mismatch")
    require(Path(command_option(command, "--state-dir")).resolve() == (STATE_DIR / "notify").resolve(),
            "live watcher state directory mismatch")
    require(not any("hadd" in word.lower() or "watch_job_merge" in word for word in command),
            "merge command found in watcher")
    return {"pid": pid, "live": True, "host": HOST, "expectedSubjobs": JOBS,
            "emailContent": "minimal", "script": watcher, "eventLog": settings["eventLog"],
            "automaticMerge": False}


def audit(state, config):
    require(state.get("jobid") == JOBID and state.get("actualJobs") == JOBS
            and state.get("inputFiles") == INPUTS, "launcher job ID/count mismatch")
    require(state.get("host") == HOST and state.get("anaName") == ANA
            and state.get("mainconf") == MAINCONF, "launcher identity mismatch")
    require(state.get("automaticMerge") is False and state.get("automaticResubmit") is False,
            "launcher enables merge/resubmission")
    require(state.get("submitReturnCode") == 0, "launcher submission did not return successfully")
    require(Path(state.get("runmeta", "")).resolve() == RUNMETA.resolve(), "launcher runmeta path mismatch")
    watcher = check_watcher(state, config)
    proof_path = STATE_DIR / "input_coverage.json"
    if not RUNMETA.is_file() or not proof_path.is_file():
        raise Pending("runmeta/input coverage not yet available")
    meta, proof = read_json(RUNMETA), read_json(proof_path)
    require(meta.get("schemaVersion") == 1 and meta.get("jobid") == JOBID
            and meta.get("anaName") == ANA and meta.get("jobPrefix") == ANA + JOBID,
            "runmeta identity/schema mismatch")
    require(Path(meta.get("projectRoot", "")).resolve() == ROOT, "runmeta project root mismatch")
    require(meta["submit"]["host"] == HOST and meta["submit"]["mainconf"] == MAINCONF,
            "runmeta host/mainconf mismatch")
    require(proof.get("workerJobs") == JOBS and proof.get("inputFiles") == INPUTS,
            "saved input coverage count mismatch")
    pattern = re.compile(re.escape(ANA + JOBID) + r"_(\d+)\.list$")
    found = {}
    for path in (ROOT / "job/run").glob(ANA + JOBID + "_*.list"):
        match = pattern.fullmatch(path.name)
        require(match is not None and path.is_file(), "invalid worker input-list name")
        index = int(match.group(1))
        require(index not in found and match.group(1) == str(index), "duplicate/noncanonical worker index")
        found[index] = str(path)
    require(set(found) == set(range(JOBS)), "worker lists are not contiguous 0..4057")
    frozen = CONTROL / "covered388_pico_files.list"
    require(digest_file(frozen)["sha256"] == INPUT_SHA256, "frozen input list changed")
    launcher = CONTROL / "launch.py"
    require(digest_file(launcher)["sha256"] == config["fingerprints"][str(launcher)],
            "native coverage helper changed")
    # Import definitions only; launch.py's guarded main is not called. Suppress
    # bytecode writes so this observer's only output remains its audit JSON.
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("tree_production_native_coverage", str(launcher))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    entry = config["entries"][0]
    native = module.coverage(ROOT, entry, JOBID, config)
    require(native["workerJobs"] == JOBS and native["inputFiles"] == INPUTS,
            "native coverage count mismatch")
    require(native["lists"] == [found[i] for i in range(JOBS)] == proof.get("lists"),
            "saved/native worker-list order mismatch")
    require(native["orderedGroupSignatures"] == proof.get("orderedGroupSignatures"),
            "saved/native worker grouping mismatch")
    require(1 <= native["maxFiles"] <= 25 and native["maxFiles"] == proof.get("maxFiles"),
            "worker exceeds the approved 25 input files")
    prefix = ANA + "_" + JOBID
    artifacts = {"runmeta": digest_file(RUNMETA)}
    artifacts["configlog"] = artifact(meta["submit"].get("configSnapshot"),
                                     ROOT / "job/run/configlog" / ("config_" + prefix + ".txt"))
    xml_path = ROOT / "job/run/joblistlog" / ("joblist_" + prefix + ".xml")
    artifacts["joblistlog"] = artifact(meta["submit"].get("joblistSnapshot"), xml_path)
    artifacts["runtimeBundle"] = artifact(meta["runtime"].get("runtimeBundle"),
                                         ROOT / "job/run/runmeta" / ("runtime_bundle_" + prefix + ".tar.gz"))
    artifacts["sumsArtifacts"] = artifact(meta["runtime"].get("sumsArtifacts"),
                                         ROOT / "job/run/runmeta" / ("sums_artifacts_" + prefix + ".tar.gz"))
    xml = ET.parse(xml_path).getroot()
    require(xml.get("name") == ANA and xml.get("maxFilesPerProcess") == "25"
            and xml.get("copyInputLocally") == "true", "saved joblist identity/staging mismatch")
    require(xml.find("input").get("URL") == "filelist:" + str(frozen), "saved joblist input mismatch")
    require(xml.find("SandBox").get("installer") == "ZIP", "saved joblist is not ZIP-frozen")
    outputs = xml.findall("output")
    expected_stem = "auau3p85_anaFemtoLambda_4He_tree"
    require(len(outputs) == 1 and outputs[0].get("fromScratch")
            == expected_stem + "/" + expected_stem + "_$JOBID.root", "tree and QA output must be one ROOT")
    require(outputs[0].get("toURL") == "file:" + str(ROOT / "rootfile" / expected_stem / "production20261117") + "/",
            "saved joblist output directory mismatch")
    command = xml.findtext("command") or ""
    require("analysis/run_anaFemtoLambda_4He_tree.C" in command and MAINCONF in command
            and "hadd" not in command.lower(), "saved command runner/mainconf/merge mismatch")
    return {"coverage": {"passed": True, "nativeHelper": str(launcher) + ":coverage",
                         "contiguousLists": JOBS, "workerJobs": native["workerJobs"],
                         "inputFiles": native["inputFiles"], "maxFiles": native["maxFiles"],
                         "frozenInputSha256": INPUT_SHA256, "savedGroupingMatches": True},
            "artifacts": artifacts, "watcher": watcher, "automaticMerge": False,
            "automaticResubmit": False, "productionRootDataRead": False}


def main():
    started, deadline = now(), time.monotonic() + TIMEOUT_SECONDS
    base = {"schemaVersion": 1, "startedAt": started, "observerPid": os.getpid(),
            "host": socket.getfqdn(), "jobid": JOBID, "expectedJobs": JOBS, "expectedInputFiles": INPUTS}
    last_message = None
    try:
        require(CONTROL.name == "femto_lambda_tree_production_20261117_001", "unexpected control directory")
        config = read_json(CONTROL / "launch_config.json")
        require(base["host"] == HOST == config.get("submitHost"), "observer is on the wrong submit host")
        require(Path(config.get("projectRoot", "")).resolve() == ROOT
                and Path(config.get("directory", "")).resolve() == CONTROL, "configuration root mismatch")
        require(config.get("expectedFiles") == INPUTS and config.get("email") == EMAIL
                and config.get("emailContent") == "minimal", "configuration count/email mismatch")
        require(Path(config.get("inputList", "")).resolve() == CONTROL / "covered388_pico_files.list",
                "configuration frozen-input path mismatch")
        require(len(config.get("entries", [])) == 1 and config["entries"][0].get("key") == "tree"
                and config["entries"][0].get("anaName") == ANA
                and config["entries"][0].get("mainconf") == MAINCONF
                and Path(config["entries"][0].get("stateDirectory", "")).resolve() == STATE_DIR,
                "configuration is not the approved single tree entry")
        while True:
            state_path = STATE_DIR / "launch_state.json"
            state = read_json(state_path) if state_path.is_file() else {}
            phase = state.get("phase", "state_not_available")
            base["launchPhase"] = phase
            require(phase != "attention_required", "launcher attention_required: " + str(state.get("error", "unspecified")))
            require(phase != "notification_completed", "notification already completed before live-watcher audit; inspect saved mail records")
            reason = "waiting for submission and completion watcher"
            if phase in ("monitor_starting", "monitor_started"):
                try:
                    result = audit(state, config)
                except Pending as pending:
                    reason = str(pending)
                else:
                    record = dict(base, phase="passed", checkedAt=now(), **result)
                    save_report(record)
                    print(json.dumps({"phase": "passed", "jobid": JOBID, "workerJobs": JOBS,
                                      "inputFiles": INPUTS, "report": str(REPORT)}), flush=True)
                    return 0
            require(time.monotonic() < deadline, "three-hour observer timeout: " + reason)
            record = dict(base, phase="waiting", checkedAt=now(), reason=reason)
            save_report(record)
            message = (phase, reason)
            if message != last_message:
                print(json.dumps({"phase": "waiting", "launchPhase": phase, "reason": reason}), flush=True)
                last_message = message
            time.sleep(POLL_SECONDS)
    except Exception as error:
        record = dict(base, phase="attention_required", checkedAt=now(), error=str(error)[:1500],
                      automaticResubmit=False, automaticMerge=False, emailSentByObserver=False)
        save_report(record)
        print(json.dumps({"phase": "attention_required", "jobid": JOBID,
                          "error": record["error"], "report": str(REPORT)}), flush=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
