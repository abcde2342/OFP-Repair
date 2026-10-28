"""End-to-end checks of emitted diffs, compiled repairs, and failure behavior."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys

import mpmath as mp


def validate_pipeline(output, compiler, newline, run_step, flags):
    directory = output / ("source-patch-" + newline)
    project = directory / "project"
    project.mkdir(parents=True)
    tool = Path(__file__).resolve().parent
    spec = json.loads((tool / "examples/cos_x2.input.json").read_text(encoding="utf-8"))
    spec["target"]["path"] = "original.cpp"
    original = (tool / "examples/original.cpp").read_text(encoding="utf-8")
    if newline == "crlf":
        original = original.replace("\n", "\r\n")
    elif newline == "mixed":
        original = original.replace("\n", "\r\n", 1).replace("double cos_x2(double x) {\n", "double cos_x2(double x) {\r\n")
    source = project / "original.cpp"
    source.write_bytes(original.encode("utf-8"))
    source_bytes = source.read_bytes()
    spec_path = directory / "input.json"
    spec_path.write_text(json.dumps(spec), encoding="utf-8")
    bundle = directory / "bundle"
    record = {"name": "source_patch_" + newline, "status": "failed", "steps": [], "checks": {}}

    def step(name, command):
        result, stdout = run_step(directory / name, command)
        record["steps"].append(result)
        if result["exit_code"] != 0:
            raise RuntimeError(f"{name} failed; see {directory / name / 'stderr.log'}")
        return stdout

    command = [sys.executable, str(tool / "repair.py"), "--spec", str(spec_path),
               "--project-root", str(project), "--output", str(bundle)]
    step("generate", command)
    checks = record["checks"]
    checks["generation_did_not_modify_source"] = source.read_bytes() == source_bytes
    metadata = json.loads((bundle / "cos_x2.json").read_text(encoding="utf-8"))
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    header_hash = hashlib.sha256((bundle / "cos_x2.hpp").read_bytes()).hexdigest()
    checks["header_hash_matches_actual_bytes"] = header_hash == metadata["header_sha256"] == manifest["header_sha256"]
    result, _ = run_step(directory / "reject-existing-output", command)
    record["steps"].append(result)
    checks["existing_output_rejected"] = result["exit_code"] == 2
    # Git is used only in an isolated project under the requested output folder.
    git = shutil.which("git")
    if git is None:
        raise RuntimeError("git is required to verify and apply the generated unified diff")
    step("git-init", [git, "-C", str(project), "init", "--quiet"])
    step("git-check", [git, "-C", str(project), "-c", "core.autocrlf=false", "apply", "--check", str(bundle / "patch.diff")])
    suffix = ".exe" if os.name == "nt" else ""
    # Relative compiler paths also work under MinGW in non-ASCII project paths.
    step("compile-original", compiler + flags + ["../project/original.cpp", "-o", "../original" + suffix])
    before = step("run-original", [str(directory / ("original" + suffix))])
    step("git-apply", [git, "-C", str(project), "-c", "core.autocrlf=false", "apply", str(bundle / "patch.diff")])
    checks["applied_source_matches_preview"] = source.read_bytes() == (bundle / "preview/original.cpp").read_bytes()
    checks["applied_header_matches_generated"] = (project / "ofp_generated/cos_x2.hpp").read_bytes() == (bundle / "cos_x2.hpp").read_bytes()
    step("compile-patched", compiler + flags + ["../project/original.cpp", "-o", "../patched" + suffix])
    after = step("run-patched", [str(directory / ("patched" + suffix))])
    before_rows = [list(map(float.fromhex, line.split())) for line in before.splitlines()]
    after_rows = [list(map(float.fromhex, line.split())) for line in after.splitlines()]
    if len(before_rows) != 5 or len(after_rows) != 5:
        raise ValueError("Expected five actual C++ outputs before and after applying the patch")
    checks["same_inputs"] = all(a[0] == b[0] for a, b in zip(before_rows, after_rows))
    checks["removable_point"] = math.isnan(before_rows[0][1]) and after_rows[0][1] == 0.5
    checks["outside_fallback"] = before_rows[-1][1].hex() == after_rows[-1][1].hex()
    record["numerical_comparison"] = []
    with mp.workdps(100):
        for (x, old), (_, new) in zip(before_rows, after_rows):
            xx = mp.mpf(x)
            reference = mp.mpf("0.5") if x == 0 else 2 * (mp.sin(xx / 2) / xx) ** 2
            old_error = abs(mp.mpf(old) - reference) / abs(reference)
            new_error = abs(mp.mpf(new) - reference) / abs(reference)
            record["numerical_comparison"].append({"input_hex": x.hex(), "original_hex": old.hex(),
                "repaired_hex": new.hex(), "reference": mp.nstr(reference, 100),
                "original_relative_error": str(old_error), "repaired_relative_error": str(new_error)})
            if x == 1e-8:
                checks["cancellation_improvement"] = bool(old_error > mp.mpf("0.1") and new_error < mp.mpf("1e-14"))
    checks["generated_status_does_not_claim_certification"] = manifest["status"] == "generated_unvalidated_interval"
    # Stale specs and existing evidence must not be silently overwritten.
    stale_project = directory / "stale-project"
    stale_project.mkdir()
    (stale_project / "original.cpp").write_bytes(source_bytes.replace(b"std::cos(x)", b"std::cos(2*x)"))
    stale_command = [sys.executable, str(tool / "repair.py"), "--spec", str(spec_path),
                     "--project-root", str(stale_project), "--output", str(directory / "stale-output")]
    result, _ = run_step(directory / "reject-stale-source", stale_command)
    checks["stale_source_rejected"] = result["exit_code"] == 2 and not (directory / "stale-output/patch.diff").exists()
    record["steps"].append(result)
    checks["rejections_preserve_evidence"] = hashlib.sha256((bundle / "cos_x2.hpp").read_bytes()).hexdigest() == header_hash
    record["status"] = "passed" if checks and all(checks.values()) else "failed"
    return record
