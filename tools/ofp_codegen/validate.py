#!/usr/bin/env python3
"""Compile and numerically validate generated C++ repairs on bounded examples.

Run in an environment with SymPy, mpmath and a C++17 compiler:
  python tools/ofp_codegen/validate.py --cxx c++ --output /tmp/ofp-validation
The output directory must be absent or empty. Each command, stdout and stderr
is retained. Sampled agreement is not a whole-interval error certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import platform
import shlex
import struct
import subprocess
import sys
import time
import traceback

import mpmath as mp
import sympy as sp


GENERATOR = Path(__file__).with_name("generate.py").resolve()
FLAGS = ["-std=c++17", "-O2", "-Wall", "-Wextra", "-Wconversion",
         "-Wdouble-promotion", "-Werror", "-fno-fast-math", "-ffp-contract=off"]
FAMILIES = [
    ("cos_x2", "(1-cos(x))/x**2", "0", "-0.01", "0.01"),
    ("x_x2", "(x-1)/(x**2-1)", "1", "0.99", "1.01"),
    ("sin_shifted", "sin(x)", "0.1", "0.09", "0.11"),
]
REJECT = [(expression, 10) for expression in
          ["1/x", "log(x)", "sqrt(x)", "exp(-1/x**2)", "unknown(x)",
           "__import__('math').sin(x)", "x.__class__", "[x for x in [1]][0]"]]
# The requested truncation must not conceal a nonanalytic higher-order term.
REJECT += [("1+x**(3/2)", 0), ("1+Abs(x)**3", 2)]


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")


def run_step(directory, command, timeout=120):
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory / "command.json", {"argv": list(map(str, command)), "cwd": str(directory)})
    started = time.monotonic()
    try:
        result = subprocess.run(list(map(str, command)), cwd=directory, text=True,
                                encoding="utf-8", errors="replace", capture_output=True, timeout=timeout)
        code, stdout, stderr = result.returncode, result.stdout, result.stderr
        timed_out = False
    except subprocess.TimeoutExpired as exc:
        code, timed_out = None, True
        stdout = exc.stdout.decode("utf-8", errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode("utf-8", errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
    except OSError as exc:
        code, stdout, stderr, timed_out = None, "", repr(exc), False
    (directory / "stdout.log").write_text(stdout, encoding="utf-8")
    (directory / "stderr.log").write_text(stderr, encoding="utf-8")
    result = {"command": list(map(str, command)), "directory": str(directory), "exit_code": code,
              "timed_out": timed_out, "elapsed_seconds": round(time.monotonic() - started, 6)}
    write_json(directory / "result.json", result)
    return result, stdout


def generation_command(expression, center, left, right, dtype, name, header, degree=10):
    return [sys.executable, str(GENERATOR), "--expression", expression, "--variable", "x",
            "--center", center, "--degree", str(degree), "--interval", left, right,
            "--dtype", dtype, "--name", name, "--output", str(header)]


def literal_value(literal, dtype):
    value = str(literal).rstrip("fFlL")
    parsed = float.fromhex(value) if "p" in value.lower() else float(value)
    return struct.unpack("f", struct.pack("f", parsed))[0] if dtype == "float" else parsed


def harness(name, dtype, metadata, header_name):
    left, right = metadata["interval_literals"]
    center = metadata["center_literal"]
    return f'''#include "{header_name}"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>
using T = {dtype};
static int calls = 0;
static T original(T) {{ ++calls; return T{{12345}}; }}
static_assert(std::is_same_v<decltype(ofp_codegen::{name}_series(T{{0}})), T>);
static_assert(std::is_same_v<decltype(ofp_codegen::{name}_repaired(T{{0}}, &original)), T>);
int main() {{
  const T left = {left};
  const T right = {right};
  const T center = {center};
  const T inf = std::numeric_limits<T>::infinity();
  const std::vector<T> inside = {{left, right, center,
    std::nextafter(left, right), std::nextafter(right, left),
    std::nextafter(center, left), std::nextafter(center, right),
    (left + center)/T{{2}}, (center + right)/T{{2}}}};
  std::cout << std::hexfloat;
  for (T x : inside) {{
    calls = 0;
    const T value = ofp_codegen::{name}_repaired(x, &original);
    const T direct = ofp_codegen::{name}_series(x);
    std::cout << "inside\\t" << x << '\\t' << value << '\\t' << calls << '\\t' << direct << '\\n';
  }}
  const std::vector<T> outside = {{std::nextafter(left, -inf), std::nextafter(right, inf),
    std::numeric_limits<T>::quiet_NaN(), inf, -inf}};
  for (T x : outside) {{
    calls = 0;
    const T value = ofp_codegen::{name}_repaired(x, &original);
    std::cout << "outside\\t" << x << '\\t' << value << '\\t' << calls << '\\n';
  }}
}}
'''


def oracle(family, x):
    if family == "cos_x2":
        # Algebraically identical identity preserves the reference at smallest
        # subnormals, where a direct 100-digit subtraction 1-cos(x) rounds to 0.
        return mp.mpf("0.5") if x == 0 else 2 * (mp.sin(x / 2) / x) ** 2
    if family == "x_x2":
        return 1 / (x + 1)  # removable extension at x=1, valid on this interval
    return mp.sin(x)


def verify_output(family, dtype, stdout):
    checks, failures = [], []
    tolerance = mp.mpf("2e-6" if dtype == "float" else "1e-14")
    for line in stdout.splitlines():
        fields = line.split("\t")
        if len(fields) not in (4, 5):
            raise ValueError(f"Unexpected harness output: {line!r}")
        kind, input_hex, result_hex, count = fields[:4]
        value = float.fromhex(result_hex)
        passed = int(count) == (0 if kind == "inside" else 1)
        check = {"kind": kind, "input_hex": input_hex, "result_hex": result_hex, "callback_calls": int(count)}
        if kind == "inside":
            x = mp.mpf(float.fromhex(input_hex))  # exact promotion of float/double input
            reference = oracle(family, x)
            absolute_error = abs(mp.mpf(value) - reference)
            relative_error = absolute_error / abs(reference) if reference else absolute_error
            passed = passed and mp.isfinite(value) and relative_error <= tolerance and fields[4] == result_hex
            check.update(oracle=mp.nstr(reference, 100), absolute_error=mp.nstr(absolute_error, 20),
                         relative_error=mp.nstr(relative_error, 20), relative_error_limit=str(tolerance),
                         series_hex=fields[4])
        elif kind == "outside":
            passed = passed and value == 12345.0
        else:
            raise ValueError(f"Unexpected row kind {kind!r}")
        check["passed"] = bool(passed)
        checks.append(check)
        if not passed:
            failures.append(check)
    if sum(c["kind"] == "inside" for c in checks) != 9 or sum(c["kind"] == "outside" for c in checks) != 5:
        raise ValueError("Harness must emit exactly nine inside and five fallback samples")
    return checks, failures


def validate_case(output, compiler, family, expression, center, left, right, dtype):
    name = f"{family}_{dtype}"
    directory = output / name
    directory.mkdir()
    header = directory / (name + ".hpp")
    record = {"name": name, "dtype": dtype, "status": "failed", "steps": []}
    step, _ = run_step(directory / "generate", generation_command(expression, center, left, right, dtype, name, header))
    record["steps"].append(step)
    if step["exit_code"] != 0:
        record["reason"] = "generator did not complete successfully"
        return record
    metadata = json.loads(header.with_suffix(".json").read_text(encoding="utf-8"))
    represented = literal_value(metadata["center_literal"], dtype)
    if Fraction(str(metadata["center_exact"])) != Fraction.from_float(represented):
        raise ValueError(f"{name}: center_exact does not equal the emitted center literal")
    if metadata["dtype"] != dtype or metadata["degree"] != 10:
        raise ValueError(f"{name}: incorrect dtype/degree metadata")
    cpp = directory / "harness.cpp"
    cpp.write_text(harness(name, dtype, metadata, header.name), encoding="utf-8")
    executable = directory / ("harness.exe" if os.name == "nt" else "harness")
    step, _ = run_step(directory / "compile", compiler + FLAGS + ["../harness.cpp", "-o", "../" + executable.name])
    record["steps"].append(step)
    if step["exit_code"] != 0:
        record["reason"] = "generated header/harness failed strict compilation"
        return record
    step, stdout = run_step(directory / "execute", [str(executable)])
    record["steps"].append(step)
    if step["exit_code"] != 0:
        record["reason"] = "compiled harness did not complete successfully"
        return record
    record["checks"], failures = verify_output(family, dtype, stdout)
    record.update(status="passed" if not failures else "failed", generated_header=str(header),
                  metadata=str(header.with_suffix(".json")), failures=failures)
    return record


def validate_header_guards(output, compiler):
    """Two function names differing only in case must coexist in one TU."""
    directory = output / "header-case-collision"
    directory.mkdir()
    record = {"name": "case_sensitive_header_guards", "status": "failed", "steps": []}
    for name, expression, filename in [("foo", "sin(x)", "lowercase_name.hpp"),
                                        ("Foo", "cos(x)", "uppercase_name.hpp")]:
        header = directory / filename
        step, _ = run_step(directory / ("generate-" + filename[:-4]),
                           generation_command(expression, "0", "-0.01", "0.01", "double", name, header))
        record["steps"].append(step)
        if step["exit_code"] != 0:
            record["reason"] = "generator failed while preparing the header coexistence test"
            return record
    cpp = directory / "harness.cpp"
    cpp.write_text('''#include "lowercase_name.hpp"
#include "uppercase_name.hpp"
#include <type_traits>
static_assert(std::is_same_v<decltype(ofp_codegen::foo_series(0.0)), double>);
static_assert(std::is_same_v<decltype(ofp_codegen::Foo_series(0.0)), double>);
int main() {
  const double lower = ofp_codegen::foo_series(0.0);
  const double upper = ofp_codegen::Foo_series(0.0);
  return lower == 0.0 && upper == 1.0 ? 0 : 1;
}
''', encoding="utf-8")
    executable = directory / ("harness.exe" if os.name == "nt" else "harness")
    step, _ = run_step(directory / "compile", compiler + FLAGS + ["../harness.cpp", "-o", "../" + executable.name])
    record["steps"].append(step)
    if step["exit_code"] != 0:
        record["reason"] = "both generated headers could not compile in one translation unit"
        return record
    step, _ = run_step(directory / "execute", [str(executable)])
    record["steps"].append(step)
    record["status"] = "passed" if step["exit_code"] == 0 else "failed"
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("--output must be absent or empty; existing artifacts will not be overwritten")
    output.mkdir(parents=True, exist_ok=True)
    compiler = [args.cxx] if Path(args.cxx).is_file() else shlex.split(args.cxx)
    mp.mp.dps = 100
    report = {"status": "running", "python": sys.version, "platform": platform.platform(),
              "mpmath": mp.__version__, "sympy": sp.__version__, "compiler": compiler, "flags": FLAGS,
              "oracle_precision_decimal_digits": 100, "cases": [], "rejections": [], "integration_checks": [],
              "scope": "Six explicitly selected generated repairs, 54 in-interval samples, 30 fallback checks, foo/Foo header coexistence and LF/CRLF/mixed source-patch pipelines; no whole-domain error guarantee.",
              "oracle_notes": ["Inputs are decoded from compiled C++ hexadecimal output and promoted exactly to mpmath.",
                               "cos_x2 reference uses 2*(sin(x/2)/x)^2 and limit 1/2 at zero to resolve subnormal samples.",
                               "x_x2 reference uses the equivalent 1/(x+1), including the removable limit at one."]}
    started = time.monotonic()
    report["compiler_version_step"], _ = run_step(output / "compiler-version", compiler + ["--version"])
    for family, expression, center, left, right in FAMILIES:
        for dtype in ("double", "float"):
            try:
                row = validate_case(output, compiler, family, expression, center, left, right, dtype)
            except Exception as exc:
                row = {"name": f"{family}_{dtype}", "status": "failed", "reason": repr(exc), "traceback": traceback.format_exc()}
            report["cases"].append(row)
            write_json(output / "run-report.json", report)
    for index, (expression, degree) in enumerate(REJECT):
        directory = output / f"reject-{index + 1}"
        header = directory / "rejected.hpp"
        step, _ = run_step(directory, generation_command(expression, "0", "-0.01", "0.01", "double", "rejected", header, degree))
        passed = step["exit_code"] not in (0, None) and not header.exists() and not header.with_suffix(".json").exists()
        report["rejections"].append({"expression": expression, "degree": degree,
                                     "status": "passed" if passed else "failed", "step": step})
        write_json(output / "run-report.json", report)
    try:
        guard_check = validate_header_guards(output, compiler)
    except Exception as exc:
        guard_check = {"name": "case_sensitive_header_guards", "status": "failed",
                       "reason": repr(exc), "traceback": traceback.format_exc()}
    report["integration_checks"].append(guard_check)
    from pipeline_checks import validate_pipeline
    for newline in ("lf", "crlf", "mixed"):
        try:
            pipeline = validate_pipeline(output, compiler, newline, run_step, FLAGS)
        except Exception as exc:
            pipeline = {"name": f"source_patch_{newline}", "status": "failed",
                        "reason": repr(exc), "traceback": traceback.format_exc()}
        report["integration_checks"].append(pipeline)
        write_json(output / "run-report.json", report)
    report["elapsed_seconds"] = round(time.monotonic() - started, 6)
    report["status"] = "passed" if all(r["status"] == "passed" for r in report["cases"] + report["rejections"] + report["integration_checks"]) else "failed"
    write_json(output / "run-report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output / "run-report.json")}))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
