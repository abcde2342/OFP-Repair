#!/usr/bin/env python3
"""Turn an explicit JSON repair specification into a reviewable source patch.

This command never modifies the target project. Apply the emitted patch.diff
after compiling/testing the preview. Source localization is an explicit input,
not an inference from the mathematical expression.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
from pathlib import Path, PurePosixPath
import posixpath
import re
import sys
import textwrap

from generate import GenerationError, generate


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def fields(value, required, optional=()):
    if not isinstance(value, dict):
        raise GenerationError("Specification and target must be JSON objects.")
    missing = set(required) - value.keys()
    unknown = value.keys() - set(required) - set(optional)
    if missing or unknown:
        raise GenerationError(f"Invalid fields: missing={sorted(missing)}, unknown={sorted(unknown)}")


def portable_path(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_./-]+", value):
        raise GenerationError("Paths must be relative portable paths using / and ASCII letters/digits/_.-.")
    path = PurePosixPath(value)
    if path.is_absolute() or any(p in ("..", ".git") for p in path.parts) or str(path) != value:
        raise GenerationError(f"Unsafe or noncanonical project path: {value!r}")
    return path


def project_path(root, relative):
    path = root.joinpath(*relative.parts)
    if not path.resolve().is_relative_to(root):
        raise GenerationError("Target path escapes the project root through a symbolic link.")
    return path


def load_spec(path):
    spec = json.loads(path.read_text(encoding="utf-8"))
    fields(spec, ("schema_version", "name", "expression", "variable", "center", "degree",
                  "interval", "dtype", "target"))
    if type(spec["schema_version"]) is not int or spec["schema_version"] != 1:
        raise GenerationError("Only schema_version 1 is supported.")
    for key in ("name", "expression", "variable", "center", "dtype"):
        if not isinstance(spec[key], str):
            raise GenerationError(f"{key} must be a string; use strings for exact mathematical constants.")
    if type(spec["degree"]) is not int or not 0 <= spec["degree"] <= 64:
        raise GenerationError("degree must be an integer between 0 and 64.")
    if spec["dtype"] not in ("float", "double"):
        raise GenerationError("dtype must be float or double.")
    if (not isinstance(spec["interval"], list) or len(spec["interval"]) != 2
            or not all(isinstance(x, str) for x in spec["interval"])):
        raise GenerationError("interval must contain exactly two mathematical strings.")
    fields(spec["target"], ("path", "function", "definition"))
    if not all(isinstance(v, str) for v in spec["target"].values()):
        raise GenerationError("All target fields must be strings.")
    for key, value in (("name", spec["name"]), ("variable", spec["variable"]),
                       ("target.function", spec["target"]["function"])):
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", value):
            raise GenerationError(f"{key} must be an ASCII identifier starting with a letter.")
    return spec


def prepare_source(spec, root):
    target = spec["target"]
    relative = portable_path(target["path"])
    source = project_path(root, relative)
    if source.suffix not in (".cpp", ".cc", ".cxx"):
        raise GenerationError("Target must be a C++ source file (.cpp/.cc/.cxx).")
    raw = source.read_bytes()
    old = raw.decode("utf-8")
    if old.startswith("\ufeff") or not old.endswith("\n"):
        raise GenerationError("Target must be UTF-8 without BOM and end with a newline.")
    definition = target["definition"].replace("\r\n", "\n").strip("\n")
    dtype, variable, function = spec["dtype"], spec["variable"], target["function"]
    pattern = (rf"({dtype}\s+{re.escape(function)}\s*\(\s*{dtype}\s+"
               rf"{re.escape(variable)}\s*\)\s*\{{)(.*)\}}")
    match = re.fullmatch(pattern, definition, re.DOTALL)
    if match is None:
        raise GenerationError("target.definition must be a complete one-argument float/double free-function definition matching dtype, variable and target.function.")
    # Require an exact, unique, whole-line match. No heuristic source matching.
    lines = r"\r?\n".join(re.escape(line) for line in definition.split("\n"))
    anchored = re.compile(r"(?m)^" + lines + r"(?=\r?\n)")
    matches = list(anchored.finditer(old))
    if len(matches) != 1:
        raise GenerationError(f"Expected exactly one matching target definition; found {len(matches)}. Source may have changed.")
    body = textwrap.indent(textwrap.dedent(match[2]).strip("\n"), "        ")
    replacement = (match[1] + "\n"
                   f"    return ofp_codegen::{spec['name']}_repaired({variable},\n"
                   f"        +[]({dtype} {variable}) -> {dtype} {{\n"
                   + body + "\n        });\n}")
    found = matches[0]
    matched_source = found.group()
    newline = "\r\n" if "\r\n" in matched_source else "\n"
    edited = old[:found.start()] + replacement.replace("\n", newline) + old[found.end():]
    header_relative = portable_path(f"ofp_generated/{spec['name']}.hpp")
    if project_path(root, header_relative).exists():
        raise GenerationError(f"Generated header already exists: {header_relative}; choose a fresh name.")
    include = posixpath.relpath(str(header_relative), str(relative.parent))
    first_newline = old.index("\n")
    include_newline = "\r\n" if first_newline > 0 and old[first_newline - 1] == "\r" else "\n"
    new = f'#include "{include}"' + include_newline + edited
    return relative, header_relative, old, new, digest(raw)


def unified(old, new, before, after):
    return "".join(difflib.unified_diff(old.splitlines(keepends=True), new.splitlines(keepends=True),
                                        fromfile=before, tofile=after))


def repair(spec_path, root, output):
    root, output = root.resolve(), output.resolve()
    spec = load_spec(spec_path)
    relative, header_relative, old, new, source_hash = prepare_source(spec, root)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise GenerationError("Output directory must be absent or empty; evidence will not be overwritten.")
    output.mkdir(parents=True, exist_ok=True)
    header = output / (spec["name"] + ".hpp")
    args = argparse.Namespace(**{k: spec[k] for k in
        ("name", "expression", "variable", "center", "degree", "interval", "dtype")}, output=header)
    generate(args)
    header_bytes = header.read_bytes()
    header_text = header_bytes.decode("utf-8")
    patch = (unified("", header_text, "/dev/null", "b/" + str(header_relative))
             + unified(old, new, "a/" + str(relative), "b/" + str(relative)))
    # Write bytes to preserve target line endings, including on Windows.
    (output / "patch.diff").write_bytes(patch.encode("utf-8"))
    for path, content in ((relative, new.encode("utf-8")), (header_relative, header_bytes)):
        preview = output / "preview" / str(path)
        preview.parent.mkdir(parents=True, exist_ok=True)
        preview.write_bytes(content)
    write_json(output / "input.json", spec)
    report = {
        "schema_version": 1, "status": "generated_unvalidated_interval",
        "source": str(relative), "source_sha256": source_hash,
        "patched_source_sha256": digest(new.encode("utf-8")),
        "header": header.name, "header_destination": str(header_relative),
        "header_sha256": digest(header_bytes), "metadata": header.with_suffix(".json").name,
        "patch": "patch.diff", "patch_sha256": digest(patch.encode("utf-8")),
        "input_sha256": digest((output / "input.json").read_bytes()),
        "driver_sha256": digest(Path(__file__).read_bytes()),
        "source_modified": False,
        "manual_inputs": ["mathematical expression and its equivalence to the selected source function",
                          "expansion variable, center, degree and applicability interval",
                          "source file, function name and exact function definition"],
        "automatic_steps": ["symbolic expansion and coefficient extraction", "target-type rounding",
                            "same-type Horner code and interval guard emission",
                            "unique source matching and original-body fallback integration",
                            "unified diff, preview and provenance emission"],
        "validation": "Generation only. Compile/test preview or apply to a disposable copy before use. No whole-interval accuracy certificate.",
    }
    write_json(output / "manifest.json", report)
    return {"status": report["status"], "manifest": str(output / "manifest.json"),
            "patch": str(output / "patch.diff"), "source_modified": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(repair(args.spec, args.project_root, args.output), ensure_ascii=True))
    except (GenerationError, OSError, ValueError, TypeError, OverflowError) as exc:
        print(f"Repair generation failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
