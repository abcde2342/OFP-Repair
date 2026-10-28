#!/usr/bin/env python3
"""Derive a local symbolic polynomial and emit same-type C++17 Horner code.

This portable tool reads no OFP-Repair source files or existing patches. Its
inputs are a mathematical expression and explicit expansion/interval settings.
It does NOT infer a repair strategy or certify the requested interval.
"""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import sys

import sympy as sp

FUNCTIONS = {name: getattr(sp, name) for name in (
    "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh",
    "asinh", "acosh", "atanh", "exp", "log", "sqrt", "erf", "erfc", "Abs")}
CONSTANTS = {"pi": sp.pi, "E": sp.E}


class GenerationError(ValueError):
    pass


def parse_math(text, variable=None):
    """A small mathematical grammar, not Python eval/sympify of user code."""
    if len(text) > 4096:
        raise GenerationError("Expression is too long (maximum 4096 characters).")
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as exc:
        raise GenerationError(f"Invalid mathematical expression: {exc.msg}") from exc

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            token = ast.get_source_segment(text, node).replace("_", "")
            try:
                result = sp.Rational(token)
            except (ValueError, TypeError) as exc:
                raise GenerationError("Use finite decimal numeric literals.") from exc
            if result.is_finite is not True:
                raise GenerationError("Numeric literals must be finite.")
            return result
        if isinstance(node, ast.Name):
            if variable is not None and node.id == variable.name:
                return variable
            if node.id in CONSTANTS:
                return CONSTANTS[node.id]
            raise GenerationError(f"Unknown symbol: {node.id}")
        if isinstance(node, ast.UnaryOp):
            value = visit(node.operand)
            if isinstance(node.op, ast.UAdd):
                return value
            if isinstance(node.op, ast.USub):
                return -value
        if isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return left+right
            if isinstance(node.op, ast.Sub):
                return left-right
            if isinstance(node.op, ast.Mult):
                return left*right
            if isinstance(node.op, ast.Div):
                return left/right
            if isinstance(node.op, ast.Pow):
                return left**right
            raise GenerationError("Supported operators: +, -, *, /, ** (use **, not ^).")
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id in FUNCTIONS and len(node.args) == 1 and not node.keywords):
            return FUNCTIONS[node.func.id](visit(node.args[0]))
        raise GenerationError("Only arithmetic and the documented one-argument math functions are accepted.")

    value = visit(tree.body)
    if value.has(sp.nan, sp.zoo, sp.oo, -sp.oo):
        raise GenerationError("Expression contains an undefined or infinite constant.")
    return value


def round_constant(expression, dtype):
    """Round a finite real constant to the selected binary format.

    Exact rationals use integer nearest-even rounding. Other constants are
    evaluated to 120 and 180 decimal digits and both must round identically.
    This consistency check is not a formal correctly-rounded transcendental proof.
    """
    if expression.free_symbols or expression.is_real is not True or expression.is_finite is not True:
        raise GenerationError(f"Expected a finite real constant, got {expression}.")
    p, emin, emax = (24, -126, 127) if dtype == "float" else (53, -1022, 1023)

    def nearest(q):
        if q == 0:
            return 0.0
        sign = -1 if q < 0 else 1
        q = abs(q)
        numerator, denominator = int(q.p), int(q.q)
        exponent = numerator.bit_length()-denominator.bit_length()
        def power_two(e):
            return sp.Rational(2)**e
        if q < power_two(exponent):
            exponent -= 1
        step_exponent = max(exponent, emin)-(p-1)
        scaled = q / power_two(step_exponent)
        integer, remainder = divmod(int(scaled.p), int(scaled.q))
        twice = 2*remainder
        if twice > scaled.q or (twice == scaled.q and integer % 2):
            integer += 1
        result = sign*sp.Integer(integer)*power_two(step_exponent)
        maximum = (2-power_two(1-p))*power_two(emax)
        if abs(result) > maximum:
            raise GenerationError(f"Constant overflows {dtype}: {expression}")
        value = float(result)  # Result is exactly representable in the target type.
        return value

    if expression.is_Rational:
        return nearest(expression)
    a = nearest(sp.Rational(str(sp.N(expression, 120))))
    b = nearest(sp.Rational(str(sp.N(expression, 180))))
    if a != b:
        raise GenerationError("Constant rounding is unresolved at 120/180 digits.")
    return a


def next_value(value, upward, dtype):
    if dtype == "double":
        return math.nextafter(value, math.inf if upward else -math.inf)
    bits = struct.unpack("!I", struct.pack("!f", value))[0]
    if value == 0:
        bits = 1 if upward else 0x80000001
    else:
        bits += 1 if (value > 0) == upward else -1
    return struct.unpack("!f", struct.pack("!I", bits))[0]


def literal(value, dtype):
    text = format(value, ".9g" if dtype == "float" else ".17g")
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text+("f" if dtype == "float" else "")


def bit_pattern(value, dtype):
    return struct.pack("!f" if dtype == "float" else "!d", value).hex()


def rounded_bound(expression, lower, dtype):
    value = round_constant(expression, dtype)
    difference = sp.Rational(value)-expression
    inward_test = difference.is_negative if lower else difference.is_positive
    if inward_test is None:
        raise GenerationError("Cannot establish inward rounding for this interval bound.")
    if inward_test:
        value = next_value(value, upward=lower, dtype=dtype)
    if not math.isfinite(value):
        raise GenerationError("Interval has no finite inward-rounded bound.")
    return value


def derive(expression, variable, center, degree):
    # Inspect the untruncated expression: a low degree can otherwise hide a
    # fractional-power branch or a nonanalytic term such as Abs(x)**3.
    # Meromorphic eligibility permits removable quotients; Poly below rejects
    # their nonremovable poles. Unknown eligibility is conservatively rejected.
    if expression.is_meromorphic(variable, center) is not True:
        raise GenerationError("Local meromorphic eligibility could not be established before truncation; a specialized analytic/branch/asymptotic strategy is required.")
    u = sp.Symbol("u", real=True)
    shifted = expression.subs(variable, u+center)
    try:
        right = sp.series(shifted, u, 0, degree+1, dir="+")
        left = sp.series(shifted, u, 0, degree+1, dir="-")
        polynomial = sp.Poly(sp.expand(right.removeO()), u)
        left_poly = sp.Poly(sp.expand(left.removeO()), u)
    except (ValueError, NotImplementedError, sp.PolynomialError, sp.PoleError) as exc:
        raise GenerationError(f"An ordinary two-sided polynomial could not be generated: {exc}") from exc
    if sp.simplify(polynomial.as_expr()-left_poly.as_expr()) != 0:
        raise GenerationError("Left and right local expansions differ; choose an explicit branch-specific strategy.")
    if polynomial.is_zero and sp.simplify(shifted) != 0:
        raise GenerationError("The requested truncation is identically zero; increase the degree or supply a different strategy.")
    coefficients = [sp.simplify(polynomial.nth(k)) for k in range(degree+1)]
    for coefficient in coefficients:
        if coefficient.free_symbols or coefficient.has(sp.Derivative, sp.Integral):
            raise GenerationError("An unresolved symbolic coefficient remains.")
    return right, coefficients


def generate(args):
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", args.name):
        raise GenerationError("Name must be an ASCII identifier beginning with a letter.")
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", args.variable) or args.variable in CONSTANTS or args.variable in FUNCTIONS:
        raise GenerationError("Choose a variable name distinct from mathematical constants/functions.")
    if not 0 <= args.degree <= 64:
        raise GenerationError("Degree must be between 0 and 64.")
    variable = sp.Symbol(args.variable, real=True)
    expression = parse_math(args.expression, variable)
    requested_center = parse_math(args.center)
    center = round_constant(requested_center, args.dtype)
    center_exact = sp.Rational(center)
    left_input, right_input = (parse_math(value) for value in args.interval)
    if (right_input-left_input).is_positive is not True:
        raise GenerationError("Interval lower bound must be provably less than upper bound.")
    lower, upper = rounded_bound(left_input, True, args.dtype), rounded_bound(right_input, False, args.dtype)
    if not lower <= center <= upper or not lower < upper:
        raise GenerationError("The rounded center must lie inside a nonempty representable interval.")
    series, exact_coefficients = derive(expression, variable, center_exact, args.degree)
    rounded = [round_constant(c, args.dtype) for c in exact_coefficients]
    last = args.degree
    while last > 0 and exact_coefficients[last] == 0:
        last -= 1
    center_literal = literal(center, args.dtype)
    lower_literal, upper_literal = literal(lower, args.dtype), literal(upper, args.dtype)
    name_hash = hashlib.sha256(args.name.encode("ascii")).hexdigest()[:16].upper()
    guard = "OFP_CODEGEN_"+args.name.upper()+"_"+name_hash+"_"+args.dtype.upper()+"_HPP"
    lines = ["// Generated by ofp_codegen: symbolic series -> same-type Horner arithmetic.",
             "// The interval is supplied by the caller and is NOT an accuracy certificate.",
             "// At a removable singularity the polynomial supplies its continuous extension.",
             "// Use ordinary IEEE arithmetic; do not compile with fast-math.",
             f"#ifndef {guard}", f"#define {guard}", "", "namespace ofp_codegen {", "",
             f"inline {args.dtype} {args.name}_series({args.dtype} x) {{",
             f"    const {args.dtype} u = x - {center_literal};",
             f"    {args.dtype} value = {literal(rounded[last], args.dtype)};"]
    for k in range(last-1, -1, -1):
        lines.append(f"    value = value * u + {literal(rounded[k], args.dtype)};")
    if last == 0:
        lines.append("    (void)u;")
    lines += ["    return value;", "}", "",
              f"inline {args.dtype} {args.name}_repaired({args.dtype} x,",
              f"        {args.dtype} (*original)({args.dtype})) {{",
              f"    if (x >= {lower_literal} && x <= {upper_literal}) {{",
              f"        return {args.name}_series(x);", "    }",
              "    return original(x);", "}", "", "} // namespace ofp_codegen", "", f"#endif // {guard}", ""]
    header = "\n".join(lines)
    out = args.output.resolve()
    if out.suffix != ".hpp":
        raise GenerationError("Output must have the .hpp extension.")
    meta = {
        "name": args.name, "expression": args.expression, "symbolic_expression": str(expression),
        "variable": args.variable, "requested_center": args.center, "center_literal": center_literal,
        "center_exact": str(center_exact), "center_bits": bit_pattern(center, args.dtype),
        "centering_policy": "Derive at the exact representable value of the emitted center, not at an unrepresentable mathematical constant.",
        "requested_interval": args.interval, "interval_literals": [lower_literal, upper_literal],
        "interval_policy": "Bounds rounded inward; inclusive tests; NaN/infinity/outside inputs call the provided original function.",
        "degree": args.degree, "dtype": args.dtype, "series": str(series),
        "series_variable": "u = input - exact emitted center", "generation_method": "SymPy pre-truncation meromorphic eligibility and two-sided local series; exact coefficient extraction; same-type Horner C++ emission",
        "coefficients": [{"degree": i, "symbolic": str(c), "literal": literal(v,args.dtype),
                          "bits": bit_pattern(v,args.dtype), "underflowed": bool(c != 0 and v == 0)}
                         for i,(c,v) in enumerate(zip(exact_coefficients,rounded))],
        "validation_status": "generated_unvalidated_interval",
        "limitations": ["The caller supplies expression, center, degree and applicability interval; the tool does not infer them.",
                        "Only supported one-variable, real, ordinary polynomial expansions are handled. No arbitrary C++ parsing, multivariable extraction, general asymptotic solver or automatic special-function mapping.",
                        "SymPy must establish local meromorphic eligibility before truncation; unknown eligibility is rejected, including some analytic functions such as erf/erfc. This is not a certificate for the whole supplied interval or the unsimplified expression domain.",
                        "Coefficient constants are evaluated offline; the emitted polynomial uses only the selected float/double type.",
                        "Non-rational coefficient rounding must agree at 120 and 180 decimal digits; this is not a certified transcendental rounding proof.",
                        "The callback preserves the original implementation outside the interval. Caller must supply a non-null function of the same input/output type.",
                        "The local approximation may supply a removable-singularity limit; no claim of preserving errno, floating exceptions or signed-zero behavior.",
                        "Sampling/compilation validation is separate; a supplied interval does not prove an error bound."],
        "header_sha256": hashlib.sha256(header.encode("utf-8")).hexdigest(),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sympy_version": sp.__version__, "python_version": sys.version,
        "generated_at": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(header, encoding="utf-8", newline="\n")
    out.with_suffix(".json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8", newline="\n")
    return {"header": str(out), "metadata": str(out.with_suffix('.json')), "status": meta["validation_status"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--variable", default="x")
    parser.add_argument("--center", required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--interval", nargs=2, required=True, metavar=("LOWER", "UPPER"))
    parser.add_argument("--dtype", choices=("float", "double"), default="double")
    parser.add_argument("--name", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(generate(args), ensure_ascii=True))
    except (GenerationError, TypeError, OverflowError) as exc:
        print(f"Generation failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
