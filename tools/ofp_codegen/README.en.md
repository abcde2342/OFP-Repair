# Automated symbolic patch generation

`repair.py` takes a JSON specification and produces a C++ repair, a source-code
diff, and a machine-readable derivation record. It automatically derives the
polynomial coefficients using SymPy; it does not copy the coefficients from the
existing OFP-Repair implementations. The generated numerical code uses the
specified `float` or `double` type throughout.

This is **automated synthesis given an explicit expansion and localization
specification**. The user supplies the expression, expansion variable/center,
degree, interval, and source function. These choices are not inferred. This
tool is an additional reproducible artifact; it does not establish that the
paper's existing 57 repairs were originally synthesized automatically, or that
all 57 cases are supported by this implementation.

## Install and run

Requirements: Python >= 3.10. Generation needs SymPy and mpmath; validation also
needs Git and a C++17 compiler (GCC or Clang). MSVC `cl` is not supported by the
validation script. Run these commands from the repository root:

```bash
python -m pip install -r tools/ofp_codegen/requirements.txt
python tools/ofp_codegen/repair.py --spec tools/ofp_codegen/examples/cos_x2.input.json --project-root . --output runs/cos_x2-patch
```

The output directory must be absent or empty. The target source is left intact.
The supplied example can be compiled directly from the generated preview:

```bash
c++ -std=c++17 -O2 -fno-fast-math -ffp-contract=off runs/cos_x2-patch/preview/tools/ofp_codegen/examples/original.cpp -o runs/cos_x2-patch/demo
./runs/cos_x2-patch/demo
```

To check and apply the diff to a disposable checkout of the original project:

```bash
git apply --check runs/cos_x2-patch/patch.diff
git apply runs/cos_x2-patch/patch.diff
```

The diff adds `ofp_generated/cos_x2.hpp` and rewrites the selected example
function. It retains the original function body in a noncapturing lambda that
is called outside the repair interval, including for NaN and infinities.
`examples/project_cos_x2.input.json` instead targets the existing
`cos_x2` function in `inputs/targets.cpp`. Use that specification with a new
output directory to generate a patch for the actual benchmark source. Applying
that patch changes the benchmark's *original* function, so use a separate
checkout for it rather than changing the baseline of an evaluation run.

## Explicit input contract

The complete [example specification](examples/cos_x2.input.json) is:

```json
{
  "schema_version": 1,
  "name": "cos_x2",
  "expression": "(1-cos(x))/x**2",
  "variable": "x",
  "center": "0",
  "degree": 10,
  "interval": ["-0.01", "0.01"],
  "dtype": "double",
  "target": {
    "path": "tools/ofp_codegen/examples/original.cpp",
    "function": "cos_x2",
    "definition": "double cos_x2(double x) {\n    return (1.0 - std::cos(x)) / (x * x);\n}"
  }
}
```

| Field | Meaning |
| --- | --- |
| `schema_version` | Integer `1`; unknown or missing fields are rejected. |
| `name` | Identifier used for generated functions and header filename. |
| `expression` | Real mathematical expression; decimal literals are exact rationals. |
| `variable` | Single scalar expansion variable and C++ parameter name. |
| `center` | Mathematical string; rounded to the selected type before expansion. |
| `degree` | Maximum polynomial degree, integer from 0 through 64. |
| `interval` | Two mathematical strings specifying the closed repair interval. |
| `dtype` | `float` (binary32) or `double` (binary64). |
| `target.path` | Project-relative C++ source path, resolved under `--project-root`. |
| `target.function` | Name of the selected scalar free function. |
| `target.definition` | Exact function definition, including signature and body. |

The selected definition must match exactly once on whole source lines. LF and
CRLF (including mixed files) are supported; untouched source bytes are preserved.
The replacement uses the selected function's line-ending style. Source must be UTF-8 without BOM and end with
a newline. Paths use `/` and ASCII letters/digits/`_.-`, with no traversal or
`.git` components. An existing generated destination header is rejected.

Source integration supports a simple one-argument free function with the same
input/output type (`double f(double x)` or `float f(float x)`). It assumes a
pure scalar computation: no recursion, parameter mutation with external effects,
macros changing the function structure, function-local identity dependence
(`__func__`), or requirements to preserve floating-point exception flags. The
caller establishes that the mathematical expression models this source body.
There is no C++ parser or proof of source/expression equivalence.

## Explicit output contract

| Output | Contents |
| --- | --- |
| `<name>.hpp` | Generated same-type polynomial, interval guard and original-function callback. |
| `<name>.json` | Exact symbolic coefficients, rounded literals and IEEE bits, expansion point, interval, versions and code hashes. |
| `patch.diff` | Unified diff adding the header and replacing the selected source definition. |
| `preview/` | New header and modified source at their project-relative paths; other project dependencies are not copied. |
| `input.json` | Snapshot of the complete validated specification. |
| `manifest.json` | Input/source/patch/header hashes, output paths, manual inputs, automatic steps and validation status. |

The command prints a JSON summary and exits `0` on successful generation; invalid
input or source mismatch exits `2`. Success is labeled
`generated_unvalidated_interval`: it is not a claim of a whole-interval accuracy
bound. Validation reports are produced separately by `validate.py`.

## Generation algorithm

1. Validate the JSON and locate the exact, unique source definition.
2. Parse the mathematical expression using an allowlisted AST grammar (no eval).
3. Round the expansion center to the target type; expand about the exact value
   of that represented number. Round interval bounds inward.
4. Require SymPy to establish local meromorphic eligibility before truncation.
   Compute the left and right series to degree + 1, require matching ordinary
   polynomials, and reject poles, fractional powers and unresolved coefficients.
5. Extract exact coefficients and round them to binary32/binary64. Rationals
   use integer nearest-even rounding; transcendental constants must agree when
   evaluated at 120 and 180 decimal digits.
6. Emit Horner evaluation, the interval guard, and the source wrapper. Preserve
   the original body for fallback. Save the diff, preview and derivation evidence.

For the supplied input, the polynomial is generated as
`1/2 - x**2/24 + x**4/720 - x**6/40320 + x**8/3628800 - x**10/479001600`.
At `x = 1e-8`, the original double expression can return zero from cancellation,
while the generated result is approximately `0.5`. At `x = 0`, the polynomial
returns the removable-singularity limit `0.5` rather than the original NaN.
Outside the inward-rounded interval the original computation is used.

Offline symbolic derivation, constant evaluation and the test oracle can use
high precision. The emitted numerical implementation does not use a higher
precision type. Compile without fast-math; the tests disable FMA contraction.

## Reproducible validation

```bash
python tools/ofp_codegen/validate.py --cxx c++ --output runs/codegen-validation
```

This executes the actual generator, compiler and generated binaries. It checks:

- Six generated implementations: three expressions (`cos_x2`, `x_x2`, and
  `sin` at a nonzero center), each in `float` and `double`.
- 54 in-interval numerical samples against independent mpmath references at
  100 decimal digits; 30 boundary/NaN/infinity fallback checks.
- Strict compilation, return types, exact represented centers and include guards.
- Rejection of ten invalid/nonanalytic cases, including branch terms hidden by
  a low requested degree.
- End-to-end source patches with LF, CRLF and mixed line endings: JSON input, generation, Git diff
  check/application, compilation of original and patched sources, actual error
  improvement, fallback, artifact hashes, stale-source and output-overwrite rejection.

`run-report.json` records results, oracle values and numerical errors. Each
external command retains its arguments, stdout/stderr, exit status and timing.
The test output includes generated sources and binaries. `.github/workflows/codegen.yml`
runs this same validation on Linux for changes to the generator.

On Windows/MinGW, put the compiler's `bin` directory on `PATH` so its runtime
DLLs can be found, and pass the path to `g++.exe` as `--cxx` if necessary.

## Scope and limitations

The expression grammar supports `+ - * / **`, `pi`, `E`, and one-argument
`sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh exp log sqrt erf
erfc Abs`. Acceptance still depends on the local symbolic eligibility check;
some analytic functions such as `erf` can be conservatively rejected.

The tool does not infer the unstable subexpression, expansion strategy, degree
or valid repair radius. It does not provide general multivariate synthesis,
arbitrary special-function mappings, variable substitutions, asymptotic
expansions, or Newton-based repairs. These require explicit additional strategies.
The legacy `script/write_patches.py` exports ACESO fitted coefficients and is a
different path from this symbolic generator.

Local eligibility, a generated polynomial, and finite validation samples do not
prove accuracy across the full interval, preservation of the original expression's
domain, errno, signed zero, or floating-point exception behavior. In particular,
removable singularities are extended continuously. Changing the interval or
compiler settings requires appropriate validation. The two-precision constant
check is not a formal proof of correctly rounded transcendental coefficients.

For header-only generation without source integration, use `generate.py` as
described in the [Chinese documentation](README.md).
