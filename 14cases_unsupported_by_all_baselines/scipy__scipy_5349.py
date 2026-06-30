from __future__ import annotations

import math

from typing import Any

import mpmath as mp

import numpy as np

import scipy.special as sc



mp.mp.dps = 100

CASE_TITLE = 'hyp1f1 explodes for large arguments. Recurrence relation unstable.'

CASE_FUNCTION = 'scipy.special.hyp1f1'

CASE_ID = 'scipy/scipy#5349'

FINITE_DIFF_FIXED_STEP = mp.mpf("1e-5")

def relerr(v, ref):
    is_complex = isinstance(v, complex) or isinstance(ref, (complex, mp.ctx_mp_python.mpc))
    v = mp.mpc(v) if is_complex else mp.mpf(v)
    ref = mp.mpc(ref) if is_complex else mp.mpf(ref)
    if abs(ref) == 0:
        return float(abs(v - ref))
    return float(abs((v - ref) / ref))

def orders(before, after):
    if not math.isfinite(before) or before <= 0:
        return 'inf' if after == 0 or math.isfinite(after) else 'nan'
    if after == 0:
        return 'inf'
    if not math.isfinite(after):
        return '-inf'
    return f'{math.log10(before / after):.2f}'

def hyp1f1_taylor(a, b, z, tol=1e-16, max_terms=10000):
    z = complex(z)
    term = 1.0 + 0j
    total = term
    for n in range(1, max_terms + 1):
        term *= (float(a) + n - 1) / (float(b) + n - 1) * z / n
        total += term
        if abs(term) <= tol * max(1.0, abs(total)):
            return total
    return total

def _mp_value(x):
    if isinstance(x, (complex, mp.ctx_mp_python.mpc)):
        return mp.mpc(x)
    return mp.mpf(str(float(x))) if isinstance(x, np.generic) else mp.mpf(str(x))

def _finite_difference_step(x):
    return FINITE_DIFF_FIXED_STEP

def _condition_from_vector_finite_difference(f, xs):
    xs = [_mp_value(x) for x in xs]
    fx = f(xs)
    values = []
    steps = []
    for i, xi in enumerate(xs):
        h = _finite_difference_step(xi)
        for step in (h, -h):
            try:
                xh = list(xs)
                xh[i] = xh[i] + step
                d = (f(xh) - fx) / step
                if mp.isfinite(abs(d)):
                    values.append(abs(d) if abs(fx) == 0 else abs(xi * d / fx))
                    steps.append(step)
                    break
            except Exception:
                continue
    if not values:
        return (math.nan, math.nan)
    return (float(max(values)), max(steps, key=lambda s: abs(s)))

def _format_value(value: Any) -> str:
    try:
        if isinstance(value, (mp.ctx_mp_python.mpf, mp.ctx_mp_python.mpc)):
            return mp.nstr(value, 17)
        if isinstance(value, complex):
            return str(value)
        if isinstance(value, np.generic):
            return repr(value.item())
        return repr(value)
    except Exception:
        return str(value)

def case_scipy_scipy_5349() -> dict[str, Any]:
    """Compute and return reviewer evidence for scipy/scipy#5349."""
    case = 'scipy/scipy#5349'
    expressions = ('scipy.special.hyp1f1(a,b,z)', 'mpmath hyp1f1(a,b,z)', 'Kummer Taylor series sum')
    rows = []
    if sc is not None:
        a = 30
        b = 70
        z = 20j
        ref = mp.hyp1f1(a, b, z)
        original_value = sc.hyp1f1(a, b, z)
        taylor_value = hyp1f1_taylor(a, b, z)
        before = relerr(original_value, ref)
        after = relerr(taylor_value, ref)
        condition = _condition_from_vector_finite_difference(lambda xs: mp.hyp1f1(xs[0], xs[1], xs[2]), [30, 70, mp.mpc(0, 20)])
        rows.append({'condition_number': condition[0], 'original_expression': expressions[0], 'oracle_expression': expressions[1], 'taylor_expression': expressions[2], 'original_value': _format_value(original_value), 'oracle_value': _format_value(ref), 'taylor_value': _format_value(taylor_value), 'original_relative_error': before, 'taylor_relative_error': after, 'taylor_orders_improved': orders(before, after)})
    return (rows, expressions[2])

def print_case_report() -> None:
    rows, patch_idea = case_scipy_scipy_5349()
    print(f"{CASE_ID}: {CASE_TITLE}")
    print(f"  function: {CASE_FUNCTION}")
    print(f"  patch: {patch_idea}")
    for row in rows:
        print("  evidence:")
        print(f"    condition_number={float(row['condition_number']):.6e}")
        print(f"    original_expression={row['original_expression']}")
        print(f"    oracle_expression={row['oracle_expression']}")
        print(f"    taylor_expression={row['taylor_expression']}")
        print(f"    original_value={row['original_value']}")
        print(f"    oracle_value={row['oracle_value']}")
        print(f"    taylor_value={row['taylor_value']}")
        print(f"    original_relative_error={float(row['original_relative_error']):.6e}")
        print(f"    taylor_relative_error={float(row['taylor_relative_error']):.6e}")
        print(f"    taylor_orders_improved={row['taylor_orders_improved']}")

if __name__ == "__main__":
    print_case_report()
