from __future__ import annotations

import math

from typing import Any

import mpmath as mp

import numpy as np
import scipy

mp.mp.dps = 100

CASE_TITLE = 'BUG: stats.tukeylambda: Bad behavior of the `cdf()` and `sf()` methods in the tails.'

CASE_FUNCTION = 'scipy.stats.tukeylambda.ppf / quantile Q(p, lambda)'

CASE_ID = 'scipy/scipy#21370'

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

def tukey_q_direct(p, lam):
    return scipy.stats.tukeylambda.ppf(p, lam)

def tukey_q_taylor(p, lam, terms=6):
    u = 2.0 * (p - 0.5)
    coeff = 2.0 ** (1.0 - lam)
    s = 0.0
    prod = 1.0
    fact = 1.0
    for k in range(terms):
        if k == 0:
            term_coeff = 1.0
            power = u
        else:
            for j in range(2 * k - 1, 2 * k + 1):
                prod *= lam - j
            fact *= 2 * k * (2 * k + 1)
            term_coeff = prod / fact
            power = u ** (2 * k + 1)
        s += term_coeff * power
    return coeff * s

def tukey_q_mp(p, lam):
    p = mp.mpf(float(p))
    lam = mp.mpf(str(lam))
    return (p ** lam - (1 - p) ** lam) / lam

def _mp_value(x):
    if isinstance(x, (complex, mp.ctx_mp_python.mpc)):
        return mp.mpc(x)
    return mp.mpf(str(float(x))) if isinstance(x, np.generic) else mp.mpf(str(x))

def _finite_difference_step(x):
    return FINITE_DIFF_FIXED_STEP

def _finite_difference_derivative(f, x):
    x = _mp_value(x)
    fx = f(x)
    h = _finite_difference_step(x)

    def candidate(step):
        return ((f(x + step) - fx) / step, step)
    for step in (h, -h):
        try:
            d, used = candidate(step)
            if mp.isfinite(abs(d)):
                return (d, used, fx)
        except Exception:
            continue
    return (mp.nan, h, fx)

def _condition_from_finite_difference(f, x):
    d, h, fx = _finite_difference_derivative(f, x)
    x = _mp_value(x)
    if abs(fx) == 0:
        value = abs(d)
    else:
        value = abs(x * d / fx)
    return (float(value), h)

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

def case_scipy_scipy_21370() -> dict[str, Any]:
    case = 'scipy/scipy#21370'
    expressions = ('(p**lambda - (1-p)**lambda)/lambda', 'mpmath same quantile formula', 'Taylor of Q(p, lambda) around p=0.5')
    rows = []
    for p, lam in [(0.50001, 0.2)]:
        ref = tukey_q_mp(p, lam)
        original_value = tukey_q_direct(float(p), lam)
        taylor_value = tukey_q_taylor(float(p), lam)
        before = relerr(original_value, ref)
        after = relerr(taylor_value, ref)
        condition = _condition_from_finite_difference(lambda pp, ll=lam: tukey_q_mp(pp, ll), p)
        rows.append({'condition_number': condition[0], 'original_expression': expressions[0], 'oracle_expression': expressions[1], 'taylor_expression': expressions[2], 'original_value': _format_value(original_value), 'oracle_value': _format_value(ref), 'taylor_value': _format_value(taylor_value), 'original_relative_error': before, 'taylor_relative_error': after, 'taylor_orders_improved': orders(before, after)})
    return (rows, expressions[2])

def print_case_report() -> None:
    rows, patch_idea = case_scipy_scipy_21370()
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
