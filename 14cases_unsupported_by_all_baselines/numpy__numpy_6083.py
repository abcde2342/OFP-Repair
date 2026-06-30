from __future__ import annotations

import math

from typing import Any

import mpmath as mp

import numpy as np

mp.mp.dps = 100

CASE_TITLE = 'BUG: arctan inaccurate for small complex numbers'

CASE_FUNCTION = 'numpy.arctan(complex)'

CASE_ID = 'numpy/numpy#6083'

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

def arctan_taylor(z):
    return z - z ** 3 / 3 + z ** 5 / 5

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

def print_case_report() -> None:
    z = 0.01 + 1e-14j
    expressions = (
        'np.arctan(z)',
        'mp.atan(mp.mpc(z.real, z.imag))',
        'arctan_taylor(z)',
    )

    oracle_value = mp.atan(mp.mpc(z.real, z.imag)).imag
    original_value = (0.0099996666866652394+9.9920072216263095e-15j).imag
    taylor_value = arctan_taylor(z).imag

    original_relative_error = relerr(original_value, oracle_value)
    taylor_relative_error = relerr(taylor_value, oracle_value)
    condition, _ = _condition_from_finite_difference(mp.atan, mp.mpc(z.real, z.imag))
    patch_idea = expressions[2]

    print(f"{CASE_ID}: {CASE_TITLE}")
    print(f"  function: {CASE_FUNCTION}")
    print(f"  patch: {patch_idea}")
    print("  evidence:")
    print(f"    condition_number={condition:.6e}")
    print(f"    original_expression={expressions[0]}")
    print(f"    oracle_expression={expressions[1]}")
    print(f"    taylor_expression={expressions[2]}")
    print(f"    original_value={_format_value(original_value)}")
    print(f"    oracle_value={_format_value(oracle_value)}")
    print(f"    taylor_value={_format_value(taylor_value)}")
    print(f"    original_relative_error={original_relative_error:.6e}")
    print(f"    taylor_relative_error={taylor_relative_error:.6e}")
    print(f"    taylor_orders_improved={orders(original_relative_error, taylor_relative_error)}")

if __name__ == "__main__":
    print_case_report()
