from __future__ import annotations

import math

from typing import Any

import mpmath as mp

import numpy as np

mp.mp.dps = 100

CASE_TITLE = '[Bug] softplus returns inf for large beta while other frameworks return finite value'

CASE_FUNCTION = 'torch.nn.functional.softplus'

CASE_ID = 'pytorch/pytorch#171249'

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

def log1p_taylor_scalar(s, terms=8):
    total = 0.0
    power = 1.0
    for n in range(1, terms + 1):
        power *= s
        total += (-1.0) ** (n + 1) * power / n
    return total

def softplus_taylor(x, beta, terms=8):
    y = beta * x
    if y <= 0:
        t = math.exp(y)
        return log1p_taylor_scalar(t, terms=terms) / beta
    t = math.exp(-y)
    return x + log1p_taylor_scalar(t, terms=terms) / beta

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

def _mp_softplus(x, beta):
    y = _mp_value(beta) * _mp_value(x)
    if y >= 0:
        return _mp_value(x) + mp.log(1 + mp.e ** (-y)) / beta
    return mp.log(1 + mp.e ** y) / beta

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

def case_pytorch_pytorch_171249() -> dict[str, Any]:
    case = 'pytorch/pytorch#171249'
    expressions = ('log1p(exp(beta*x))/beta', 'mpmath softplus', 'x + (q - q**2/2 + q**3/3 - ...)/beta, q=exp(-beta*x)')
    rows = []
    beta = 1e+30
    # for val in [0.5, -1.0, 2.0]:
    for val in [2.0]:
        ref = mp.log(1 + mp.e ** (mp.mpf(str(beta)) * mp.mpf(str(val)))) / mp.mpf(str(beta))
        bval = math.inf
        before = relerr(bval, ref) if math.isfinite(bval) else math.inf
        taylor_value = softplus_taylor(val, beta)
        after = relerr(taylor_value, ref)
        condition = _condition_from_finite_difference(lambda xx: _mp_softplus(xx, beta), val)
        rows.append({'condition_number': condition[0], 'original_expression': expressions[0], 'oracle_expression': expressions[1], 'taylor_expression': expressions[2], 'original_value': _format_value(bval), 'oracle_value': _format_value(ref), 'taylor_value': _format_value(taylor_value), 'original_relative_error': before, 'taylor_relative_error': after, 'taylor_orders_improved': orders(before, after)})
    return (rows, expressions[2])

def print_case_report() -> None:
    rows, patch_idea = case_pytorch_pytorch_171249()
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
