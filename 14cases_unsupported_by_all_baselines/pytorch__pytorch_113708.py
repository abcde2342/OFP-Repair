from __future__ import annotations

import math

from typing import Any
import torch
import mpmath as mp

import numpy as np

mp.mp.dps = 100

CASE_TITLE = '`log_softmax` could be `2**124` to `2**1021` times more accurate on small outputs'

CASE_FUNCTION = 'torch.log_softmax'

CASE_ID = 'pytorch/pytorch#113708'

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

def log_softmax_naive(x):
    x_np = np.asarray(x)
    x_torch = torch.tensor(x_np, dtype=torch.float32)
    return torch.log_softmax(x_torch, dim=0).detach().cpu().numpy()

def log1p_taylor_scalar(s, terms=8):
    total = 0.0
    power = 1.0
    for n in range(1, terms + 1):
        power *= s
        total += (-1.0) ** (n + 1) * power / n
    return total

def log_softmax_taylor(x, terms=8):
    x = np.asarray(x)
    i = int(np.argmax(x))
    off = x - x[i]
    ex = np.exp(off)
    ex[i] = 0
    s = float(np.sum(ex))
    return off - log1p_taylor_scalar(s, terms=terms)

def log_softmax_mp(x):
    vals = [mp.mpf(str(float(v))) for v in x]
    m = max(vals)
    s = mp.fsum([mp.e ** (v - m) for v in vals])
    return [v - m - mp.log(s) for v in vals]

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

def _mp_log_softmax_component(xs, index=0):
    vals = [_mp_value(v) for v in xs]
    m = max(vals)
    return vals[index] - m - mp.log(mp.fsum([mp.e ** (v - m) for v in vals]))

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

def _log_softmax_dominant_case(case: str, expressions: tuple[str, str, str]) -> dict[str, Any]:
    x = np.array([1 - np.log(2 * np.finfo(np.float32).eps), 0], dtype=np.float32)
    ref = log_softmax_mp(x)
    original_value = float(log_softmax_naive(x)[0])
    taylor_value = float(log_softmax_taylor(x)[0])
    before = relerr(original_value, ref[0])
    after = relerr(taylor_value, ref[0])
    condition = _condition_from_vector_finite_difference(lambda xs: _mp_log_softmax_component(xs, 0), [1 - np.log(2 * np.finfo(np.float32).eps), 0])
    row = {'condition_number': condition[0], 'original_expression': expressions[0], 'oracle_expression': expressions[1], 'taylor_expression': expressions[2], 'original_value': _format_value(original_value), 'oracle_value': _format_value(ref[0]), 'taylor_value': _format_value(taylor_value), 'original_relative_error': before, 'taylor_relative_error': after, 'taylor_orders_improved': orders(before, after)}
    return ([row], expressions[2])

def case_pytorch_pytorch_113708() -> dict[str, Any]:
    return _log_softmax_dominant_case('pytorch/pytorch#113708', ('x_i - max(x) - log(sum(exp(x-max(x))))', 'mpmath shifted logsumexp', 'x_i-max(x)-(s - s**2/2 + s**3/3 - ...), s=sum(exp(nonmax))'))

def print_case_report() -> None:
    rows, patch_idea = case_pytorch_pytorch_113708()
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
