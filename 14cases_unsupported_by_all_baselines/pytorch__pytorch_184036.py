from __future__ import annotations

import math

from typing import Any

import mpmath as mp

import numpy as np
import torch

mp.mp.dps = 100

CASE_TITLE = 'torch.nn.PairwiseDistance returns +Inf for moderately large p (overflow in fp32, no log-space rewrite)'

CASE_FUNCTION = 'torch.nn.PairwiseDistance / vector_norm'

CASE_ID = 'pytorch/pytorch#184036'

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

def pnorm_naive(vals, p):
    arr = torch.as_tensor(vals, dtype=torch.float64)
    x1 = arr.reshape(1, -1)
    x2 = torch.zeros_like(x1)
    out = torch.nn.PairwiseDistance(p=p)(x1, x2)
    return float(out)

def pnorm_taylor(vals, p, terms=6):
    arr = np.abs(np.asarray(vals, dtype=np.float64))
    
    scale = float(np.max(arr))
    if scale == 0:
        return 0.0
    ratios = arr / scale
    max_index = int(np.argmax(ratios))
    small = np.delete(ratios, max_index) ** float(p)
    s = float(np.sum(small))
    alpha = 1.0 / float(p)
    total = 1.0
    coeff = 1.0
    power = 1.0
    for n in range(1, terms + 1):
        coeff *= (alpha - (n - 1)) / n
        power *= s
        total += coeff * power
    return scale * total

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

def _mp_pnorm(xs, p):
    return mp.fsum([abs(_mp_value(v)) ** p for v in xs]) ** (mp.mpf(1) / p)

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

def case_pytorch_pytorch_184036() -> dict[str, Any]:
    case = 'pytorch/pytorch#184036'
    expressions = ('(sum(abs(x_i)**p))**(1/p)', 'mpmath p-norm', 'max(abs(x))*(1+s)**(1/p) binomial Taylor')
    vals = [2.5577, 2.5867, 2.5688, 3.2271]
    p = 100.0
    ref = mp.fsum([mp.mpf(str(abs(v))) ** p for v in vals]) ** (mp.mpf(1) / p)
    original_value = pnorm_naive(vals, p)
    taylor_value = pnorm_taylor(vals, p)
    before = relerr(original_value, ref)
    after = relerr(taylor_value, ref)
    condition = _condition_from_vector_finite_difference(lambda xs: _mp_pnorm(xs, mp.mpf('100')), vals)
    row = {'condition_number': condition[0], 'original_expression': expressions[0], 'oracle_expression': expressions[1], 'taylor_expression': expressions[2], 'original_value': _format_value(original_value), 'oracle_value': _format_value(ref), 'taylor_value': _format_value(taylor_value), 'original_relative_error': before, 'taylor_relative_error': after, 'taylor_orders_improved': orders(before, after)}
    return ([row], expressions[2])

def print_case_report() -> None:
    rows, patch_idea = case_pytorch_pytorch_184036()
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
