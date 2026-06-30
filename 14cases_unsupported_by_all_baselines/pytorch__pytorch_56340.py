from __future__ import annotations

import math

from typing import Any

import mpmath as mp

import numpy as np
import torch

mp.mp.dps = 100

CASE_TITLE = 'log_softmax(x) != x - logsumexp(x)'

CASE_FUNCTION = 'torch.nn.functional.log_softmax / logsumexp-based manual formula'

CASE_ID = 'pytorch/pytorch#56340'

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

def log_softmax_mp(x):
    vals = [mp.mpf(str(float(v))) for v in x]
    m = max(vals)
    s = mp.fsum([mp.e ** (v - m) for v in vals])
    return [v - m - mp.log(s) for v in vals]

def log_softmax(x, dim=-1):
    return x - torch.logsumexp(x, dim=dim, keepdim=True)

def log_softmax_original_component(x, index=0):
    a = torch.tensor(np.asarray(x, dtype=np.float32))
    pt = log_softmax(a, dim=-1)
    return float(pt[index].item())

def log1p_taylor_tensor(s, terms=8):
    total = torch.zeros((), dtype=s.dtype, device=s.device)
    power = torch.ones((), dtype=s.dtype, device=s.device)
    for n in range(1, terms + 1):
        power *= s
        sign = 1.0 if n % 2 == 1 else -1.0
        total = total + sign * power / n
    return total

def log_softmax_taylor_equalmax_component(x, index=0, terms=8):
    a = torch.tensor(np.asarray(x, dtype=np.float32))
    m = torch.max(a)
    shifted = a - m
    max_mask = shifted == 0
    k = torch.sum(max_mask).to(dtype=a.dtype)
    tail = torch.sum(torch.exp(shifted[~max_mask]))
    z = tail / k
    return float((shifted[index] - torch.log(k) - log1p_taylor_tensor(z, terms=terms)).item())

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

def case_pytorch_pytorch_56340() -> dict[str, Any]:
    case = 'pytorch/pytorch#56340'
    expressions = ('log_softmax(a): a - torch.logsumexp(a, dim=-1, keepdim=True)', 'mpmath shifted logsumexp', 'torch tensor shifted equal-max log(k)+(z - z**2/2 + z**3/3 - ...), z=tail/k')
    xvec = np.array([800, 800, 400], dtype=np.float32)
    ref0 = log_softmax_mp(xvec)[0]
    original_value = log_softmax_original_component(xvec, index=0)
    taylor_value = log_softmax_taylor_equalmax_component(xvec, index=0)
    before = relerr(original_value, ref0)
    after = relerr(taylor_value, ref0)
    condition = _condition_from_vector_finite_difference(lambda xs: _mp_log_softmax_component(xs, 0), [800, 800, 400])
    row = {'condition_number': condition[0], 'original_expression': expressions[0], 'oracle_expression': expressions[1], 'taylor_expression': expressions[2], 'original_value': _format_value(original_value), 'oracle_value': _format_value(ref0), 'taylor_value': _format_value(taylor_value), 'original_relative_error': before, 'taylor_relative_error': after, 'taylor_orders_improved': orders(before, after)}
    return ([row], expressions[2])

def print_case_report() -> None:
    rows, patch_idea = case_pytorch_pytorch_56340()
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
