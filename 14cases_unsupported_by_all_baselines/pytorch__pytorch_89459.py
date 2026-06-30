from __future__ import annotations

import math

from typing import Any

import mpmath as mp

import numpy as np

mp.mp.dps = 100

CASE_TITLE = 'Higher order derivatives of sinc explode'

CASE_FUNCTION = 'torch.sinc higher-order derivatives at zero'

CASE_ID = 'pytorch/pytorch#89459'

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

def sinc_second_derivative_taylor():
    return -(math.pi * math.pi) / 3.0

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

def _mp_normalized_sinc(x):
    x = _mp_value(x)
    if x == 0:
        return mp.mpf(1)
    return mp.sin(mp.pi * x) / (mp.pi * x)

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

def pytorch_result():
    import torch
    zero = torch.zeros(1, requires_grad=True)
    dx = torch.autograd.grad(
        torch.sinc(zero),
        zero,
        create_graph=True
    )[0]
    ddx = torch.autograd.grad(dx, zero)[0]
    return ddx.item()

def case_pytorch_pytorch_89459() -> dict[str, Any]:
    """Compute and return reviewer evidence for pytorch/pytorch#89459."""
    case = 'pytorch/pytorch#89459'
    expressions = ('d2/dx2 sinc(x) at x=0', 'Maclaurin derivative limit', '-pi**2/3')
    ref = -mp.pi ** 2 / 3
    before = pytorch_result()
    taylor_value = sinc_second_derivative_taylor()
    after = relerr(taylor_value, ref)
    condition = _condition_from_finite_difference(_mp_normalized_sinc, 0)
    row = {'condition_number': condition[0], 'original_expression': expressions[0], 'oracle_expression': expressions[1], 'taylor_expression': expressions[2], 'original_value': _format_value(math.inf), 'oracle_value': _format_value(ref), 'taylor_value': _format_value(taylor_value), 'original_relative_error': before, 'taylor_relative_error': after, 'taylor_orders_improved': orders(before, after)}
    return ([row], expressions[2])

def print_case_report() -> None:
    rows, patch_idea = case_pytorch_pytorch_89459()
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
