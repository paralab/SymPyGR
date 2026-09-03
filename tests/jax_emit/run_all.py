#!/usr/bin/env python
"""Gate B -- JAX printer self-checks.

    python tests/jax_emit/run_all.py

Printer machinery only, not any solver: printed source must be numerically
faithful to its sympy expression, valid Python, and must fire the four
transliteration rules recovered from the bssneqs_SSL_HD_dxsq.cpp diff. Plain
asserts, same shape as tests/numgate/run_all.py. numpy stands in for jnp.
"""

import ast
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import numpy as np  # noqa: E402
import sympy as sym  # noqa: E402

import dendrosym.derivs  # noqa: E402
from dendrosym.jax_printer import (  # noqa: E402
    DendroJaxPrinter,
    jax_symbol_name,
    safe_name,
)

FAILURES = []


def check(name, fn):
    try:
        fn()
    except AssertionError as exc:
        FAILURES.append(name)
        print(f"FAIL  {name}\n      {exc}")
    except Exception as exc:  # a raising gate is a failing gate
        FAILURES.append(name)
        print(f"FAIL  {name}\n      {type(exc).__name__}: {exc}")
    else:
        print(f"PASS  {name}")


# ---------------------------------------------------------------------------
# the expression battery: one entry per emission rule we care about
# ---------------------------------------------------------------------------
def _battery():
    a, b, c = sym.symbols("alpha[pp] chi[pp] gt0[pp]")
    return [
        ("plain mul/add", 2 * a + 3 * b * c),
        ("rational coeff", sym.Rational(3, 4) * a + sym.Rational(2, 3) * b),
        ("int pow 2", a ** 2),
        ("int pow 6", a ** 6),
        ("int pow -1", b ** -1),
        ("int pow -2", b ** -2),
        ("sqrt", sym.sqrt(b)),
        ("exp", sym.exp(-a * a)),
        ("log", sym.log(b)),
        ("nested", sym.sqrt(b) * (a - sym.sqrt(b)) * sym.exp(-a / 2)),
        ("max", sym.Max(b, sym.Float(1e-4))),
        ("min", sym.Min(a, sym.Float(2.0))),
        ("max 3-ary", sym.Max(a, b, c)),
        ("mixed deep", (a * b + c ** 3) / (sym.sqrt(b) + a ** -2)),
    ]


def test_numeric_fidelity():
    """Printed source, evaluated, must equal sympy's evaluation."""
    p = DendroJaxPrinter()
    rng = np.random.default_rng(20260903)
    syms = sorted({s for _, e in _battery() for s in e.free_symbols}, key=str)
    worst = 0.0
    for trial in range(8):
        # positive: chi/alpha are positive-definite; sqrt/log of a negative
        # would compare nan to nan.
        vals = {s: float(rng.uniform(0.25, 2.0)) for s in syms}
        env = {"jnp": np, **{jax_symbol_name(s.name): v for s, v in vals.items()}}
        for name, expr in _battery():
            src = p.doprint(expr)
            got = eval(src, {"__builtins__": {}}, env)  # noqa: S307
            want = float(expr.subs(vals).evalf(30))
            denom = max(abs(want), 1.0)
            rel = abs(got - want) / denom
            worst = max(worst, rel)
            assert rel < 1e-14, (
                f"{name!r} trial {trial}: printed {src!r} -> {got!r}, "
                f"sympy -> {want!r} (rel {rel:.3e})"
            )
    print(f"      worst relative error over {8 * len(_battery())} evals: {worst:.3e}")


def test_output_is_valid_python():
    p = DendroJaxPrinter()
    for name, expr in _battery():
        src = p.doprint(expr)
        try:
            ast.parse(src, mode="eval")
        except SyntaxError as exc:
            raise AssertionError(f"{name!r} printed unparseable source {src!r}: {exc}")


def test_no_numpy_leakage():
    """NumPyPrinter bakes 'numpy.' into 29 targets; none may survive."""
    p = DendroJaxPrinter()
    leaked = {k: v for k, v in p.known_functions.items()
              if isinstance(v, str) and "numpy." in v}
    assert not leaked, f"known_functions still targets numpy: {leaked}"
    for name, expr in _battery():
        src = p.doprint(expr)
        assert "numpy." not in src, f"{name!r} emitted a numpy call: {src!r}"


def test_point_index_stripped():
    p = DendroJaxPrinter()
    src = p.doprint(sym.Symbol("alpha[pp]") * sym.Symbol("grad_0_gt0[pp]"))
    assert "[pp]" not in src, f"point index survived: {src!r}"
    assert "alpha" in src and "grad_0_gt0" in src, src


def test_keyword_mangling():
    """R3: `lambda` is declared by bssn, ccz4 and emda."""
    assert safe_name("lambda") == "lambda_param"
    assert safe_name("alpha") == "alpha"
    # the point index goes, a real subscript stays, identifier gets mangled
    assert jax_symbol_name("lambda[0]") == "lambda_param[0]"
    assert jax_symbol_name("lambda_f[1]") == "lambda_f[1]"
    assert jax_symbol_name("alpha[pp]") == "alpha"
    p = DendroJaxPrinter()
    src = p.doprint(sym.Symbol("lambda[0]") * sym.Symbol("alpha[pp]"))
    assert "lambda_param[0]" in src, src
    ast.parse(src, mode="eval")  # would raise on a bare `lambda[0]`


def test_keyword_scan_covers_every_config():
    """General guard, not a `lambda` special case."""
    import keyword
    for kw in ("lambda", "class", "is", "in", "and", "not", "None", "True"):
        assert keyword.iskeyword(kw)
        assert safe_name(kw) == kw + "_param"
        ast.parse(safe_name(kw) + " = 1")


def test_integer_pow_lowering():
    """R1: no ** or pow() in the ranges the C printer lowers."""
    p = DendroJaxPrinter()
    b = sym.Symbol("chi[pp]")
    for e in list(range(2, 7)) + list(range(-6, 0)):
        src = p.doprint(b ** e)
        assert "**" not in src, f"exp {e} kept a power operator: {src!r}"
        assert "pow" not in src, f"exp {e} emitted a pow call: {src!r}"
    # non-integer exponents fall through to the base printer
    src = p.doprint(b ** sym.Symbol("p_expo"))
    assert "**" in src or "power" in src, src


def test_max_min_are_binary():
    """The inherited amax(asarray([...])) form allocates."""
    p = DendroJaxPrinter()
    src = p.doprint(sym.Max(sym.Symbol("chi[pp]"), sym.Float(1e-4)))
    assert "jnp.maximum" in src, src
    assert "asarray" not in src and "amax" not in src, src
    src = p.doprint(sym.Min(sym.Symbol("chi[pp]"), sym.Float(1.0)))
    assert "jnp.minimum" in src, src
    assert "asarray" not in src and "amin" not in src, src


def test_rational_prints_float_pair():
    p = DendroJaxPrinter()
    src = p.doprint(sym.Rational(3, 4) * sym.Symbol("alpha[pp]"))
    assert "3.0/4.0" in src, src


def test_grad2_mixed_indices_do_not_collapse():
    """Guard for the grad2 index swap fixed in 015339c."""
    p = DendroJaxPrinter()
    f = sym.Function("alpha")(dendrosym.derivs.xx,
                              dendrosym.derivs.yy,
                              dendrosym.derivs.zz)
    # d/dy d/dx -- decreasing order on input, must come out as (0, 1)
    src = p.doprint(sym.Derivative(f, dendrosym.derivs.yy, dendrosym.derivs.xx))
    assert "grad2(0, 1," in src, f"mixed deriv collapsed or mis-ordered: {src!r}"
    src = p.doprint(sym.Derivative(f, dendrosym.derivs.xx, dendrosym.derivs.zz))
    assert "grad2(0, 2," in src, src
    # pure second derivative keeps the repeated index
    src = p.doprint(sym.Derivative(f, dendrosym.derivs.yy, dendrosym.derivs.yy))
    assert "grad2(1, 1," in src, src


def test_first_derivative_form():
    p = DendroJaxPrinter()
    f = sym.Function("chi")(dendrosym.derivs.xx,
                            dendrosym.derivs.yy,
                            dendrosym.derivs.zz)
    src = p.doprint(sym.Derivative(f, dendrosym.derivs.zz))
    assert "grad(2," in src, src


def test_matches_c_printer_where_it_should():
    """Printers agree once the known rule deltas are undone.

    Same normalization the oracle diff used.
    """
    import re
    from dendrosym.code_printer import DendroCPrinter

    cp, jp = DendroCPrinter(), DendroJaxPrinter()

    def norm(s):
        s = s.replace(" ", "").replace("[pp]", "")
        s = s.replace("jnp.", "").replace("numpy.", "")
        s = s.replace("lambda_param", "lambda")
        s = re.sub(r"\((\d+\.\d+/\d+\.\d+)\)", r"\1", s)   # R2
        s = re.sub(r"\(([A-Za-z_0-9]+)\*\1\)", r"\1*\1", s)  # R1
        return s

    mismatched = []
    for name, expr in _battery():
        if expr.has(sym.Max) or expr.has(sym.Min):
            continue  # the C printer has no jnp.maximum analogue; checked above
        a, b = norm(cp.doprint(expr)), norm(jp.doprint(expr))
        if a != b:
            mismatched.append((name, a, b))
    assert not mismatched, "printers diverge beyond the known rules:\n" + "\n".join(
        f"  {n}\n    c  : {a}\n    jax: {b}" for n, a, b in mismatched
    )


if __name__ == "__main__":
    check("numeric fidelity vs sympy", test_numeric_fidelity)
    check("output is valid python", test_output_is_valid_python)
    check("no numpy.* leakage", test_no_numpy_leakage)
    check("point index stripped", test_point_index_stripped)
    check("keyword mangling (R3)", test_keyword_mangling)
    check("keyword guard is general", test_keyword_scan_covers_every_config)
    check("integer pow lowering (R1)", test_integer_pow_lowering)
    check("max/min are binary", test_max_min_are_binary)
    check("rational prints float pair (R2)", test_rational_prints_float_pair)
    check("grad2 mixed indices (015339c)", test_grad2_mixed_indices_do_not_collapse)
    check("first derivative form", test_first_derivative_form)
    check("agrees with C printer mod rules", test_matches_c_printer_where_it_should)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        sys.exit(1)
    print("all jax_emit checks passed")
