#!/usr/bin/env python
"""run_all.py -- LTS codegen checks.

    python tests/lts/run_all.py            # fast checks + BSSN known answer
    python tests/lts/run_all.py --ccz4     # also audit CCZ4 from ../ccz4-gr
    python tests/lts/run_all.py --fast     # skip the BSSN equations (~1 min)

Plain asserts, no pytest. Each check asserts both the finding and its absence on
a clean input, so it fails if the piece it covers stops detecting (or starts
over-reporting).
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, ROOT)

import sympy as sym  # noqa: E402

from dendrosym import lts_audit  # noqa: E402

FAILURES = []


def check(name, fn):
    try:
        fn()
    except AssertionError as exc:
        FAILURES.append(name)
        print(f"FAIL  {name}\n      {exc}")
    else:
        print(f"PASS  {name}")


class _StubConfig:
    """Just enough of a DendroConfiguration for the audit."""

    idx_str = "[pp]"

    def __init__(self, rhs, constraints=()):
        self._rhs = rhs
        self.all_var_names = {"constraint": list(constraints)}

    def _extract_rhs_expressions(self, var_type, append_rhs_to_var=True):
        names = list(self._rhs)
        return [self._rhs[n] for n in names], names, [], [], len(names)


def _sym(name):
    return sym.Symbol(f"{name}[pp]")


grad = sym.Function("grad")
k, dt = sym.Symbol("k"), sym.Symbol("dt")
u, v, w, ham = _sym("u"), _sym("v"), _sym("w"), _sym("ham")


def test_stiff_rate_found():
    a = lts_audit.audit(_StubConfig({"u": -k * u + grad(0, v), "v": grad(1, u)}))
    assert [r.var for r in a.stiff_rates] == ["u"], a.summary()
    assert a.stiff_rates[0].expr == -k, a.stiff_rates[0].expr
    assert a.rate_params == ["k"]


def test_derivative_terms_are_not_rates():
    a = lts_audit.audit(_StubConfig({"u": k * grad(0, u), "v": u * k}))
    assert not a.stiff_rates, a.summary()


def test_constraint_read_found():
    a = lts_audit.audit(
        _StubConfig({"u": -k * u + ham, "v": grad(0, ham), "w": u}, ["ham"])
    )
    assert sorted(a.constraint_reads) == ["u", "v"], a.constraint_reads
    assert not a.clean
    clean = lts_audit.audit(_StubConfig({"u": -k * u + v, "w": u}, ["ham"]))
    assert clean.clean and not clean.constraint_reads, clean.summary()


def test_dt_term_found():
    a = lts_audit.audit(_StubConfig({"u": dt * u + v, "v": u}))
    assert sorted(a.dt_terms) == ["u"], a.dt_terms
    assert not a.clean
    assert lts_audit.audit(_StubConfig({"u": k * u})).clean


def test_rate_cpp_names_and_magnitude():
    a = lts_audit.audit(_StubConfig({"u": -2 * k * u}))
    code = lts_audit.rate_cpp(a.stiff_rates[0], lambda n: f"NS::P_{n.upper()}")
    assert code == "std::fabs(-2*NS::P_K)", code


def test_bssn_eta_is_the_known_stiff_rate():
    sys.path.insert(0, os.path.join(ROOT, "tests", "jax_emit"))
    import bssn_eqns

    a = lts_audit.audit(bssn_eqns.dendroConfigs)
    assert [r.var for r in a.stiff_rates] == ["B0", "B1", "B2"], a.summary()
    assert a.rate_params == ["eta"], a.summary()
    assert all(r.expr == -sym.Symbol("eta") for r in a.stiff_rates)
    assert a.clean, a.summary()


def test_ccz4_audit():
    sys.path.insert(0, os.path.join(ROOT, "..", "ccz4-gr"))
    import ccz4

    a = lts_audit.audit(ccz4.dendroConfigs)
    assert a.clean, a.summary()
    assert not a.indeterminate, a.indeterminate
    assert {"eta", "kappa[0]"} <= set(a.rate_params), a.rate_params


check("audit finds a stiff rate", test_stiff_rate_found)
check("audit ignores derivative terms", test_derivative_terms_are_not_rates)
check("audit finds constraint reads", test_constraint_read_found)
check("audit finds dt terms", test_dt_term_found)
check("rate_cpp prints names and magnitude", test_rate_cpp_names_and_magnitude)
if "--fast" not in sys.argv:
    check("BSSN: eta in B is the stiff rate", test_bssn_eta_is_the_known_stiff_rate)
if "--ccz4" in sys.argv:
    check("CCZ4 audit: eta and kappa1 rates, no reads", test_ccz4_audit)

if FAILURES:
    print(f"\n{len(FAILURES)} failed: {', '.join(FAILURES)}")
    sys.exit(1)
print("\nall LTS checks passed")
