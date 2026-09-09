#!/usr/bin/env python
"""Gate E -- the emitted JAX body must be the C++ body, statement for statement.

    python tests/jax_emit/gate_e_cpp_parity.py <config.py> [config.py ...]

Both backends share one cse_list, so the bodies must correspond exactly once the
struct/naming rewrites are undone. Catches the two paths seeing different
*expressions*, not just different text.
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import dendrosym  # noqa: E402
from dendrosym.codegen_jax import change_deriv_names_jax  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from gate_c_fullbody import find_config_object, load_config  # noqa: E402


def cpp_statements(src):
    """`double X = expr;` -> (lhs, rhs)."""
    out = []
    for line in src.splitlines():
        line = line.strip()
        if not line or line.startswith("//"):
            continue
        if not line.endswith(";") or " = " not in line:
            raise AssertionError(f"unparsed C++ statement: {line[:120]}")
        lhs, rhs = line[:-1].split(" = ", 1)
        out.append((lhs.strip(), rhs.strip()))
    return out


def norm_cpp(lhs, rhs, idx_str):
    """C++ statement -> the form the jax body prints."""
    def strip(s):
        s = s.replace(idx_str, "")
        s = re.sub(r"\b(?:in|out)\.", "", s)
        s = re.sub(r"^double\s+", "", s)
        s = s.replace("lambda[", "lambda_param[")
        return s

    lhs = strip(lhs)
    if not lhs.startswith("DENDRO_"):       # out.alpha is the jax alpha_rhs
        lhs += "_rhs"
    rhs = change_deriv_names_jax(strip(rhs))
    # accepted rewrite: x**(-1/2) is 1/sqrt(x) in jnp, 1 ULP apart at worst
    rhs = re.sub(r"pow\(([^,]+), -1\.0/2\.0\)", r"1/sqrt(\1)", rhs)
    return lhs, rhs


def norm_jax(lhs, rhs):
    return lhs, rhs.replace("jnp.", "")


def compare(cfg, var_type):
    cfg.find_derivatives(var_type)
    cpu = cfg.generate_rhs_code(var_type, arc_type="cpu")
    jax = cfg.generate_rhs_code(var_type, arc_type="jax")

    C = cpp_statements(cpu)
    J = list(jax.statements)
    idx_str = dendrosym.derivs.idx_str

    diffs = []
    if len(C) != len(J):
        diffs.append(f"statement count: cpp {len(C)} vs jax {len(J)}")
    for i, (c, j) in enumerate(zip(C, J)):
        a, b = norm_cpp(*c, idx_str), norm_jax(*j)
        if a != b:
            diffs.append(f"stmt {i}\n    cpp: {a[0]} = {a[1][:220]}\n    jax: {b[0]} = {b[1][:220]}")

    # a derivative of a CSE temp is never a buffer the driver can fill
    orphan = sorted({
        m for _lhs, rhs in J
        for m in re.findall(r"\b(?:a?grad|grad2|kograd)(?:_\d)+_(DENDRO_\w+)", rhs)
    })
    if orphan:
        diffs.append(f"deriv buffers named after CSE temps (no driver can fill these): {orphan}")
    return len(C), diffs


def main():
    paths = sys.argv[1:]
    if not paths:
        print("usage: gate_e_cpp_parity.py <config.py> [...]")
        return 2

    failures = 0
    for path in paths:
        cfg = find_config_object(load_config(path))
        name = os.path.basename(path)
        for vt in cfg.all_var_names:
            if cfg.all_rhs_functions.get(vt) is None:
                continue
            n, diffs = compare(cfg, vt)
            if diffs:
                failures += 1
                print(f"FAIL  {name}:{vt}  ({n} statements)")
                for d in diffs[:10]:
                    print("      " + d.replace("\n", "\n      "))
                if len(diffs) > 10:
                    print(f"      ... {len(diffs) - 10} more")
            else:
                print(f"PASS  {name}:{vt}  {n} statements identical")

    print("\nall gate E checks passed" if not failures
          else f"\n{failures} var_type(s) diverge from the C++ backend")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
