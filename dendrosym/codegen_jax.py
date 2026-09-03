"""JAX emission -- the jnp twin of the scalar C++ path in :mod:`dendrosym.codegen`.

The important property: this consumes the **same** ``cse_list`` the C++ path
consumes. ``construct_expression_list`` / ``construct_cse_from_list`` are
backend-agnostic, so the ``DENDRO_n`` numbering and statement order are identical
between the two backends by construction. The physics is never re-derived and
never re-CSE'd for JAX.

Two string passes the C++ path runs are deliberately skipped here:
``apply_input_struct`` (there is no ``in.`` struct in JAX -- the array is the
value) and the output-struct rewrite. The derivative rename is *reimplemented*
rather than reused: ``codegen.change_deriv_names`` requires a literal ``[pp]``
in both of its patterns, which this backend has stripped by then, so calling it
would silently leave ``grad(0, alpha)`` in the output as a call to a
nonexistent function.
"""

import heapq
import re
from typing import List, NamedTuple, Tuple

import sympy as sym

from dendrosym.jax_printer import DendroJaxPrinter, jax_symbol_name

# the derivative call forms the printer emits
DERIV_FUNCS = {"grad": "grad", "grad2": "grad2", "agrad": "agrad", "kograd": "kograd"}

# codegen.change_deriv_names hardcodes `\[pp\]` in both patterns, so it matches
# nothing here -- the JAX printer strips the point index. Same rewrite, index
# optional. Kept local rather than loosening the shared C patterns, which are
# under the byte-identity gate.
_DERIV1_JAX = re.compile(r"\b(agrad|grad|kograd)\((\d),\s*(\w+)\)")
_DERIV2_JAX = re.compile(r"\bgrad2\((\d),\s*(\d),\s*(\w+)\)")


def _d1(m):
    return f"{m.group(1)}_{m.group(2)}_{m.group(3)}"


def _d2(m):
    a, b = int(m.group(1)), int(m.group(2))
    if a > b:                       # keep the 015339c ordering fix
        a, b = b, a
    return f"grad2_{a}_{b}_{m.group(3)}"


def change_deriv_names_jax(src: str) -> str:
    """``grad(i, var)`` -> ``grad_i_var``; ``grad2(i, j, var)`` -> ``grad2_i_j_var``."""
    return _DERIV2_JAX.sub(_d2, _DERIV1_JAX.sub(_d1, src))


class JaxBody(NamedTuple):
    """An emitted RHS body, kept structured so gates can evaluate it.

    ``statements`` is every assignment in dependency order -- CSE temporaries
    then outputs -- as ``(lhs identifier, rhs source)``. ``outputs`` names the
    subset that are RHS outputs, in registered order. ``exprs`` is the sympy
    expression each statement was printed *from*, in the same order: a gate that
    re-parses the emitted text is checking the text against itself, so the only
    honest differential compares the text against these.
    """

    statements: List[Tuple[str, str]]
    outputs: List[str]
    exprs: List[sym.Expr] = ()

    def render(self, indent: str = "") -> str:
        return "\n".join(f"{indent}{lhs} = {rhs}" for lhs, rhs in self.statements)

    @property
    def n_temps(self) -> int:
        return len(self.statements) - len(self.outputs)


def _identifier(name) -> str:
    """Normalise an emitted assignment target to a bare Python identifier.

    Strips a C declaration (``double x``), a struct qualifier (``out.alpha``)
    and the point index, then keyword-mangles what is left.
    """
    name = str(name).strip()
    if " " in name:                       # "double DENDRO_STAGED_VAR_0"
        name = name.split()[-1]
    if "." in name:                       # "out.alpha"
        name = name.rsplit(".", 1)[-1]
    return jax_symbol_name(name)


def atomize_derivs(expr, printer=None):
    """Replace ``grad``/``grad2``/``agrad``/``kograd`` applications with Symbols
    named exactly as the emitter names them.

    Lets a differential gate evaluate the ORIGINAL sympy tree using the same leaf
    names the emitted source reads, so the two sides share only the leaf values
    -- not the printed text.
    """
    if printer is None:
        printer = DendroJaxPrinter(additional_user_funcs=DERIV_FUNCS)
    repl = {}
    for app in expr.atoms(sym.Function):
        if getattr(app.func, "__name__", "") in DERIV_FUNCS:
            repl[app] = sym.Symbol(change_deriv_names_jax(printer.doprint(app)))
    return expr.xreplace(repl) if repl else expr


def build_jax_body(
    cse_list,
    rhs_var_names,
    fields=None,
    interleave_outputs: bool = False,
) -> JaxBody:
    """Turn a ``(cse_temps, output_exprs)`` pair into JAX statements.

    ``interleave_outputs`` mirrors the C++ path: a staged block defines its
    quantities in terms of each other, so temporaries can reference outputs and
    the temps-then-outputs layout would emit a read above its definition. Emit
    in dependency order instead.
    """
    printer = DendroJaxPrinter(additional_user_funcs=DERIV_FUNCS, fields=fields)

    def rhs_src(expr) -> str:
        return change_deriv_names_jax(printer.doprint(expr))

    out_names = [_identifier(n) for n in rhs_var_names]

    if not interleave_outputs:
        statements, exprs = [], []
        for v1, v2 in cse_list[0]:
            statements.append((_identifier(v1), rhs_src(v2)))
            exprs.append(v2)
        for i, e in enumerate(cse_list[1]):
            statements.append((out_names[i], rhs_src(e)))
            exprs.append(e)
        return JaxBody(statements, out_names, exprs)

    # dependency-ordered emission (same algorithm as generate_cpu_preextracted)
    stmts = []  # (defined symbol, lhs, rhs source, expr, referenced symbols)
    for v1, v2 in cse_list[0]:
        stmts.append((v1, _identifier(v1), rhs_src(v2), v2, v2.free_symbols))
    for i, e in enumerate(cse_list[1]):
        bare = str(rhs_var_names[i]).split()[-1]
        stmts.append((sym.Symbol(bare), out_names[i], rhs_src(e), e, e.free_symbols))

    defined = {}
    for k, (dsym, _l, _t, _e, _r) in enumerate(stmts):
        if dsym is not None:
            defined.setdefault(dsym, k)
    children = [[] for _ in stmts]
    indeg = [0] * len(stmts)
    for k, (_d, _l, _t, _e, refs) in enumerate(stmts):
        deps = {defined[r] for r in refs if r in defined and defined[r] != k}
        indeg[k] = len(deps)
        for j in deps:
            children[j].append(k)

    ready = [k for k in range(len(stmts)) if indeg[k] == 0]
    heapq.heapify(ready)
    ordered, ordered_exprs = [], []
    while ready:
        k = heapq.heappop(ready)
        ordered.append((stmts[k][1], stmts[k][2]))
        ordered_exprs.append(stmts[k][3])
        for c in children[k]:
            indeg[c] -= 1
            if indeg[c] == 0:
                heapq.heappush(ready, c)
    if len(ordered) != len(stmts):
        stuck = [stmts[k][1] for k in range(len(stmts)) if indeg[k] > 0]
        raise ValueError(
            "cyclic dependency among the emitted statements; "
            f"{len(stmts) - len(ordered)} unresolved, e.g. {stuck[:8]}"
        )
    return JaxBody(ordered, out_names, ordered_exprs)


def generate_jax_preextracted(
    cse_list,
    rhs_var_names,
    fields=None,
    interleave_outputs: bool = False,
    return_stats: bool = False,
):
    """Text-emitting wrapper, shaped like ``generate_cpu_preextracted``."""
    body = build_jax_body(
        cse_list, rhs_var_names, fields=fields, interleave_outputs=interleave_outputs
    )
    text = body.render()
    if not return_stats:
        return text
    reduced_ops = sum(sym.count_ops(e) for e in cse_list[1])
    reduced_ops += sum(sym.count_ops(e) for _s, e in cse_list[0])
    return text, reduced_ops
