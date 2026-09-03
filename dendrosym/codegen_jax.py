"""JAX emission -- the jnp twin of the scalar C++ path in :mod:`dendrosym.codegen`.

Consumes the same ``cse_list`` the C++ path does, so DENDRO_n numbering and
order match between backends and the physics is never re-CSE'd.

Skips the C++ string passes ``apply_input_struct`` and the output-struct
rewrite (JAX has neither). The derivative rename is reimplemented below because
``codegen.change_deriv_names`` requires a literal ``[pp]``.
"""

import heapq
import re
from typing import List, NamedTuple, Tuple

import sympy as sym

from dendrosym.jax_printer import DendroJaxPrinter, jax_symbol_name

DERIV_FUNCS = {"grad": "grad", "grad2": "grad2", "agrad": "agrad", "kograd": "kograd"}

# Same rewrite as codegen.change_deriv_names, index optional. Local rather than
# loosening the shared C patterns, which are under the byte-identity gate.
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
    """Emitted body, kept structured so gates can evaluate it.

    statements: (lhs, rhs source) in dependency order, temps then outputs.
    exprs: the sympy expression each was printed from -- a gate that re-parses
    the emitted text is only checking it against itself.
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
    """``double x`` / ``out.alpha`` / ``alpha[pp]`` -> a bare identifier."""
    name = str(name).strip()
    if " " in name:                       # "double DENDRO_STAGED_VAR_0"
        name = name.split()[-1]
    if "." in name:                       # "out.alpha"
        name = name.rsplit(".", 1)[-1]
    return jax_symbol_name(name)


def atomize_derivs(expr, printer=None):
    """Derivative applications -> Symbols named as the emitter names them.

    Lets a gate evaluate the original sympy tree against the emitted source
    sharing only leaf values, not text.
    """
    if printer is None:
        printer = DendroJaxPrinter(additional_user_funcs=DERIV_FUNCS)
    repl = {}
    for app in expr.atoms(sym.Function):
        if getattr(app.func, "__name__", "") in DERIV_FUNCS:
            repl[app] = sym.Symbol(change_deriv_names_jax(printer.doprint(app)))
    return expr.xreplace(repl) if repl else expr


_DERIV_LEAF = re.compile(r"^(agrad|grad|kograd)_(\d)_(\w+)$")
_DERIV2_LEAF = re.compile(r"^grad2_(\d)_(\d)_(\w+)$")


def classify_leaves(body, field_names=(), param_names=()):
    """Split what a body reads into fields / derivs / params / unknown.

    `unknown` is what matters: anything the equations use but the config never
    declared, surfaced here instead of as a NameError inside a jitted kernel.
    """
    fields = {str(f).split("[")[0] for f in field_names}
    fields = {jax_symbol_name(f) for f in fields}
    params = {jax_symbol_name(str(p)) for p in param_names}

    defined = set()
    read = []
    for (lhs, _rhs), e in zip(body.statements, body.exprs or []):
        for s in atomize_derivs(e).free_symbols:
            n = jax_symbol_name(s.name)
            if n not in defined:
                read.append(n)
        defined.add(lhs)
    read = sorted(set(read))

    out = {"field": [], "grad": [], "grad2": [], "agrad": [], "kograd": [],
           "param": [], "unknown": []}
    for n in read:
        base = n.split("[")[0]
        m2 = _DERIV2_LEAF.match(base)
        m1 = _DERIV_LEAF.match(base)
        if m2:
            out["grad2"].append(n)
        elif m1:
            out[m1.group(1)].append(n)
        elif base in fields:
            out["field"].append(n)
        elif base in params:
            out["param"].append(base)
        else:
            out["unknown"].append(n)
    out["param"] = sorted(set(out["param"]))
    return out


def build_jax_body(
    cse_list,
    rhs_var_names,
    fields=None,
    interleave_outputs: bool = False,
) -> JaxBody:
    """``(cse_temps, output_exprs)`` -> JAX statements.

    interleave_outputs: staged blocks define quantities in terms of each other,
    so emit in dependency order or a temp lands above the output it reads.
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
    """Text wrapper, shaped like ``generate_cpu_preextracted``."""
    body = build_jax_body(
        cse_list, rhs_var_names, fields=fields, interleave_outputs=interleave_outputs
    )
    text = body.render()
    if not return_stats:
        return text
    reduced_ops = sum(sym.count_ops(e) for e in cse_list[1])
    reduced_ops += sum(sym.count_ops(e) for _s, e in cse_list[0])
    return text, reduced_ops
