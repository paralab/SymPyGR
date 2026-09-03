"""jnp-targeting printer -- the JAX twin of :mod:`dendrosym.code_printer`.

The emission rules here were derived, not chosen. DendroJAX's hand-written BSSN
RHS (``examples/bssn_solver.py``) is a transliteration of
``dendrogr_dfvk/CodeGen/bssneqs_SSL_HD_dxsq.cpp``; a statement-level diff of the
two -- 850 statements, identical LHS sequence -- reduces to five mechanical
rules with zero residue. Four are implemented here; the fifth (grad -> agrad) is
a config flag, not a printer concern.

Note on where the work happens: by print time fields are plain
``Symbol("alpha[pp]")``, not ``Function``. :class:`DendroCPrinter` also carries a
``known_functions`` entry per field name, but that mapping never fires --
``general_configs`` *rebinds* ``dendrosym.derivs.variable_strs`` while
``code_printer`` captured the original empty list at import, so the C printer's
copy is permanently ``[]``. The live mechanism is :meth:`_print_Symbol`, and
``idx_str`` is read through the module on every call so this printer cannot
inherit the same staleness.
"""

import keyword

import dendrosym.derivs
from sympy.printing.numpy import NumPyPrinter
from sympy.printing.precedence import precedence

# Suffix appended to a declared name that collides with a Python keyword.
# Not hypothetical: `lambda` is declared by bssn_eqns.py, ccz4_eqns.py AND
# emda_configs.py, so an unmangled `lambda[0]` is a SyntaxError, not a warning.
KEYWORD_SUFFIX = "_param"


def safe_name(name):
    """Mangle `name` if it would be a Python keyword. Leaves everything else."""
    return name + KEYWORD_SUFFIX if keyword.iskeyword(name) else name


def jax_symbol_name(name, idx_str=None):
    """Turn a Dendro symbol name into a valid Python identifier reference.

    Strips the point index (``alpha[pp]`` -> ``alpha``; in JAX the whole array is
    the value) and keyword-mangles the identifier while preserving a real
    subscript (``lambda[0]`` -> ``lambda_param[0]``).
    """
    if idx_str is None:
        idx_str = dendrosym.derivs.idx_str
    if idx_str and name.endswith(idx_str):
        name = name[: -len(idx_str)]
    # a surviving trailing [..] is a genuine parameter subscript, not the point
    # index -- mangle only the identifier in front of it.
    if name.endswith("]") and "[" in name:
        base, sub = name[: name.index("[")], name[name.index("["):]
        return safe_name(base) + sub
    return safe_name(name)


class DendroJaxPrinter(NumPyPrinter):
    """Prints Dendro expressions as ``jax.numpy`` source.

    Mirrors :class:`dendrosym.code_printer.DendroCPrinter` where it can:
    ``grad``/``grad2`` are emitted in call form so
    :func:`dendrosym.codegen.change_deriv_names` -- string-level and
    backend-agnostic -- rewrites them unchanged, and integer powers get the same
    multiply/reciprocal lowering.

    Where it must differ:

    * symbols drop the point index and are keyword-mangled
      (:func:`jax_symbol_name`).
    * ``Rational`` prints ``3.0/4.0`` rather than ``(3/4)``, matching the
      reference and keeping the literal unambiguously float.
    * ``Max``/``Min`` print ``jnp.maximum``/``jnp.minimum``. NumPyPrinter's
      default is ``amax(asarray([a, b]), axis=0)``, which materialises an array
      per call -- unusable in a per-point kernel.
    * every inherited ``numpy.*`` target is retargeted to ``jnp.*``. Setting
      ``_module`` alone is not enough: 29 entries in NumPyPrinter's function
      table carry a hardcoded ``numpy.`` prefix, so ``exp`` would otherwise
      emit ``numpy.exp`` into a traced kernel.
    """

    _module = "jnp"

    def __init__(self, additional_user_funcs={}, fields=None):
        super().__init__()

        # retarget the hardcoded numpy.* table (exp, log, sin, ... 29 entries)
        self.known_functions = {
            k: (self._module + v[len("numpy"):]
                if isinstance(v, str) and v.startswith("numpy.") else v)
            for k, v in self.known_functions.items()
        }
        self.known_functions = dict(self.known_functions, **additional_user_funcs)

        # Registered field names. Passed explicitly by the emitter; the module
        # global is only a fallback and is usually stale (see module docstring).
        self._fields = set(
            fields if fields is not None else dendrosym.derivs.variable_strs
        )

    # ------------------------------------------------------------------
    # leaves
    # ------------------------------------------------------------------
    def _print_Symbol(self, expr):
        return jax_symbol_name(expr.name)

    def _print_Function(self, expr):
        """A field as a function of position prints as the bare array name.

        ``alpha(xx, yy, zz)`` *is* the field alpha, so it prints as ``alpha``.
        Matched on the (xx, yy, zz) signature rather than on the registered-name
        list, because that list cannot be relied on -- which is why the C
        printer raises here instead. Anything else defers to the base printer.
        """
        pos = (dendrosym.derivs.xx, dendrosym.derivs.yy, dendrosym.derivs.zz)
        name = getattr(expr.func, "__name__", "")
        if name not in self.known_functions and (
            expr.args == pos or name in self._fields
        ):
            return jax_symbol_name(name)
        return super()._print_Function(expr)

    def _print_im(self, expr):
        # all GR quantities are real; sympy emits im() only when it cannot prove
        # this (e.g. psi**p_expo with unconstrained p_expo). Mirrors the C printer.
        return "0.0"

    def _print_re(self, expr):
        return self._print(expr.args[0])

    def _print_Rational(self, expr):
        return f"{float(expr.p)}/{float(expr.q)}"

    # ------------------------------------------------------------------
    # operators
    # ------------------------------------------------------------------
    def _print_Max(self, expr):
        return self._nested_binary("maximum", expr.args)

    def _print_Min(self, expr):
        return self._nested_binary("minimum", expr.args)

    def _nested_binary(self, fn, args):
        """Fold an n-ary Max/Min into nested two-arg jnp calls."""
        call = self._module_format(f"{self._module}.{fn}")
        out = self._print(args[0])
        for a in args[1:]:
            out = f"{call}({out}, {self._print(a)})"
        return out

    def _print_Pow(self, expr):
        """Integer powers become multiply / reciprocal chains, never ``**``.

        Same lowering as the C99 printer, for the same reason: keeps the
        arithmetic explicit and avoids a pow call for small exponents. The
        negative branch matters -- ``pow(x, -2)`` is the one the cascade found
        surviving as a real libm call in the compiled C++ object.
        """
        PREC = precedence(expr)
        if expr.exp in range(2, 7):
            inner = "*".join([self.parenthesize(expr.base, PREC)] * int(expr.exp))
            return "(" + inner + ")"
        elif expr.exp in range(-6, 0):
            inner = "*".join([self.parenthesize(expr.base, PREC)] * int(-expr.exp))
            return "(1.0/(" + inner + "))"
        else:
            return super()._print_Pow(expr)

    def _print_Derivative(self, expr):
        """Emit ``grad(i, f)`` / ``grad2(i, j, f)``; codegen rewrites the names.

        Identical index logic to the C printer, including the mixed-index swap
        fixed in 015339c -- do not collapse ``idx2`` onto ``idx1``.
        """
        differand, *wrt_counts = expr.args
        d_order = (dendrosym.derivs.xx, dendrosym.derivs.yy, dendrosym.derivs.zz)

        if len(wrt_counts) > 2 or wrt_counts[0][1] > 2:
            raise NotImplementedError(
                "Currently only first and second order derivs are supported!"
            )

        if wrt_counts[0][1] == 2:
            ((wrt, _count),) = wrt_counts
            idx = d_order.index(wrt)
            return "grad2(%d, %d, %s)" % (idx, idx, self._print(differand))

        if len(wrt_counts) == 2:
            (wrt1, _c1), (wrt2, _c2) = wrt_counts
            idx1, idx2 = d_order.index(wrt1), d_order.index(wrt2)
            if idx1 > idx2:
                idx1, idx2 = idx2, idx1
            return "grad2(%d, %d, %s)" % (idx1, idx2, self._print(differand))

        ((wrt, _count),) = wrt_counts
        return "grad(%d, %s)" % (d_order.index(wrt), self._print(differand))
