"""jnp printer -- the JAX twin of :mod:`dendrosym.code_printer`.

Rules derived from a statement-level diff of DendroJAX's hand-written BSSN RHS
against dendrogr_dfvk/CodeGen/bssneqs_SSL_HD_dxsq.cpp (850 statements, zero
residue).

Fields arrive as ``Symbol("alpha[pp]")``, so :meth:`_print_Symbol` is the live
mechanism. (DendroCPrinter's per-field ``known_functions`` map is dead code:
general_configs rebinds ``derivs.variable_strs`` after code_printer captured
it.) ``idx_str`` is read through the module so this cannot go stale the same way.
"""

import keyword

import dendrosym.derivs
from sympy.printing.numpy import NumPyPrinter
from sympy.printing.precedence import precedence

# `lambda` is declared by bssn_eqns.py, ccz4_eqns.py and emda_configs.py.
KEYWORD_SUFFIX = "_param"


def safe_name(name):
    """Mangle `name` if it is a Python keyword."""
    return name + KEYWORD_SUFFIX if keyword.iskeyword(name) else name


def jax_symbol_name(name, idx_str=None):
    """``alpha[pp]`` -> ``alpha``; ``lambda[0]`` -> ``lambda_param[0]``.

    Strips the point index (in JAX the array is the value) and keyword-mangles
    the identifier while keeping a real subscript.
    """
    if idx_str is None:
        idx_str = dendrosym.derivs.idx_str
    if idx_str and name.endswith(idx_str):
        name = name[: -len(idx_str)]
    if name.endswith("]") and "[" in name:      # a parameter subscript, not [pp]
        base, sub = name[: name.index("[")], name[name.index("["):]
        return safe_name(base) + sub
    return safe_name(name)


class DendroJaxPrinter(NumPyPrinter):
    """Prints Dendro expressions as ``jax.numpy`` source.

    Emits ``grad``/``grad2`` in call form so codegen_jax's rename handles them,
    and lowers integer powers like the C99 printer.
    """

    _module = "jnp"

    def __init__(self, additional_user_funcs={}, fields=None):
        super().__init__()
        # 29 inherited entries hardcode a "numpy." prefix, so _module alone
        # would still emit numpy.exp into a traced kernel.
        self.known_functions = {
            k: (self._module + v[len("numpy"):]
                if isinstance(v, str) and v.startswith("numpy.") else v)
            for k, v in self.known_functions.items()
        }
        self.known_functions = dict(self.known_functions, **additional_user_funcs)
        self._fields = set(
            fields if fields is not None else dendrosym.derivs.variable_strs
        )

    def _print_Symbol(self, expr):
        return jax_symbol_name(expr.name)

    def _print_Function(self, expr):
        """``alpha(xx, yy, zz)`` is the field alpha, so print it bare.

        Matched on the position signature, since the registered-name list is
        unreliable -- which is why the C printer raises here.
        """
        pos = (dendrosym.derivs.xx, dendrosym.derivs.yy, dendrosym.derivs.zz)
        name = getattr(expr.func, "__name__", "")
        if name not in self.known_functions and (
            expr.args == pos or name in self._fields
        ):
            return jax_symbol_name(name)
        return super()._print_Function(expr)

    def _print_Float(self, flt):
        """17 digits like the C printer; the default drops a ULP."""
        num = str(flt.evalf(17))
        if "e" not in num and "." not in num:
            num += ".0"
        head, *tail = num.split("e")
        head = head.rstrip("0")
        if head.endswith("."):
            head += "0"
        return "e".join([head] + tail)

    def _print_im(self, expr):
        return "0.0"        # GR quantities are real; sympy just can't prove it

    def _print_re(self, expr):
        return self._print(expr.args[0])

    def _print_Rational(self, expr):
        return f"{float(expr.p)}/{float(expr.q)}"

    def _print_Max(self, expr):
        return self._nested_binary("maximum", expr.args)

    def _print_Min(self, expr):
        return self._nested_binary("minimum", expr.args)

    def _nested_binary(self, fn, args):
        """n-ary Max/Min as nested 2-arg calls; the inherited form allocates."""
        call = self._module_format(f"{self._module}.{fn}")
        out = self._print(args[0])
        for a in args[1:]:
            out = f"{call}({out}, {self._print(a)})"
        return out

    def _print_Pow(self, expr):
        """Integer powers as multiply/reciprocal chains, never ``**``."""
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
        """``grad(i, f)`` / ``grad2(i, j, f)``; codegen_jax renames them."""
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
            if idx1 > idx2:                 # keep the 015339c ordering fix
                idx1, idx2 = idx2, idx1
            return "grad2(%d, %d, %s)" % (idx1, idx2, self._print(differand))

        ((wrt, _count),) = wrt_counts
        return "grad(%d, %s)" % (d_order.index(wrt), self._print(differand))
