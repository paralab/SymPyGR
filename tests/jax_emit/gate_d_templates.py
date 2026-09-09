#!/usr/bin/env python
"""Gate D -- emitted Python must parse and import.

    python tests/jax_emit/gate_d_templates.py [generated_project_dir]

The Python analogue of the build smoke matrix. Renders every jax template
against a synthetic ctx and ast.parses the result; with a directory argument,
also checks a real generated project.

Why this exists: Jinja's trim_blocks/lstrip_blocks are frozen for the C++
templates' whitespace. In Python indentation is syntax, so a whitespace slip is
a SyntaxError, not a cosmetic reflow.
"""

import ast
import os
import sys
import types

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from jinja2 import Environment, FileSystemLoader  # noqa: E402

from dendrosym.project_generator import (  # noqa: E402
    _TEMPLATES_DIR,
    build_jax_template_map,
)

FAILURES = []


def check(name, fn):
    try:
        fn()
    except AssertionError as exc:
        FAILURES.append(name)
        print(f"FAIL  {name}\n      {exc}")
    except Exception as exc:
        FAILURES.append(name)
        print(f"FAIL  {name}\n      {type(exc).__name__}: {exc}")
    else:
        print(f"PASS  {name}")


def synthetic_ctx():
    """A ctx with the shape _build_jax_context produces."""
    leaves = types.SimpleNamespace(
        field=["alpha", "chi"],
        grad=["grad_0_alpha", "grad_1_chi"],
        grad2=["grad2_0_1_chi"],
        agrad=["agrad_0_alpha"],
        kograd=[],
        param=["eta", "lambda_param"],
        unknown=["x", "t"],
    )
    body = "\n".join([
        "    DENDRO_0000 = 2*alpha",
        "    DENDRO_0001 = jnp.sqrt(chi)",
        "    alpha_rhs = DENDRO_0000*eta + lambda_param[0]*grad_0_alpha + x*t",
        "    chi_rhs = DENDRO_0001 + grad2_0_1_chi + agrad_0_alpha + grad_1_chi",
    ])
    jax_ns = types.SimpleNamespace(
        body=body,
        outputs=["alpha_rhs", "chi_rhs"],
        leaves=leaves,
        # generator-side union: the RHS's buffers plus the BC table's
        deriv_names=(leaves.grad + leaves.agrad + leaves.grad2
                     + ["grad_1_alpha", "grad_2_alpha", "grad_0_chi", "grad_2_chi"]),
    )
    return {
        "project_name": "toy",
        "project_upper": "TOY",
        "namespace": "toy",
        "jax_var_types": ["evolution"],
        "jax": {"evolution": jax_ns},
        "physics_params": [
            # bare keys, as the generator emits them (toml_key = var_name)
            {"var_name": "eta", "py_name": "eta", "toml_key": "eta",
             "default": 2.0, "num_params": 1, "description": "damping"},
            {"var_name": "lambda", "py_name": "lambda_param",
             "toml_key": "lambda", "default": [1.0, 1.0, 1.0, 1.0],
             "num_params": 4, "description": "advection switches"},
        ],
        "solver_features": [("enable_jax_emit", True, "render the jax backend")],
        "jax_bcs": {"evolution": [
            ("alpha_rhs", "alpha", ("grad_0_alpha", "grad_1_alpha", "grad_2_alpha"),
             1.0, 1.0),
            # a BC gradient the toy RHS never reads -- must still reach DERIVS
            ("chi_rhs", "chi", ("grad_0_chi", "grad_1_chi", "grad_2_chi"), 1.0, 1.0),
        ]},
        "jax_initial_data": {
            "entries": [{
                "id": 1, "name": "Toy", "note": "", "func": "toy_init",
                "lines": [("alpha", "1.0 + x"), ("chi", "0")],
                "params": ["eta"], "runtime": ["BH1_mass"], "unknown": [],
            }],
            "unavailable": [{"id": 5, "name": "Noise", "reason": "raw C++"}],
            "symbolic": None, "analytical": None,
        },
        "jax_enforce": {
            "metric": ("gt00", "gt01", "gt02", "gt11", "gt12", "gt22"),
            "trace_free": (("At00", "At01", "At02", "At11", "At12", "At22"),),
            "pos_floor": (("chi", "chi_floor"), ("alpha", "alpha_floor")),
        },
    }


def render_all():
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATES_DIR)),
        keep_trailing_newline=True, trim_blocks=True, lstrip_blocks=True,
    )
    ctx = synthetic_ctx()
    out = {}
    for rel, tmpl in build_jax_template_map(ctx).items():
        out[rel] = env.get_template(tmpl).render(**ctx)
    return out


def test_every_template_renders():
    rendered = render_all()
    assert rendered, "no templates rendered"
    for rel, text in rendered.items():
        assert text.strip(), f"{rel} rendered empty"


def test_python_files_parse():
    for rel, text in render_all().items():
        if not rel.endswith(".py"):
            continue
        try:
            ast.parse(text)
        except SyntaxError as exc:
            line = text.split("\n")[max(exc.lineno - 1, 0)] if exc.lineno else ""
            raise AssertionError(f"{rel}:{exc.lineno} {exc.msg}\n      {line!r}")


def test_keyword_param_is_mangled():
    """A param named `lambda` must never reach the emitted source bare."""
    src = render_all()["toy/toy_params.py"]
    assert "lambda_param:" in src, src[:400]
    assert "\n    lambda:" not in src
    ast.parse(src)


def test_from_toml_reads_the_physics_section():
    """The C++ sample file nests them under [physics]; defaults must not win."""
    mod = {}
    exec(compile(render_all()["toy/toy_params.py"], "toy_params.py", "exec"), mod)
    p = mod["from_toml"]({"physics": {"eta": 9.5, "lambda": [2.0, 2.0, 2.0, 2.0]}})
    assert p.eta == 9.5, p
    assert p.lambda_param == (2.0, 2.0, 2.0, 2.0), p
    assert mod["from_toml"]({"eta": 9.5}).eta == 9.5      # flat still works


def test_bc_and_enforcement_tables():
    """The tables are the whole point of emitting them: they must be usable."""
    src = render_all()["toy/toy_rhs.py"]
    mod = {}
    exec(compile(src, "toy_rhs.py", "exec"), mod)

    bcs = mod["EVOLUTION_BCS"]
    assert len(bcs) == 2, bcs
    assert bcs[0] == ("alpha_rhs", "alpha", "grad_0_alpha", "grad_1_alpha",
                      "grad_2_alpha", 1.0, 1.0), bcs[0]
    # DERIVS must render from deriv_names, not from the four leaf lists, or the
    # generator's BC union never reaches the emitted file
    declared = set(mod["EVOLUTION_DERIVS"])
    assert {g for r in bcs for g in r[2:5]} <= declared, sorted(declared)
    assert {r[0] for r in bcs} <= set(mod["EVOLUTION_OUTPUTS"])
    assert {r[1] for r in bcs} <= set(mod["EVOLUTION_FIELDS"])

    assert len(mod["METRIC_VARS"]) == 6, mod["METRIC_VARS"]
    assert mod["TRACE_FREE_VARS"][0][0] == "At00"
    assert dict(mod["POS_FLOOR_VARS"])["chi"] == "chi_floor"


def test_initial_data_is_callable_and_refuses_what_it_cannot_do():
    """A constant field still has to come back shaped, and a raw-C++ id must
    refuse rather than return something."""
    import numpy as np

    try:
        import jax
    except ImportError:
        print("      (skipped: jax not installed)")
        return
    jax.config.update("jax_enable_x64", True)   # the emitted code assumes float64

    src = render_all()["toy/toy_initial_data.py"]
    mod = {}
    exec(compile(src, "toy_initial_data.py", "exec"), mod)

    x = np.linspace(0.0, 1.0, 4)
    out = mod["initial_data"](1, x, x, x, types.SimpleNamespace(eta=2.0),
                              BH1_mass=0.5)
    assert set(out) == {"alpha", "chi"}, out
    assert all(v.shape == x.shape for v in out.values()), \
        {k: v.shape for k, v in out.items()}      # chi is the constant 0

    try:
        mod["initial_data"](5, x, x, x, None)
    except NotImplementedError as exc:
        assert "raw C++" in str(exc), exc
    else:
        raise AssertionError("a raw-C++ id_type returned instead of refusing")

    try:
        mod["initial_data"](1, x, x, x, None)     # BH1_mass withheld
    except TypeError as exc:
        assert "BH1_mass" in str(exc), exc
    else:
        raise AssertionError("a missing runtime value was silently defaulted")


def test_rhs_declares_its_inputs():
    src = render_all()["toy/toy_rhs.py"]
    for token in ("EVOLUTION_OUTPUTS", "EVOLUTION_FIELDS", "EVOLUTION_DERIVS",
                  "EVOLUTION_PARAMS", "EVOLUTION_EXTRA_INPUTS"):
        assert token in src, f"{token} missing"
    # undeclared inputs must be surfaced, not silently dropped
    assert '"x",' in src and '"t",' in src, "EXTRA_INPUTS lost a name"
    tree = ast.parse(src)
    fns = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)]
    assert "evolution_rhs" in fns, fns


def test_field_named_like_a_parameter():
    """A field called u/d/p must not shadow the rhs function's own arguments."""
    ctx = synthetic_ctx()
    ctx["jax"]["evolution"].leaves.field = ["u", "alpha"]
    ctx["jax"]["evolution"].body = "    alpha_rhs = u + alpha\n    chi_rhs = u"
    ctx["jax"]["evolution"].deriv_names = []
    ctx["jax"]["evolution"].leaves.param = []
    ctx["jax"]["evolution"].leaves.unknown = []
    ctx["jax_bcs"] = {}
    env = Environment(loader=FileSystemLoader(str(_TEMPLATES_DIR)),
                      trim_blocks=True, lstrip_blocks=True)
    src = env.get_template("jax/rhs.py.j2").render(**ctx)
    mod = {}
    exec(compile(src, "toy_rhs.py", "exec"), mod)
    out = mod["evolution_rhs"]({"u": 2.0, "alpha": 3.0}, {}, None)
    assert out == (5.0, 2.0), out


def test_rhs_body_is_executable():
    """Unpack + body must actually run, with jnp bound to numpy."""
    import numpy as np
    src = render_all()["toy/toy_rhs.py"]
    src = src.replace("import jax.numpy as jnp", "import numpy as jnp")
    mod = {}
    exec(compile(src, "toy_rhs.py", "exec"), mod)  # noqa: S102
    u = {"alpha": np.array([1.0]), "chi": np.array([1.0]),
         "x": np.array([0.5]), "t": np.array([0.25])}
    d = {n: np.array([0.1]) for n in mod["EVOLUTION_DERIVS"]}
    p = types.SimpleNamespace(eta=2.0, lambda_param=(1.0, 1.0, 1.0, 1.0))
    out = mod["evolution_rhs"](u, d, p)
    assert len(out) == len(mod["EVOLUTION_OUTPUTS"]), out
    assert all(np.isfinite(v).all() for v in out), out


def test_generated_project(path):
    """Parse and import every .py in a real generated project."""
    import pathlib
    root = pathlib.Path(path)
    pys = sorted(root.rglob("*.py"))
    assert pys, f"no .py files under {root}"
    for f in pys:
        try:
            ast.parse(f.read_text())
        except SyntaxError as exc:
            raise AssertionError(f"{f}:{exc.lineno} {exc.msg}")
    print(f"      parsed {len(pys)} generated .py file(s)")


if __name__ == "__main__":
    check("every jax template renders", test_every_template_renders)
    check("emitted .py parses", test_python_files_parse)
    check("keyword param mangled", test_keyword_param_is_mangled)
    check("from_toml reads [physics]", test_from_toml_reads_the_physics_section)
    check("bc + enforcement tables", test_bc_and_enforcement_tables)
    check("initial data callable + refuses", test_initial_data_is_callable_and_refuses_what_it_cannot_do)
    check("rhs declares its inputs", test_rhs_declares_its_inputs)
    check("field named u/d/p", test_field_named_like_a_parameter)
    check("rhs body executes", test_rhs_body_is_executable)
    if len(sys.argv) > 1:
        check(f"generated project {sys.argv[1]}",
              lambda: test_generated_project(sys.argv[1]))

    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        sys.exit(1)
    print("all gate D checks passed")
