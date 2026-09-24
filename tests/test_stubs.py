"""The type stubs describe the extension that is actually built.

The stubs are hand written, so a function added to or renamed in the bindings
would otherwise leave them silently stale: type checkers would reject working
code, and the API reference, which takes its signatures from the stubs, would
document something that does not exist.
"""
import ast
import inspect
from pathlib import Path

from radiological_material_clearance_finder import _core

#: Next to the extension, so this checks the stub a wheel actually ships.
STUB = Path(_core.__file__).with_name("_core.pyi")


def _stub_names():
    tree = ast.parse(STUB.read_text())
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.AnnAssign):
            names.add(node.target.id)
    return names


def _runtime_names():
    return {name for name in vars(_core) if not name.startswith("_") or name == "__version__"}


def test_the_stub_names_exactly_what_the_extension_exports():
    stub, runtime = _stub_names(), _runtime_names()
    assert stub - runtime == set(), "in the stub but not the extension"
    assert runtime - stub == set(), "in the extension but not the stub"


def test_function_parameters_match():
    tree = ast.parse(STUB.read_text())
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        stubbed = [a.arg for a in node.args.posonlyargs + node.args.args + node.args.kwonlyargs]
        actual = list(inspect.signature(getattr(_core, node.name)).parameters)
        assert stubbed == actual, node.name
