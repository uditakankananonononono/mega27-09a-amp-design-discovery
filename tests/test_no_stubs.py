"""Repo-wide audit: no stub functions - every function body must do real work."""
import ast
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
SKIP = {"tests"}


def iter_functions():
    for path in ROOT.rglob("*.py"):
        if any(part in SKIP or part.startswith(".") for part in path.parts[len(ROOT.parts):]):
            continue
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                yield path, node


def test_no_stub_bodies():
    bad = []
    for path, fn in iter_functions():
        body = fn.body
        stmts = [s for s in body if not isinstance(s, ast.Expr) or not isinstance(getattr(s, "value", None), ast.Constant)]
        if not stmts:
            bad.append(f"{path}:{fn.name} (docstring/pass only)")
        for s in stmts:
            if isinstance(s, ast.Pass):
                bad.append(f"{path}:{fn.name} (pass)")
            if isinstance(s, ast.Raise) and isinstance(s.exc, ast.Call) and getattr(s.exc.func, "id", "") == "NotImplementedError":
                bad.append(f"{path}:{fn.name} (NotImplementedError)")
    assert not bad, "Stub functions found:\n" + "\n".join(bad)
