"""Static code-hygiene checks."""

from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = ROOT / "src" / "trajectory"
SCRIPTS_DIR = ROOT / "scripts"

ALL_PY = sorted(SRC_DIR.rglob("*.py")) + sorted(SCRIPTS_DIR.glob("*.py"))


def _read(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def test_no_sys_path_insert():
    """Scripts should import from the installed package, not hack sys.path."""
    pattern = re.compile(r"sys\.path\.(insert|append)\(")
    violations = []
    for p in ALL_PY:
        for i, line in enumerate(_read(p).splitlines(), 1):
            if pattern.search(line) and not line.lstrip().startswith("#"):
                violations.append(f"{p.relative_to(ROOT)}:{i}")
    assert violations == [], (
        f"sys.path manipulation found (use pip install -e . instead):\n"
        + "\n".join(f"  {v}" for v in violations)
    )


def test_no_ndarray_np_number():
    """NDArray[np.number] is invalid — use np.floating, np.integer, etc."""
    pattern = re.compile(r"NDArray\[np\.number\]")
    violations = []
    for p in ALL_PY:
        for i, line in enumerate(_read(p).splitlines(), 1):
            if pattern.search(line):
                violations.append(f"{p.relative_to(ROOT)}:{i}")
    assert violations == [], (
        f"Invalid NDArray[np.number] annotation:\n"
        + "\n".join(f"  {v}" for v in violations)
    )


def _find_unused_imports(filepath: Path) -> list[str]:
    """Return list of 'name (line N)' for imports whose names aren't used."""
    source = _read(filepath)
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []

    imports: list[tuple[str, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                name = alias.asname or alias.name.split(".")[0]
                imports.append((name, node.lineno))
        elif isinstance(node, ast.ImportFrom):
            if node.module == "__future__":
                continue
            for alias in node.names:
                name = alias.asname or alias.name
                imports.append((name, node.lineno))

    unused = []
    for name, lineno in imports:
        if name.startswith("_"):
            continue
        count = len(re.findall(rf"\b{re.escape(name)}\b", source))
        if count <= 1:
            unused.append(f"{name} (line {lineno})")
    return unused


def test_no_unused_imports():
    """Every imported name should be referenced at least once beyond the import."""
    violations = {}
    for p in ALL_PY:
        unused = _find_unused_imports(p)
        if unused:
            violations[str(p.relative_to(ROOT))] = unused
    assert violations == {}, (
        "Unused imports found:\n"
        + "\n".join(
            f"  {f}: {', '.join(names)}" for f, names in violations.items()
        )
    )
