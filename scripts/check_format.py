"""Fail when Mojo sources would be changed by `mojo format`."""

from __future__ import annotations

import filecmp
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SOURCES = (
    Path("src/graphical_counts/__init__.mojo"),
    Path("tests/test_graphical_counts.mojo"),
    Path("examples/basic.mojo"),
)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mojo-format-check-") as temp:
        root = Path(temp)
        copies: list[Path] = []
        for source in SOURCES:
            target = root / source
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            copies.append(target)

        result = subprocess.run(["mojo", "format", "--quiet", *map(str, copies)])
        if result.returncode != 0:
            return result.returncode

        changed = [
            source
            for source, formatted in zip(SOURCES, copies, strict=True)
            if not filecmp.cmp(source, formatted, shallow=False)
        ]
        if changed:
            print("Mojo formatting required:", file=sys.stderr)
            for source in changed:
                print(f"  {source}", file=sys.stderr)
            print("Run `pixi run format`.", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
