"""The run manifest: everything needed to know what produced a result."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np

from search import __version__
from search.data import ExpressionSet


def versions() -> dict[str, str]:
    out = {"search": __version__, "python": sys.version.split()[0], "platform": platform.platform()}
    for name in ("numpy", "scipy", "pandas", "pydantic", "pyarrow", "joblib"):
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            out[name] = "absent"
    return out


def file_sha256(path: Path | str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def expression_hash(expr: ExpressionSet, seeds: list[str], universe: list[str] | None) -> str:
    """Content hash of the matrix, its gene and sample names, the seeds and the universe."""
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(expr.values).tobytes())
    h.update("\n".join(str(g) for g in expr.genes).encode())
    h.update("\n".join(str(s) for s in expr.samples).encode())
    h.update("\n".join(seeds).encode())
    if universe is not None:
        h.update("\n".join(universe).encode())
    return h.hexdigest()[:16]


def now_iso() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def write_manifest(run_dir: Path | str, manifest: dict[str, Any]) -> Path:
    path = Path(run_dir) / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, default=str))
    return path


def read_manifest(run_dir: Path | str) -> dict[str, Any]:
    return json.loads((Path(run_dir) / "manifest.json").read_text())
