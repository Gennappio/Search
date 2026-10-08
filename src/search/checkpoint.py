"""Chunked, resumable storage of per-replicate statistics.

Each chunk is one parquet file (gene column plus one column per replicate) with a sidecar
``meta.json``. A chunk is reused on resume only when its metadata matches the current
configuration and inputs; anything else is recomputed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray


def chunk_ranges(n_rep: int, chunk_size: int) -> list[tuple[int, int]]:
    """1-based inclusive replicate ranges: [(1, 100), (101, 200), ...]."""
    if n_rep <= 0:
        return []
    size = max(1, chunk_size)
    return [(start, min(start + size - 1, n_rep)) for start in range(1, n_rep + 1, size)]


def rep_column(rep: int) -> str:
    return f"rep_{rep:05d}"


@dataclass
class ChunkMeta:
    stage: str
    rep_start: int
    rep_end: int
    seed_base: int
    config_hash: str
    input_hash: str
    extra: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage": self.stage,
            "rep_start": self.rep_start,
            "rep_end": self.rep_end,
            "seed_base": self.seed_base,
            "config_hash": self.config_hash,
            "input_hash": self.input_hash,
            "extra": self.extra,
        }

    def matches(self, other: dict[str, Any]) -> bool:
        mine = self.to_dict()
        return all(
            mine[key] == other.get(key)
            for key in ("stage", "rep_start", "rep_end", "seed_base", "config_hash", "input_hash")
        )


class ChunkStore:
    """Parquet chunks of one stage, for one cohort, under ``directory``."""

    def __init__(self, directory: Path | str, stage: str, dtype: str = "float64") -> None:
        self.directory = Path(directory)
        self.stage = stage
        self.dtype = np.float32 if dtype == "float32" else np.float64
        self.directory.mkdir(parents=True, exist_ok=True)

    def path(self, chunk_no: int) -> Path:
        return self.directory / f"{self.stage}_chunk_{chunk_no:04d}.parquet"

    def meta_path(self, chunk_no: int) -> Path:
        return self.directory / f"{self.stage}_chunk_{chunk_no:04d}.meta.json"

    def aux_path(self, chunk_no: int, name: str) -> Path:
        return self.directory / f"{self.stage}_{name}_{chunk_no:04d}.parquet"

    def is_valid(self, chunk_no: int, meta: ChunkMeta) -> bool:
        if not self.path(chunk_no).exists() or not self.meta_path(chunk_no).exists():
            return False
        try:
            stored = json.loads(self.meta_path(chunk_no).read_text())
        except (OSError, ValueError):
            return False
        return meta.matches(stored)

    def write(
        self,
        chunk_no: int,
        genes: pd.Index,
        columns: dict[int, NDArray[np.float64]],
        meta: ChunkMeta,
        aux: dict[str, pd.DataFrame] | None = None,
    ) -> None:
        frame = pd.DataFrame(
            {
                rep_column(rep): np.asarray(values, dtype=self.dtype)
                for rep, values in columns.items()
            }
        )
        frame.insert(0, "gene", np.asarray([str(g) for g in genes]))
        tmp = self.path(chunk_no).with_suffix(".tmp.parquet")
        frame.to_parquet(tmp, index=False)
        tmp.replace(self.path(chunk_no))
        for name, table in (aux or {}).items():
            table.to_parquet(self.aux_path(chunk_no, name), index=False)
        self.meta_path(chunk_no).write_text(json.dumps(meta.to_dict(), indent=1))

    def read(self, chunk_no: int) -> pd.DataFrame:
        frame = pd.read_parquet(self.path(chunk_no))
        return frame.set_index("gene")

    def read_aux(self, chunk_no: int, name: str) -> pd.DataFrame:
        return pd.read_parquet(self.aux_path(chunk_no, name))

    def read_all(self, n_chunks: int) -> pd.DataFrame:
        frames = [self.read(i) for i in range(1, n_chunks + 1)]
        return pd.concat(frames, axis=1)

    def delete(self, chunk_no: int) -> None:
        for path in (self.path(chunk_no), self.meta_path(chunk_no)):
            if path.exists():
                path.unlink()
