"""The random number generator of R, as far as SEARCH 0.0.2 uses it.

``.SEARCH.randomise.inner`` calls ``set.seed(iteration + sampling.seed - 1)`` and then
``sample(idx, replace = TRUE)`` (bootstrap) or ``sample(idx, size = n - 1, replace = FALSE)``
(jackknife). To draw the same samples as R from the same seed, this module reproduces:

- the default generator of R, the Mersenne Twister, seeded as ``set.seed`` seeds it (a linear
  congruential scrambling of the seed, ``RNG_Init`` in ``src/main/RNG.c``);
- ``unif_rand``: a 32-bit output scaled to (0, 1);
- ``sample``: ``do_sample`` and ``R_unif_index``, in the two variants R has had. Up to R 3.5 an
  index was ``floor(n * unif_rand())`` (``sample.kind = "Rounding"``); from R 3.6.0 it is drawn
  by rejection from the integers below the next power of two (``"Rejection"``, the default).
  The package was written for R 2.15 to 3.1, so its original runs used ``"Rounding"``; an R of
  today gives those draws only after ``RNGkind(sample.kind = "Rounding")``.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
from numpy.typing import NDArray

SampleKind = Literal["Rejection", "Rounding"]

_MASK32 = 0xFFFFFFFF
_I2_32M1 = 2.328306437080797e-10  # 1 / (2^32 - 1), the bound of fixup() in RNG.c


class RRandom:
    """``set.seed(seed)`` followed by calls of ``unif_rand`` and ``sample``."""

    def __init__(self, seed: int, sample_kind: SampleKind = "Rejection") -> None:
        if sample_kind not in ("Rejection", "Rounding"):
            raise ValueError(f"unknown sample kind {sample_kind!r}")
        self.sample_kind: SampleKind = sample_kind
        # RNG_Init: initial scrambling, then one linear congruential step per state word
        s = int(seed) & _MASK32
        for _ in range(50):
            s = (69069 * s + 1) & _MASK32
        words = np.empty(625, dtype=np.uint32)
        for j in range(625):
            s = (69069 * s + 1) & _MASK32
            words[j] = s
        # the first word is the position in the state (dummy[0] = mti), set to 624 by FixupSeeds
        self._bits = np.random.MT19937()
        self._bits.state = {
            "bit_generator": "MT19937",
            "state": {"key": words[1:].copy(), "pos": 624},
        }

    def unif_rand(self) -> float:
        value = float(int(self._bits.random_raw())) * 2.3283064365386963e-10
        # fixup(): keep the value inside the open interval
        if value <= 0.0:
            return 0.5 * _I2_32M1
        if 1.0 - value <= 0.0:
            return 1.0 - 0.5 * _I2_32M1
        return value

    def _rbits(self, bits: int) -> float:
        v = 0
        n = 0
        while n <= bits:
            v1 = int(math.floor(self.unif_rand() * 65536))
            v = 65536 * v + v1
            n += 16
        if bits < 64:
            v &= (1 << bits) - 1
        return float(v)

    def unif_index(self, dn: float) -> float:
        """``R_unif_index``: a random integer in 0 .. dn - 1."""
        if self.sample_kind == "Rounding":
            return float(math.floor(dn * self.unif_rand()))
        if dn <= 0:
            return 0.0
        bits = int(math.ceil(math.log2(dn)))
        while True:
            dv = self._rbits(bits)
            if dn > dv:
                return dv

    def sample(self, n: int, size: int | None = None, replace: bool = False) -> NDArray[np.int64]:
        """``sample(1:n, size, replace)``: 1-based indices, as R returns them."""
        k = n if size is None else size
        if not replace and k > n:
            raise ValueError(
                "cannot take a sample larger than the population when 'replace = FALSE'"
            )
        out = np.empty(k, dtype=np.int64)
        if replace or k < 2:
            for i in range(k):
                out[i] = int(self.unif_index(float(n)) + 1)
            return out
        x = list(range(n))
        m = n
        for i in range(k):
            j = int(self.unif_index(float(m)))
            out[i] = x[j] + 1
            m -= 1
            x[j] = x[m]
        return out
