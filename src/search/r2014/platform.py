"""Which R the port stands for.

The same R code gives different last digits in two installations of R, for two reasons, and the
package turns last digits into decisions (a correlation against a threshold, a mean against a
null value, the fifth decimal of D).

- **The processor.** On x86 R's sums run in 80-bit extended precision and are rounded to a
  double at the end; on arm64 they are rounded at every step.
- **The version.** R 3.6.0 changed how ``sample()`` draws, R 4.0.0 how ``round()`` rounds, and
  smaller things changed with them (``cor`` did not use to clamp a correlation below -1).

A platform is a named pair of the two. Only those the port has been compared with, case by
case, are defined:

- :data:`R_4_ARM64`: R 4.4.3 on macOS arm64, the default.
- :data:`R_3_1_X86`: R 3.1.0 on Linux x86-64, the R released three days after the package was
  dated (2014-04-10), run from the image ``rocker/r-ver:3.1.0``.

The platform is chosen for a block of code with :func:`use`::

    with platform.use(platform.R_3_1_X86):
        SEARCH_randomise(...)
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class RPlatform:
    name: str
    r_version: tuple[int, int, int]
    extended: bool
    """``long double`` is the 80-bit type of x86 and R's accumulators use it."""

    @property
    def sample_kind(self) -> Literal["Rejection", "Rounding"]:
        """How ``sample()`` draws an index: by rejection from R 3.6.0, by rounding before."""
        return "Rejection" if self.r_version >= (3, 6, 0) else "Rounding"

    @property
    def round_to_representable(self) -> bool:
        """``round()`` picks the nearer of the two neighbouring decimals (R 4.0.0 and later)
        instead of computing ``rint(x * 10^d) / 10^d``."""
        return self.r_version >= (4, 0, 0)

    @property
    def clamps_at_minus_one(self) -> bool:
        """``cor`` brings a result below -1 back to -1. R 3.1.0 did so only above 1; R 4.4.3
        does both. The version in between that changed it was not looked up."""
        return self.r_version >= (4, 0, 0)


R_4_ARM64 = RPlatform("R 4.4.3, macOS arm64", (4, 4, 3), extended=False)
R_3_1_X86 = RPlatform("R 3.1.0, Linux x86-64", (3, 1, 0), extended=True)

_current: ContextVar[RPlatform] = ContextVar("search_r2014_platform", default=R_4_ARM64)


def current() -> RPlatform:
    return _current.get()


@contextmanager
def use(platform: RPlatform) -> Iterator[RPlatform]:
    """Run a block as the given R would compute it."""
    token = _current.set(platform)
    try:
        yield platform
    finally:
        _current.reset(token)


def matching(r_version: str, long_double_bytes: int) -> RPlatform | None:
    """The platform for an R that describes itself with ``R.version.string`` and
    ``.Machine$sizeof.longdouble``; ``None`` when the port was never compared with such an R."""
    extended = long_double_bytes > 8
    for candidate in (R_4_ARM64, R_3_1_X86):
        major, minor, _ = candidate.r_version
        if f"R version {major}.{minor}." in r_version and candidate.extended == extended:
            return candidate
    return None
