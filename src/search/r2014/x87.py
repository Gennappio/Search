"""The 80-bit extended precision of x86, exactly.

On x86 R is built with ``long double`` accumulators of 80 bits (a 64-bit significand): sums of
squares in ``cor``, the mean in ``rowMeans``, the digits of a number being read, the scaling in
``round`` all carry eleven bits more than a double and are rounded to a double once, at the
end. On arm64 ``long double`` is a double and the same C code rounds at every step. The two
give different last digits, and the package then compares them with a threshold or with each
other, so a run on an Intel machine and one on an Apple Silicon machine can disagree.

To reproduce what R computes on x86 without x86, this module does the arithmetic on integers:
a value is ``m * 2**e`` with ``|m| < 2**64``, every operation is carried out exactly and
rounded once to 64 bits, to nearest, ties to even, which is what the x87 unit does in its
default mode. It is slow (a few operations per microsecond) and used only where R uses
``long double``.

Only finite values are represented; the callers deal with missing and infinite values before
they get here.
"""

from __future__ import annotations

import math

Ext = tuple[int, int]
"""A finite extended-precision value: (m, e) stands for m * 2**e."""

ZERO: Ext = (0, 0)
_BITS = 64


def _round(m: int, e: int, bits: int = _BITS) -> Ext:
    """m * 2**e rounded to ``bits`` significant bits, to nearest, ties to even."""
    if m == 0:
        return ZERO
    sign = -1 if m < 0 else 1
    a = -m if m < 0 else m
    shift = a.bit_length() - bits
    if shift <= 0:
        return (m, e)
    q = a >> shift
    rem = a & ((1 << shift) - 1)
    half = 1 << (shift - 1)
    if rem > half or (rem == half and q & 1):
        q += 1
        if q.bit_length() > bits:
            q >>= 1
            shift += 1
    return (sign * q, e + shift)


def from_double(x: float) -> Ext:
    """A double is an extended value without rounding."""
    if x == 0.0:
        return ZERO
    mantissa, exponent = math.frexp(x)
    return (int(mantissa * 9007199254740992.0), exponent - 53)


def from_int(n: int) -> Ext:
    return _round(n, 0)


def to_double(v: Ext) -> float:
    """``(double) v``: the second rounding, to 53 bits."""
    m, e = v
    if m == 0:
        return 0.0
    m53, e53 = _round(m, e, 53)
    return math.ldexp(float(m53), e53)


def add(a: Ext, b: Ext) -> Ext:
    (m1, e1), (m2, e2) = a, b
    if m1 == 0:
        return b
    if m2 == 0:
        return a
    e = e1 if e1 < e2 else e2
    return _round((m1 << (e1 - e)) + (m2 << (e2 - e)), e)


def neg(a: Ext) -> Ext:
    return (-a[0], a[1])


def sub(a: Ext, b: Ext) -> Ext:
    return add(a, neg(b))


def mul(a: Ext, b: Ext) -> Ext:
    return _round(a[0] * b[0], a[1] + b[1])


def div(a: Ext, b: Ext) -> Ext:
    """a / b for b != 0: the quotient to more bits than are kept, with a sticky bit."""
    (m1, e1), (m2, e2) = a, b
    if m2 == 0:
        raise ZeroDivisionError("division by zero in extended precision")
    if m1 == 0:
        return ZERO
    sign = -1 if (m1 < 0) != (m2 < 0) else 1
    a1, a2 = abs(m1), abs(m2)
    k = max(0, _BITS + 4 + a2.bit_length() - a1.bit_length())
    q, r = divmod(a1 << k, a2)
    e = e1 - e2 - k
    if r:
        q = (q << 1) | 1
        e -= 1
    return _round(sign * q, e)


def sqrt(a: Ext) -> Ext:
    """``sqrtl``: correctly rounded."""
    m, e = a
    if m < 0:
        raise ValueError("square root of a negative number")
    if m == 0:
        return ZERO
    if e & 1:
        m <<= 1
        e -= 1
    k = max(0, (2 * (_BITS + 4) - m.bit_length() + 1) // 2)
    m <<= 2 * k
    s = math.isqrt(m)
    e_out = (e - 2 * k) // 2
    if s * s != m:
        s = (s << 1) | 1
        e_out -= 1
    return _round(s, e_out)


def less(a: Ext, b: Ext) -> bool:
    return sub(a, b)[0] < 0


def nearbyint(a: Ext) -> Ext:
    """``nearbyintl``: to the nearest integer, ties to even."""
    m, e = a
    if e >= 0 or m == 0:
        return a
    shift = -e
    sign = -1 if m < 0 else 1
    am = abs(m)
    q = am >> shift
    rem = am & ((1 << shift) - 1)
    half = 1 << (shift - 1)
    if rem > half or (rem == half and q & 1):
        q += 1
    return _round(sign * q, 0)


def power_of_ten(n: int) -> Ext:
    """10**n for n >= 0 as R builds it, by repeated squaring in extended precision
    (``for (...; n; n >>= 1, p10 *= p10) if (n & 1) fac *= p10``). Exact up to 10**27."""
    fac = from_int(1)
    p10 = from_int(10)
    while n:
        if n & 1:
            fac = mul(fac, p10)
        n >>= 1
        if n:
            p10 = mul(p10, p10)
    return fac
