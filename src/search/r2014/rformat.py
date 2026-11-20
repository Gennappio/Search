"""How R prints a number in ``write.table``, and how it reads one back.

``write.table`` formats every number on its own with up to 15 significant digits
(``formatReal`` and ``scientific`` in ``src/main/format.c``, called with ``R_print.digits =
DBL_DIG``): the fewest digits that give the same 15-digit value, in fixed notation unless the
scientific one is shorter. So 0.0004 is written ``4e-04`` and 0.00047 is written ``0.00047``.
The files of the package are read back by ``summarise.iterations``, and are what a user of the
package has on disk, so the port writes them the same way, character for character.
"""

from __future__ import annotations

import math

import numpy as np

from search.r2014 import platform, x87
from search.r2014.rsem import NA_REAL, UnportedState, is_na_real

_DIGITS = 15
_KP_MAX = 22
_R_DEC_MIN_EXPONENT = -308


def _table_power(k: int) -> x87.Ext:
    """10^k as the table of ``format.c`` holds it in an R built with 80-bit ``long double``.
    The table is said to hold the powers a 64-bit significand represents exactly, up to
    10^27, but its entries are written as double constants: from 10^23 on they are the nearest
    double, not the power. Beyond the table R calls ``powl``; here the power is built by
    multiplication."""
    if k <= 27:
        return x87.from_double(float(10**k))
    return x87.power_of_ten(k)


def _scientific(r: float, digits: int) -> tuple[int, int, bool]:
    """``scientific()`` for r > 0: the exponent of 10, the number of significant digits (at
    most ``digits``), and whether rounding to them widens the fixed notation (9996 with three
    digits is 1e+04).

    R does not ask ``sprintf`` for the digits: it scales the number to an integer of
    ``digits`` digits with a power of ten taken from a table, rounds it, and counts the
    trailing zeros. The scaling is one multiplication or division in the precision of the
    platform, so the 15th digit of a number that needs all of them can differ from the
    correctly rounded one, and between a platform with extended precision and one without.

    In extended precision R takes the power of ten from a table up to 10^27 and from ``powl``
    of the C library beyond; here it is always built by multiplication, so a number smaller
    than 1e-13 or larger than 1e41 may be printed with a different last digit than R prints.
    """
    kp = int(math.floor(math.log10(r))) - digits + 1
    if platform.current().extended:
        scaled = x87.from_double(r)
        if kp > 0:
            scaled = x87.div(scaled, _table_power(kp))
        elif kp < 0:
            scaled = x87.mul(scaled, _table_power(-kp))
        if x87.less(scaled, x87.power_of_ten(digits - 1)):
            scaled = x87.mul(scaled, x87.from_int(10))
            kp -= 1
        alpha = x87.to_double(x87.nearbyint(scaled))
        kp_max = 27
    else:
        r_prec = r
        if abs(kp) <= _KP_MAX:
            if kp > 0:
                r_prec /= 10.0**kp
            elif kp < 0:
                r_prec *= 10.0 ** (-kp)
        elif kp <= _R_DEC_MIN_EXPONENT:
            r_prec = (r * 1e303) / math.pow(10.0, kp + 303)
        else:
            r_prec /= math.pow(10.0, kp)
        if r_prec < 10.0 ** (digits - 1):
            r_prec *= 10.0
            kp -= 1
        alpha = float(round(r_prec))  # nearbyint: to nearest, ties to even
        kp_max = _KP_MAX
    nsig = digits
    for _ in range(digits):
        alpha /= 10.0
        if alpha == math.floor(alpha):
            nsig -= 1
        else:
            break
    if nsig == 0 and digits > 0:
        nsig = 1
        kp += 1
    kpower = kp + digits - 1
    # scientific notation may round further than the fixed one
    rgt = min(max(digits - kpower, 0), kp_max)
    fuzz = 0.5 / 10.0**rgt
    widens = False
    if 0 < kpower <= kp_max:
        if platform.current().extended:
            limit = x87.sub(_table_power(kpower), x87.from_double(fuzz))
            widens = x87.less(x87.from_double(r), limit)
        else:
            widens = r < 10.0**kpower - fuzz
    return kpower, nsig, widens


def format_real(x: float, digits: int = _DIGITS) -> str:
    """One number as ``write.table`` writes it."""
    if math.isnan(x):
        return "NA" if bool(is_na_real(np.array([x]))[0]) else "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0.0:
        return "0"
    neg = 1 if x < 0 else 0
    kpower, nsig, widens = _scientific(abs(x), digits)
    left = kpower + 1
    if widens:
        left -= 1
    sleft = neg + (1 if left <= 0 else left)
    rgt = max(nsig - left, 0)
    wide_exponent = left > 100 or left <= -99
    width_fixed = sleft + rgt + (1 if rgt != 0 else 0)
    d = nsig - 1
    width_sci = neg + (1 if d > 0 else 0) + d + 4 + (2 if wide_exponent else 1)
    if width_fixed <= width_sci:
        return f"{x:.{rgt}f}"
    text = f"{x:.{d}e}"
    mantissa, exponent = text.split("e")
    sign = exponent[0]
    power = exponent[1:].lstrip("0") or "0"
    return f"{mantissa}e{sign}{power.zfill(3 if wide_exponent else 2)}"


def r_strtod(text: str) -> float:
    """``R_strtod``: how R turns the text of a number into a double when it reads a file.

    R does not call the C library. It adds up the digits one by one, ``ans = 10 * ans +
    digit``, and scales the result by a power of ten at the end, all in a ``long double``.

    Where that is a double (arm64), every step is exact up to 15 significant digits and the
    result is the nearest double, the same ``float`` gives; with more digits the running sum
    is rounded at every step, and one number in 50 is read into a neighbouring double with 16
    digits, one in four with 17. Where it is the 80-bit type (x86) the sum holds 19 digits
    exactly, but the scaled result is rounded twice, to 64 bits and then to 53, and now and
    then a number of few digits comes out one unit in the last place away.
    """
    p = text.strip()
    sign = 1.0
    if p[:1] in ("+", "-"):
        sign = -1.0 if p[0] == "-" else 1.0
        p = p[1:]
    lowered = p.lower()
    if lowered == "nan":
        return float("nan")
    if lowered in ("inf", "infinity"):
        return sign * float("inf")
    if lowered.startswith("0x"):
        return sign * float.fromhex(p)
    digits: list[int] = []
    expn = 0
    i = 0
    while i < len(p) and p[i].isdigit():
        digits.append(ord(p[i]) - 48)
        i += 1
    if i < len(p) and p[i] == ".":
        i += 1
        while i < len(p) and p[i].isdigit():
            digits.append(ord(p[i]) - 48)
            i += 1
            expn -= 1
    if not digits:
        raise ValueError(f"not a number: {text!r}")
    if i < len(p) and p[i] in "eE":
        i += 1
        exponent_sign = 1
        if i < len(p) and p[i] in "+-":
            exponent_sign = -1 if p[i] == "-" else 1
            i += 1
        n = 0
        while i < len(p) and p[i].isdigit():
            n = n * 10 + (ord(p[i]) - 48)
            i += 1
        expn += exponent_sign * n
    if i != len(p):
        raise ValueError(f"not a number: {text!r}")
    if platform.current().extended:
        return sign * _scale_extended(digits, expn)
    ans = 0.0
    for digit in digits:
        ans = 10.0 * ans + digit
    # avoid unnecessary underflow for large negative exponents
    if expn + len(digits) < -300:
        for _ in digits:
            ans /= 10.0
        expn += len(digits)
    fac = 1.0
    p10 = 10.0
    n = abs(expn)
    if expn < -307:  # use underflow, not overflow
        while n:
            if n & 1:
                fac /= p10
            n >>= 1
            p10 *= p10
        ans *= fac
    elif expn < 0:  # positive powers are exact
        while n:
            if n & 1:
                fac *= p10
            n >>= 1
            p10 *= p10
        ans /= fac
    elif ans != 0.0:
        while n:
            if n & 1:
                fac *= p10
            n >>= 1
            p10 *= p10
        ans *= fac
    return sign * ans


def _scale_extended(digits: list[int], expn: int) -> float:
    """The digits and the exponent of a number, put together in 80-bit precision as
    ``R_strtod5`` of R 3.1.0 does."""
    if expn + len(digits) < -300 or expn < -307 or expn > 300:
        raise UnportedState("a number with an exponent near the limits of a double")
    ten = x87.from_int(10)
    ans = x87.ZERO
    for digit in digits:
        ans = x87.add(x87.mul(ten, ans), x87.from_int(digit))
    if expn < 0:
        ans = x87.div(ans, x87.power_of_ten(-expn))
    elif expn > 0:
        ans = x87.mul(ans, x87.power_of_ten(expn))
    return x87.to_double(ans)


def parse_real(token: str) -> float:
    """A field of a file as ``read.table`` reads it into a numeric column."""
    t = token.strip()
    if t == "NA":
        return NA_REAL
    return r_strtod(t)
