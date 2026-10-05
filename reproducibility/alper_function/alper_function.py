"""Certified values of the Alper function A(h), 0 <= h < 2, h outside the exceptional set E.

A(h) = 3/2 + (3/4) gamma + (1/4) log(3/2) - sum_k (B_k(h) - 3/2 + 3/(4k)),
B_k(h) = sum of P_k(n - h/2) over the block k of G_n(h) (a win with m pays 2m - h).

Write h0 = floor(h) and phi = h - h0.  The thresholds t_k(h) are found by exact sign tests for
k <= CUTOFF + 1.  Beyond, t_k(h) = t_k(h0) because theta_k < 1/(2k) (Lemma 10.1 of the paper) and
CUTOFF >= 1/(1 - phi) is asserted.  The tail sum uses the exact block expansion of the certificate
program with the offset a_eta replaced by a_eta - phi, Hurwitz zeta values and the same rigorous
truncation bounds.  All arithmetic is Arb ball arithmetic (python-flint).

Run:  python alper_function.py      (writes alper_function_values.json and table_alper_rows.tex)
"""
import json
import sys
import time
import platform
from math import comb
from pathlib import Path

import flint
from flint import arb, fmpq, fmpq_poly, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "certificate"))
import alper_constant_certificate as ac  # noqa: E402

CUTOFF, ORDER, DIGITS = 192, 40, 60
ctx.dps = DIGITS


def thresholds(h, kmax):
    """t_k(h) for 1 <= k <= kmax + 1: the first n with P_k(n - h/2) > 0 (exact rational sign test)."""
    t = {1: 3}
    for k in range(2, kmax + 2):
        n = ac.threshold(k, 0)          # t_k(h) >= t_k(0), since floor(beta_k + h/2) grows with h
        while True:
            y = fmpq(n) - h / 2
            num, den = fmpq(1), fmpq(1)
            for j in range(1, k + 1):
                num *= y - 2 * j
                den *= y - 2 * j + 1
            value = (2 * k * k + 1) * num - 2 * (k * k - 1) * den   # sign of P_k(y), since den > 0
            if value == 0:
                raise ValueError("h is exceptional (k = %d)" % k)
            if value > 0:
                break
            n += 1
        t[k] = n
    return t


def block_ball(k, c, h, first, last):
    """B_k from the digamma formula with a moving window, for real h."""
    half = arb(h) / 2
    harmonic = ((arb(last) - half).digamma() - (arb(first - 1) - half).digamma()) / 2
    weighted = arb(0)
    for x in range(k):
        weight = 2 * c[x] * (2 * k - 2 * x - 1) * c[k - 1 - x]
        weighted += arb(weight) * harmonic
        if x + 1 < k:
            harmonic += (1 / (2 * arb(first) - arb(h) - 4 * x - 6) + 1 / (2 * arb(first) - arb(h) - 4 * x - 4)
                         - 1 / (2 * arb(last) - arb(h) - 4 * x - 4) - 1 / (2 * arb(last) - arb(h) - 4 * x - 2))
    return (last - first + 1 - arb(fmpq(2 * k * k + 1, 3)) * weighted) / (2 * k + 1)


def expansion(order, length, offset, moments):
    """Coefficients b_j and the remainder polynomial for given ell_k = length(k), a = offset(k)."""
    mmax, degree = order, order - 1
    bernoulli = ac.bernoulli_numbers(mmax)
    length_powers = [length**j for j in range(mmax + 2)]
    offset_powers = [offset**j for j in range(mmax + 1)]
    powersums, centered = [], []
    for j in range(mmax + 1):
        powersums.append(sum((comb(j + 1, l) * bernoulli[l] * length_powers[j + 1 - l] for l in range(j + 1)),
                             fmpq_poly([])) / (j + 1))
        centered.append(sum((comb(j, l) * (-4)**l * moments[l] * offset_powers[j - l] for l in range(j + 1)),
                            fmpq_poly([])))
    f = [fmpq(0) for _ in range(3 * mmax + 4)]
    for d in range(length.degree() + 1):
        f[3 - d] -= length[d] / 2
    for m in range(1, mmax + 1):
        t = sum((comb(m, j) * 2**j * powersums[j] * centered[m - j] for j in range(m + 1)), fmpq_poly([]))
        factor = -fmpq(-3, 4)**m
        for d in range(t.degree() + 1):
            p = 3 * m + 1 - d
            assert p >= 0
            f[p] += factor * t[d]
            f[p + 2] += factor * t[d] / 2
    coefficients = []
    for j in range(degree + 1):
        coefficients.append((f[j] - (coefficients[-1] if j else 0)) / 2)
    g = f.copy()
    for j, coefficient in enumerate(coefficients):
        g[j] -= 2 * coefficient
        g[j + 1] -= coefficient
    assert all(v == 0 for v in g[:degree + 1])
    return coefficients, g


MOMENTS = ac.moment_polynomials(ORDER)
C = ac.central_coefficients(CUTOFF + 2)


def regular_expansions(h):
    """Expansions for the regular blocks: thresholds of h0, offset a_eta - phi; index 0 = even k."""
    h0 = 0 if h < 1 else 1
    phi = h - h0
    k = fmpq_poly([0, 1])
    coeff, remain = [], []
    for sign in (1, -1):
        eta = sign if h0 == 0 else -sign            # eta = (-1)^(k + h0)
        length = 2 * k**2 + 2 * k + fmpq(1 + eta, 2)
        offset = fmpq(5, 3) * k + fmpq(1 - eta, 2) - phi
        a, g = expansion(ORDER, length, offset, MOMENTS)
        coeff.append(a)
        remain.append(g)
    return coeff, remain


def alper(h):
    """Ball enclosure of A(h), the thresholds t_1..t_{CUTOFF+1}, and the expansion coefficients."""
    h = fmpq(h)
    assert 0 <= h < 2
    h0 = 0 if h < 1 else 1
    phi = h - h0
    assert CUTOFF * (1 - phi) >= 1, "Lemma 10.1 does not cover the tail; raise CUTOFF"
    t = thresholds(h, CUTOFF)
    for k in (CUTOFF, CUTOFF + 1):
        assert t[k] == ac.threshold(k, h0), (h, k)
    finite = arb(0)
    for k in range(1, CUTOFF):
        finite += block_ball(k, C, h, t[k], t[k + 1] - 2) - arb(fmpq(3, 2)) + arb(fmpq(3, 4 * k))
    coeff, remain = regular_expansions(h)
    tail = arb(0)
    for parity in (0, 1):
        first = CUTOFF + (parity - CUTOFF) % 2
        for j in range(2, ORDER):
            tail += arb(coeff[parity][j]) * arb(j).zeta(arb(fmpq(first, 2))) / 2**j
    geometric, polynomial = ac.truncation_bound(ORDER, CUTOFF, remain)
    total = finite + tail + arb(0, arb(geometric + polynomial).upper())
    A = arb(fmpq(3, 2)) - 3 * arb(1).digamma() / 4 + arb(fmpq(3, 2)).log() / 4 - total
    return A, t, coeff


def defect(h, t):
    """The k with t_k(h) != t_k(h0); Lemma 10.1: t_k(h) = t_k(h0) + 1 there."""
    h0 = 0 if h < 1 else 1
    out = [k for k in range(2, CUTOFF + 2) if t[k] != ac.threshold(k, h0)]
    assert all(t[k] == ac.threshold(k, h0) + 1 and k % 2 == h0 for k in out)
    return out


def main():
    start = time.perf_counter()
    result = {"python_version": platform.python_version(), "python_flint_version": flint.__version__,
              "cutoff": CUTOFF, "order": ORDER, "digits": DIGITS}
    # 1. the first coefficients of the regular-block expansion: 3/2, -3/4, 3 eta/8 - 3 phi/4
    for h in [fmpq(j, 20) for j in range(40)]:
        h0 = 0 if h < 1 else 1
        coeff, _ = regular_expansions(h)
        for parity in (0, 1):
            eta = (-1)**(parity + h0)
            assert coeff[parity][:3] == [fmpq(3, 2), fmpq(-3, 4), fmpq(3 * eta, 8) - 3 * (h - h0) / 4], (h, parity)
    print("Expansion check: B_k(h) = 3/2 - 3/(4k) + (3 eta/8 - 3 phi/4)/k^2 + ... for h = j/20.", flush=True)
    # 2. consistency with the 80-digit certificates at h = 0 and h = 1
    for h in (0, 1):
        cert = json.loads((HERE.parent / "certificate" / ("certificate_h%d.json" % h)).read_text(encoding="utf-8"))
        A, _, _ = alper(h)
        assert A.overlaps(arb(cert["certified_80_decimals"]["A"])), h
    print("A(0) and A(1) agree with the 80-digit certificates.", flush=True)
    # 3. the table h = j/10
    table = {}
    for j in range(20):
        h = fmpq(j, 10)
        A, t, _ = alper(h)
        # Preserve the certified enclosure in the file, not just its first
        # 30 displayed digits. Check the decimal -> Arb roundtrip as well.
        ball_text = A.str(DIGITS + 8, more=True)
        saved = arb(ball_text)
        assert saved.contains(A), "serialized ball must enclose the computed ball"
        assert saved.rad() < arb("3e-54"), "serialized radius must support the paper's bound"
        table["%.1f" % (j / 10)] = {"A_ball": ball_text, "radius": saved.rad().str(20, more=True),
                                     "A_20_decimals": ac.rounded_decimal(saved, 20), "A_10_decimals": ac.rounded_decimal(saved, 10),
                                     "thresholds_2_to_8": [t[k] for k in range(2, 9)], "defect_k": defect(h, t)}
        print("h = %.1f  A(h) = %s   defect k = %s" % (j / 10, ac.rounded_decimal(A, 20), defect(h, t)), flush=True)
    result["table"] = table
    # 4. the first jump: h* = 6 - 2 sqrt(7), from beta_2 = 5 + sqrt(7)
    sqrt7 = arb(7).sqrt()
    hstar = 6 - 2 * sqrt7
    scale = 10**12
    lo = fmpq(int((hstar.mid() * scale).floor().unique_fmpz()), scale)
    hi = lo + fmpq(1, scale)
    assert hstar > arb(lo) and hstar < arb(hi)
    left, t_left, _ = alper(lo)            # within 1e-12 of h*: approximates the one-sided limits
    right, t_right, _ = alper(hi)
    assert t_left[2] == 8 and t_right[2] == 9 and all(t_left[k] == t_right[k] for k in range(3, CUTOFF + 2))
    predicted = (sqrt7 - 1) / 6
    result["jump_k2"] = {"h_star": hstar.str(20), "h_left": str(lo), "h_right": str(hi),
                         "A_left": left.str(25), "A_right": right.str(25),
                         "difference": (left - right).str(20), "predicted_(sqrt7-1)/6": predicted.str(20)}
    assert abs(float((left - right - predicted).mid())) < 1e-10
    print("A(h*-) ~ %s, A(h*+) ~ %s, jump %s, predicted %s" % (left.str(15), right.str(15), (left - right).str(15),
                                                                predicted.str(15)), flush=True)
    result["seconds"] = time.perf_counter() - start
    (HERE / "alper_function_values.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    # LaTeX rows: h and A(h) to 10 decimals, two columns of ten
    rows = []
    for j in range(10):
        a, b = table["%.1f" % (j / 10)], table["%.1f" % (1 + j / 10)]
        rows.append(r"$%.1f$ & %s & $%.1f$ & %s \\" % (j / 10, a["A_10_decimals"], 1 + j / 10, b["A_10_decimals"]))
    (HERE / "table_alper_rows.tex").write_text("\n".join(rows) + "\n", encoding="utf-8")
    print("\n".join(rows))
    print("seconds: %.1f" % result["seconds"])


if __name__ == "__main__":
    main()
