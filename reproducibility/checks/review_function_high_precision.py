"""High-precision cross-check of Section 10, separate from moving-window sums.

Uses 110-digit mpmath, a separate digamma pair for every partial-fraction
pole, cutoff 128 and order 56. The published table uses 60-digit Arb,
cutoff 192 and order 40. Exact threshold and expansion generators are
shared; this is a second numerical route, not a separate proof of those
identities. The reported tail bound excludes floating-point roundoff.
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from flint import fmpq, fmpq_poly

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "reproducibility" / "alper_function"))
import alper_function as af

CUTOFF, ORDER, DIGITS = 128, 56, 110


def mpq(q):
    return mp.mpf(int(q.numerator)) / int(q.denominator)


def evaluate(h, moments):
    h0 = int(h >= 1)
    phi = h - h0
    assert CUTOFF * (1 - phi) >= 1
    thresholds = af.thresholds(h, CUTOFF)
    c = [mp.mpf(1)]
    for x in range(1, CUTOFF):
        c.append(c[-1] * (2*x - 1) / (2*x))
    half = mpq(h) / 2
    head = mp.mpf(0)
    for k in range(1, CUTOFF):
        first, last = thresholds[k], thresholds[k+1] - 2
        # Evaluate every pole separately: do not update a digamma window.
        terms = []
        for x in range(k):
            weight_over_two = c[x] * (2*k - 2*x - 1) * c[k-1-x]
            delta = mp.digamma(last - half - 2*x) - mp.digamma(first - half - 1 - 2*x)
            terms.append(weight_over_two * delta)
        B = (last-first+1 - mp.mpf(2*k*k+1) * mp.fsum(terms) / 3) / (2*k+1)
        head += B - mp.mpf(3)/2 + mp.mpf(3)/(4*k)
    k = fmpq_poly([0, 1])
    tail, remainders = mp.mpf(0), []
    for parity in (0, 1):
        eta = (-1)**(parity+h0)
        length = 2*k**2 + 2*k + fmpq(1+eta, 2)
        offset = fmpq(5, 3)*k + fmpq(1-eta, 2) - phi
        coefficients, remainder = af.expansion(ORDER, length, offset, moments)
        remainders.append(remainder)
        first = CUTOFF + (parity-CUTOFF) % 2
        for j in range(2, ORDER):
            tail += mpq(coefficients[j]) * mp.zeta(j, mp.mpf(first)/2) / 2**j
    geometric, polynomial = af.ac.truncation_bound(ORDER, CUTOFF, remainders)
    value = mp.mpf(3)/2 + 3*mp.euler/4 + mp.log(mp.mpf(3)/2)/4 - head - tail
    return value, mpq(geometric + polynomial)


def main():
    start = time.perf_counter()
    mp.mp.dps = DIGITS
    table = json.loads((ROOT / "reproducibility/alper_function/alper_function_values.json").read_text())["table"]
    moments = af.ac.moment_polynomials(ORDER)
    results = {}
    for key in ("0.3", "0.8", "1.3", "1.9"):
        fraction = Fraction(key)
        h = fmpq(fraction.numerator, fraction.denominator)
        value, tail_bound = evaluate(h, moments)
        saved = table[key]["A_ball"]
        midpoint = mp.mpf(saved.lstrip("[").split()[0])
        radius = mp.mpf(saved.split("+/-")[1].strip(" ]"))
        difference = abs(value-midpoint)
        assert tail_bound < mp.mpf("1e-80")
        assert difference < radius + tail_bound + mp.mpf("1e-95")
        results[key] = {"value": mp.nstr(value, 95),
                        "difference_from_saved_midpoint": mp.nstr(difference, 12),
                        "tail_bound": mp.nstr(tail_bound, 12)}
        print(key, results[key], flush=True)
    report = {"status": "PASS", "mpmath_digits": DIGITS, "cutoff": CUTOFF,
              "order": ORDER, "values": results,
              "shared_components": "exact threshold, moment and expansion generators",
              "seconds": round(time.perf_counter()-start, 2)}
    Path(__file__).with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n")
    print("PASS", report["seconds"], "seconds", flush=True)


if __name__ == "__main__":
    main()
