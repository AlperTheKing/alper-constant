"""Checks of the displayed identities, examples, tables and asymptotic claims of the paper
"Parity Rock--Paper--Scissors Games and the Alper Constant".

G_n(h): a win with the number m pays 2m - h.  Exact rational arithmetic (fractions.Fraction) is
used for every algebraic identity and equilibrium; mpmath (50 digits) for digamma block sums and
asymptotic statements; SymPy for the symbolic identities; SciPy, if installed, for independent
floating-point linear programming.  Finite checks support, but do not replace, the written proofs.

Run from this directory:  python verify_paper.py
"""
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from fractions import Fraction as F
from math import comb
from pathlib import Path
import functools
import json
import platform
import random
import re
import sys
import time

import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CERT = {h: json.loads((ROOT / "reproducibility" / "certificate" / ("certificate_h%d.json" % h)).read_text(encoding="utf-8"))
        for h in (0, 1)}
CPP = {0: ROOT / "reproducibility" / "figures" / "S_h0_40dp.txt"}
TEX = ROOT / "main.tex"

mp.mp.dps = 50
COUNT = {"assertions": 0}
REPORT = {}
KAPPA = mp.cbrt(mp.mpf(81) / 16)
H01 = (F(0), F(1))


def check(condition, message):
    COUNT["assertions"] += 1
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------- the games
def payoff_rule(i, j, h):
    """Payment received by i against j, straight from the verbal rule (payment 2w - h)."""
    if i == j:
        return F(0)
    w = min(i, j) if (i - j) % 2 else max(i, j)
    return (2 * w - h) if w == i else -(2 * w - h)


def sgn(x):
    return (x > 0) - (x < 0)


def a(i, j, h):
    """Lemma 2.1, valid for all integers."""
    return (i - j) + sgn(i - j) * (-1) ** ((i + j) % 2) * (i + j - h)


def t(k, h):
    """Theorem A: t_k(0) = ceil((4k^3+5k)/6)+1, t_k(1) = floor((4k^3+5k)/6)+2."""
    return (4 * k**3 + 5 * k) // 6 + 2 if h == 1 else (4 * k**3 + 5 * k + 5) // 6 + 1


def R(k, y):
    r = F(1) if isinstance(y, (int, F)) else mp.mpf(1)
    for j in range(1, k + 1):
        r *= (y - 2 * j) / (y - 2 * j + 1)
    return r


def P(k, y):
    return ((2 * k * k + 1) * R(k, y) - 2 * (k * k - 1)) / (3 * (2 * k + 1))


def cc(j):
    return F(comb(2 * j, j), 4**j)


def w(k, j):
    return 2 * (2 * k - 2 * j + 1) * cc(j - 1) * cc(k - j)


def block_of(n, h):
    k = 1
    while t(k + 1, h) <= n:
        k += 1
    return k


def payoffs(p, i, h):
    return sum((a(i, j, h) * pj for j, pj in p.items()), F(0))


def interval_solution(k, b, h):
    """Theorem 4.3: chain formulas (eq:chain1), (eq:chain2)."""
    a_ = b - 2 * k
    L = 2 * k + 1
    y = b - h / 2
    Ry = R(k, y)
    K = F(a_ + b - 1, 2) - F(1, 2 * L)
    c = (1 - h) / 3
    Q = (2 * y * (Ry - 1) + 2 * k) / 3 + F(1, L)
    q = {a_ - 2: F(0)}
    for m in range(a_ - 1, b, 2):
        prod = F(1)
        for i in range(a_, m):
            if (i - b) % 2 == 0:
                prod *= (2 * i - h) / (2 * i + 2 - h)
        q[m] = K + c - F(m, 3) - (K + c - F(a_ - 1, 3)) * prod
    for m in range(a_, b + 1, 2):
        prod = F(1)
        for i in range(m + 1, b):
            if (i - b) % 2:
                prod *= (2 * i + 2 - h) / (2 * i - h)
        q[m] = K + c + Q - F(m, 3) + (F(b, 3) - K - c) * prod
    p = {m: q[m] + q[m - 1] for m in range(a_, b + 1)}
    return p, q, K, Q


def equilibrium(n, h):
    k = block_of(n, h)
    b = n - 1 if n == t(k + 1, h) - 1 else n
    p, _, _, _ = interval_solution(k, b, h)
    return k, b, p


def gauss_rank(rows):
    m = [r[:] for r in rows]
    rank, col, nrows, ncols = 0, 0, len(m), len(m[0])
    while rank < nrows and col < ncols:
        piv = next((r for r in range(rank, nrows) if m[r][col] != 0), None)
        if piv is None:
            col += 1
            continue
        m[rank], m[piv] = m[piv], m[rank]
        for r in range(nrows):
            if r != rank and m[r][col] != 0:
                f = m[r][col] / m[rank][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[rank])]
        rank += 1
        col += 1
    return rank


def det(rows):
    m = [r[:] for r in rows]
    n, d = len(m), F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if m[r][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            m[c], m[piv] = m[piv], m[c]
            d = -d
        d *= m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            m[r] = [x - f * y for x, y in zip(m[r], m[c])]
    return d


# ---------------------------------------------------------------- sections 2-4
def check_game_basics():
    for h in (F(0), F(1), F(1, 3), F(7, 5), F(19, 10)):
        for i in range(-3, 61):
            for j in range(-3, 61):
                if i >= 1 and j >= 1:
                    check(a(i, j, h) == payoff_rule(i, j, h), "Lemma 2.1")
                check(a(j, i, h) == -a(i, j, h), "skew symmetry")
        p3 = {1: (4 - h) / (12 - 3 * h), 2: (6 - h) / (12 - 3 * h), 3: (2 - h) / (12 - 3 * h)}
        check(all(payoffs(p3, i, h) == 0 for i in (1, 2, 3)) and sum(p3.values()) == 1, "n = 3 equilibrium")
        check(interval_solution(1, 3, h)[0] == p3 and P(1, 3 - h / 2) == (2 - h) / (3 * (4 - h)), "n = 3 closed form")
    check([interval_solution(1, 3, F(0))[0][m] for m in (1, 2, 3)] == [F(1, 3), F(1, 2), F(1, 6)], "h = 0, n = 3")
    check([interval_solution(1, 3, F(1))[0][m] for m in (1, 2, 3)] == [F(1, 3), F(5, 9), F(1, 9)], "h = 1, n = 3")
    check(payoff_rule(3, 1, 1) == 5 and payoff_rule(2, 3, 1) == 3 and payoff_rule(1, 2, 1) == 1, "weighted RPS payments")
    check([t(k, 0) for k in range(1, 9)] == [3, 8, 22, 47, 89, 150, 236, 349], "first t_k(0)")
    check([t(k, 1) for k in range(1, 9)] == [3, 9, 22, 48, 89, 151, 236, 350], "first t_k(1)")
    for k in range(1, 100001):
        x = F(4 * k**3 + 5 * k, 6)
        ceil_x = -((-x.numerator) // x.denominator)
        check(t(k, 0) == ceil_x + 1 and t(k, 1) == x.numerator // x.denominator + 2, "floor/ceil forms")
        for h in (0, 1):
            check(F(t(k, h)) == F(2, 3) * k**3 + F(5, 6) * k + F(5 + 2 * h - (-1) ** (k + h), 4), "eq:tk-poly")
        check((4 * k**3 + 5 * k) % 6 == (0 if k % 2 == 0 else 3), "integer/half-integer pattern")


def check_rows_lemma():
    rng = random.Random(1)
    for h in (F(0), F(1), F(1, 2), F(3, 2)):
        for n in range(3, 25):
            for _ in range(4):
                raw = [rng.randint(0, 9) for _ in range(n)]
                if sum(raw) == 0:
                    raw[0] = 1
                p = {j + 1: F(v, sum(raw)) for j, v in enumerate(raw)}
                for i in range(1, n - 1):
                    X = sum((p[j] for j in p if payoff_rule(i, j, h) > 0), F(0))
                    lhs = payoffs(p, i + 2, h) - payoffs(p, i, h)
                    rhs = 4 * (X - p[i + 1]) + (2 * i + 4 - h) * (p[i] + p[i + 2]) - (4 * i + 2 - 2 * h) * p[i + 1]
                    check(lhs == rhs and X >= p[i + 1], "Lemma 3.1")


def check_kernel_and_even_intervals():
    for h in (F(0), F(1), F(1, 3), F(7, 5), F(9, 5)):
        for k in range(1, 6):
            for b in range(2 * k + 1, 2 * k + 21):
                I = list(range(b - 2 * k, b + 1))
                A = [[a(i, j, h) for j in I] for i in I]
                check(gauss_rank(A) == 2 * k, "Lemma 4.4: rank 2k")
                p, _, _, _ = interval_solution(k, b, h)
                check(sum(p.values()) == 1 and all(sum(A[r][c] * p[I[c]] for c in range(len(I))) == 0 for r in range(len(I))),
                      "Lemma 4.4: the Theorem 4.3 vector spans the kernel")
    # Even intervals: nonsingular for h in {0, 1}; Pfaffian formula (remark-level fact) for k <= 4
    for h in (F(0), F(1), F(1, 3), F(6, 5)):
        for k in range(1, 5):
            for a0 in range(1, 13):
                I = list(range(a0, a0 + 2 * k))
                d = det([[a(i, j, h) for j in I] for i in I])
                u = a0 - h / 2
                pf = 2**k * (F(2 * k * k + 1) * prod_(u + 2 * j for j in range(k)) - 2 * (k * k - 1) * prod_(u + 2 * j + 1 for j in range(k))) / 3
                check(d == pf * pf, "Pfaffian of an even interval")
    for h in H01:
        for k in range(1, 6):
            for a0 in range(1, 41):
                I = list(range(a0, a0 + 2 * k))
                check(det([[a(i, j, h) for j in I] for i in I]) != 0, "even intervals nonsingular for h in {0,1}")


def prod_(it):
    out = F(1)
    for v in it:
        out *= v
    return out


def check_interval_theorem():
    cases = 0
    for h in (F(0), F(1), F(1, 3), F(7, 5)):
        kmax, bmax = (10, 300) if h in H01 else (6, 60)
        for k in range(1, kmax + 1):
            L = 2 * k + 1
            for b in range(2 * k + 1, bmax + 1):
                a_ = b - 2 * k
                y = b - h / 2
                p, q, K, Q = interval_solution(k, b, h)
                check(sum(p.values()) == 1, "normalisation of the interval formulas")
                for i in range(a_, b + 1):
                    check(payoffs(p, i, h) == 0, "payoffs vanish on the interval")
                M = sum((j * pj for j, pj in p.items()), F(0))
                check(Q == sum(((-1) ** (b - j) * pj for j, pj in p.items()), F(0)), "Q = q_b")
                check(M - Q / 2 == F(a_ + b, 2) - F(1, 2 * L), "Lemma 4.1 moment equation")
                if b <= 2 * k + 40:
                    for j in range(a_, b + 1):
                        col = sum(a(i, j, h) for i in range(a_, b + 1))
                        check(col == L * F(a_ + b - 2 * j, 2) + (k if (j - b) % 2 == 0 else -(k + 1)), "Lemma 4.1 column sums")
                Ry = R(k, y)
                c = (1 - h) / 3
                check(p[b] == P(k, y), "p_b = P_k(y)")
                check(p[a_] == K + c + Q - F(a_, 3) + (F(b, 3) - K - c) * y * Ry / (y - 2 * k), "p_a formula")
                check(Ry * R(k, y + 1) == (y - 2 * k) / y, "eq:telescope")
                check(2 * y * p[b] - (2 * a_ - h) * p[a_] == (2 * k - 1) * Q + F(1, L), "eq:endpoints")
                qq = lambda i: F(0) if i < a_ else (q[i] if i <= b else (-1) ** ((i - b) % 2) * Q)
                for i in range(a_ - 3, b + 4):
                    lhs = payoffs(p, i, h) + payoffs(p, i + 1, h)
                    rhs = (2 * i + 1) - 2 * M + (-1) ** ((b - i) % 2) * Q + (2 * i + 2 - h) * qq(i + 1) - (2 * i - h) * qq(i - 1)
                    check(lhs == rhs, "Lemma 4.2")
                fa, fb = payoffs(p, a_ - 1, h), payoffs(p, b + 1, h)
                check(fa == (2 * a_ - 1) - 2 * M - Q + (2 * a_ - h) * p[a_], "f_{a-1} first form")
                check(fb == (2 * b + 1) - 2 * M - Q - 2 * y * p[b], "f_{b+1} first form")
                check(fa == 2 * y * p[b] - L * (1 + Q) and fb == L + F(1, L) - 2 * Q - 2 * y * p[b], "second forms")
                check(3 * L * fa == (2 * k * k + 4 * k + 3) * (2 * y - 4 * k - 2) - 2 * k * (k + 2) * 2 * y * Ry, "affine f_{a-1}")
                check(3 * L * fb == 2 * k * (k + 2) * (2 * y + 2) - (2 * k * k + 4 * k + 3) * 2 * y * Ry, "affine f_{b+1}")
                check(R(k + 1, y + 1) * Ry == (y - 2 * k - 1) / y and R(k + 1, y + 2) == y * Ry / (y + 1), "product identities")
                check(fa == F(2 * k + 3, L) * 2 * y * Ry * P(k + 1, y + 1), "Proposition 5.1, f_{a-1}")
                check(fb == -F(2 * k + 3, L) * (2 * y + 2) * P(k + 1, y + 2), "Proposition 5.1, f_{b+1}")
                cases += 1
    REPORT["interval_cases"] = cases


def check_example():
    h = F(0)
    p, _, _, Q = interval_solution(2, 10, h)
    check(R(2, F(10)) == F(8, 9) * F(6, 7) == F(16, 21), "Example 4.5 R_2(10)")
    check([p[m] for m in range(6, 11)] == [F(29, 315), F(92, 315), F(34, 105), F(74, 315), F(2, 35)], "Example 4.5 vector")
    M = sum((j * pj for j, pj in p.items()), F(0))
    check(Q == F(-17, 315) and M == F(496, 63) and M - Q / 2 == F(79, 10), "Example 4.5 Q, M")
    k, b, pe = equilibrium(10, h)
    check(k == 2 and b == 10 and pe == p and pe[10] == F(2, 35), "Example 4.5 is the equilibrium of G_10(0)")


def check_symbolic():
    import sympy as sp
    k, b, hs, r = sp.symbols("k b h R")
    L = 2 * k + 1
    a_ = b - 2 * k
    y = b - hs / 2
    K = (a_ + b - 1) / sp.Integer(2) - 1 / (2 * L)
    c = (1 - hs) / 3
    Qs = sp.symbols("Q")
    q_bm1 = K + c - (b - 1) / sp.Integer(3) - (K + c - (a_ - 1) / sp.Integer(3)) * r
    p_a = K + c + Qs - a_ / sp.Integer(3) + (b / sp.Integer(3) - K - c) * y * r / (y - 2 * k)
    eq = (2 * y) * (Qs + q_bm1) - (2 * a_ - hs) * p_a - (2 * k - 1) * Qs - 1 / L
    sol = sp.solve(sp.Eq(eq, 0), Qs)
    Qf = (2 * y * (r - 1) + 2 * k) / 3 + 1 / L
    check(len(sol) == 1 and sp.simplify(sol[0] - Qf) == 0, "Theorem 4.3: value of Q (symbolic)")
    check(sp.simplify(sp.diff(eq, Qs) - L) == 0, "coefficient of Q is L")
    Pk = ((2 * k**2 + 1) * r - 2 * (k**2 - 1)) / (3 * L)
    check(sp.simplify(Qf + q_bm1 - Pk) == 0, "Theorem 4.3: p_b = P_k(y) (symbolic)")
    fa = 2 * y * Pk - L * (1 + Qf)
    fb = L + 1 / L - 2 * Qf - 2 * y * Pk
    Pn = lambda rr: ((2 * (k + 1) ** 2 + 1) * rr - 2 * ((k + 1) ** 2 - 1)) / (3 * (2 * k + 3))
    check(sp.simplify(fa - (2 * k + 3) / L * 2 * y * r * Pn((y - 2 * k - 1) / (y * r))) == 0, "Prop 5.1 symbolic, f_{a-1}")
    check(sp.simplify(fb + (2 * k + 3) / L * (2 * y + 2) * Pn(y * r / (y + 1))) == 0, "Prop 5.1 symbolic, f_{b+1}")


def beta(k):
    lo, hi = mp.mpf(4 * k**3 + 5 * k + 3) / 6, mp.mpf(4 * k**3 + 5 * k + 6) / 6
    f = lambda yy: P(k, yy)
    check(f(lo) < 0 < f(hi), "Lemma 5.2 bracket")
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return (lo + hi) / 2


def check_root_bound():
    theta = {}
    for k in range(2, 401):
        bk = beta(k)
        theta[k] = bk - mp.mpf(4 * k**3 + 5 * k + 3) / 6
        check(0 < theta[k] < mp.mpf(1) / 2, "theta_k in (0, 1/2)")
        for h in H01:
            tk = t(k, int(h))
            check(P(k, tk - 1 - h / 2) < 0 < P(k, tk - h / 2), "first positive integer is t_k(h)")
            check(tk == int(mp.floor(bk + h / 2)) + 1, "t_k(h) = floor(beta_k + h/2) + 1")
        x0 = F(k * (4 * k * k - 1), 6)
        v = F(4 * k * k - 3, 12)
        check(x0 > 2 * v and x0 > k * k - k and x0 > k - F(1, 2), "x0 inequalities")
        check(sum((F(c * c) + F(1, 12) for c in [k + 1 - 2 * j for j in range(1, k + 1)]), F(0)) == k * v, "int_U t^2 dt = k v")
        check(F(6, 4 * k * k - 1) == 2 * F(3, 2 * k * k - 2) / (2 + F(3, 2 * k * k - 2)), "k/x0 = 2u/(2+u)")
        if k <= 8:
            hstar = 1 - 2 * theta[k] if k % 2 == 0 else 2 - 2 * theta[k]
            check(abs((bk + hstar / 2) - mp.nint(bk + hstar / 2)) < mp.mpf(10) ** -40, "exceptional set formula")
    check(all(P(1, F(y2, 2)) > 0 for y2 in range(5, 200)), "P_1 > 0 for y > 2")
    for k in range(2, 60):
        yv = mp.mpf(t(k, 0)) + mp.mpf(1) / 3
        x = yv - k - mp.mpf(1) / 2
        G = mp.fsum(mp.log((x + c + mp.mpf(1) / 2) / (x + c - mp.mpf(1) / 2)) for c in [k + 1 - 2 * j for j in range(1, k + 1)])
        check(abs(G + mp.log(R(k, yv))) < mp.mpf(10) ** -40, "G(x) = -log R_k(y)")
    REPORT["theta_2_10_400"] = [mp.nstr(theta[2], 6), mp.nstr(theta[10], 6), mp.nstr(theta[400], 6)]
    check(mp.nstr(theta[2], 4) == "0.1458" and mp.nstr(theta[10], 3) == "0.0372" and mp.nstr(theta[400], 2) == "0.00094",
          "theta values quoted in Remark 5.5")
    # Numerical validation of the newly written moment/mean-value proof.
    # This finite check does not prove its O(k^-3) conclusion.
    scaled_root_errors = [abs(k**3 * (theta[k] - mp.mpf(3) / (8 * k))) for k in range(2, 401)]
    check(max(scaled_root_errors) < mp.mpf(1) / 2,
          "finite check: k^3 |theta_k - 3/(8k)| < 1/2 for 2 <= k <= 400")
    REPORT["root_offset_k3_error_max_2_to_400"] = mp.nstr(max(scaled_root_errors), 12)
    for k in range(2, 401):
        delta_moment = F(k * (4 * k * k - 3), 12) - F(k**3, 12)
        check(delta_moment == F(k * (k * k - 1), 4), "root-offset second-moment difference")
    # Root bounds suffice for positivity of the first limiting mass.
    for k in range(2, 401):
        lower_beta_k = F(4 * k**3 + 5 * k + 3, 6)
        upper_beta_previous = F(4 * (k - 1)**3 + 5 * (k - 1) + 6, 6) if k >= 3 else F(2)
        check(lower_beta_k - 1 > upper_beta_previous, "root gap gives distinct exceptional equilibria")
    return theta


def check_structure_theorem(theta):
    for h in H01:
        for n in range(3, 601):
            k, b, p = equilibrium(n, h)
            a_ = b - 2 * k
            check(all(p[m] > 0 for m in range(a_, b + 1)), "positive support probabilities")
            for i in range(1, n + 1):
                f = payoffs(p, i, h)
                check(f == 0 if a_ <= i <= b else f < 0, "equilibrium conditions, strict off the support")
            expected = F(0) if n == t(k + 1, h) - 1 else P(k, n - h / 2)
            check(p.get(n, F(0)) == expected, "Theorem A value of P(n)")
    REPORT["structure_theorem_exact_n_le"] = 600
    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError:
        REPORT["independent_lp"] = "skipped (SciPy not installed)"
        return
    betas = {k: mp.mpf(4 * k**3 + 5 * k + 3) / 6 + theta[k] for k in range(2, 12)}

    def lp_matrix(n, h):
        return np.array([[float(payoff_rule(i, j, h)) for j in range(1, n + 1)] for i in range(1, n + 1)])

    def lp(n, h, objective=None):
        c = np.zeros(n) if objective is None else objective
        return linprog(c, A_ub=lp_matrix(n, h), b_ub=np.zeros(n), A_eq=np.ones((1, n)), b_eq=[1.0],
                       bounds=[(0, None)] * n, method="highs")

    worst = {}
    for hf in (0.0, 1.0, 0.2, 0.6, 1.3, 1.9):
        tk = lambda k: 3 if k == 1 else int(mp.floor(betas[k] + mp.mpf(hf) / 2)) + 1
        err = 0.0
        for n in range(3, 151):
            k = 1
            while tk(k + 1) <= n:
                k += 1
            pred = 0.0 if n == tk(k + 1) - 1 else float(P(k, mp.mpf(n) - mp.mpf(hf) / 2))
            res = lp(n, F(hf).limit_denominator(10**6))
            check(res.status == 0, "LP feasibility")
            err = max(err, abs(res.x[-1] - pred))
        check(err < 1e-7, "independent LP agrees with Theorem A for h = %s" % hf)
        worst[str(hf)] = err
    REPORT["independent_lp_n_le_150_max_abs_diff"] = worst
    # exceptional h: k = 2, h* = 1 - 2 theta_2, m = beta_2 + h*/2 = 8: G_7(h*) has a segment of equilibria
    hstar = float(1 - 2 * theta[2])
    n = 7
    def pay_float(i, j):
        if i == j:
            return 0.0
        wn = min(i, j) if (i - j) % 2 else max(i, j)
        return (2 * wn - hstar) if wn == i else -(2 * wn - hstar)

    A = np.array([[pay_float(i, j) for j in range(1, n + 1)] for i in range(1, n + 1)])
    e = np.zeros(n)
    e[-1] = 1.0
    lo = linprog(e, A_ub=A, b_ub=np.zeros(n), A_eq=np.ones((1, n)), b_eq=[1.0], bounds=[(0, None)] * n, method="highs")
    hi = linprog(-e, A_ub=A, b_ub=np.zeros(n), A_eq=np.ones((1, n)), b_eq=[1.0], bounds=[(0, None)] * n, method="highs")
    check(lo.status == 0 and hi.status == 0 and lo.x[-1] < 1e-7 and hi.x[-1] > 0.2, "non-uniqueness at an exceptional h")
    REPORT["exceptional_h_k2"] = {"h": hstar, "min_p7": float(lo.x[-1]), "max_p7": float(hi.x[-1])}


def check_block_table():
    expected = {1: (3, 6, 7, 4, "0.9055555556"), 2: (8, 20, 21, 13, "1.2051685518"), 3: (22, 45, 46, 24, "1.2449235895"),
                4: (47, 87, 88, 41, "1.3328846236"), 5: (89, 148, 149, 60, "1.3441658999"), 6: (150, 234, 235, 85, "1.3842642085")}
    for k, (lo, hi, zero, ell, Bk) in expected.items():
        B = sum((P(k, F(n)) for n in range(t(k, 0), t(k + 1, 0) - 1)), F(0))
        check((t(k, 0), t(k + 1, 0) - 2, t(k + 1, 0) - 1, t(k + 1, 0) - t(k, 0) - 1) == (lo, hi, zero, ell), "Table 1 ranges")
        check(str((Decimal(B.numerator) / Decimal(B.denominator)).quantize(Decimal(10) ** -10, rounding=ROUND_HALF_EVEN)) == Bk, "Table 1 B_k(0)")
    for h in (0, 1):
        for k in range(1, 1001):
            check(t(k + 1, h) - t(k, h) - 1 == 2 * k * k + 2 * k + (1 + (-1) ** (k + h)) // 2, "ell_k formula")


# ---------------------------------------------------------------- section 6
@functools.lru_cache(maxsize=None)
def poch(s, m):
    return F(1) if m == 0 else poch(s, m - 1) * (s + m - 1)


@functools.lru_cache(maxsize=None)
def bb_pmf(k, x):
    N = k - 1
    return comb(N, x) * poch(F(1, 2), x) * poch(F(3, 2), N - x) / poch(F(2), N)


def check_partial_fractions():
    check([w(2, j) for j in (1, 2)] == [3, 1], "w_2")
    check([w(3, j) for j in (1, 2, 3)] == [F(15, 4), F(3, 2), F(3, 4)], "w_3")
    check([w(4, j) for j in (1, 2, 3, 4)] == [F(35, 8), F(15, 8), F(9, 8), F(5, 8)], "w_4")
    for k in range(1, 31):
        check(sum(w(k, j) for j in range(1, k + 1)) == 2 * k, "sum of weights")
        for j in range(1, k + 1):
            check(w(k, j) / (2 * k) == bb_pmf(k, j - 1), "beta-binomial law")
        for y in [F(2 * k + 1), F(4 * k + 7, 2), F(2 * k + 10, 3)] + [F(n) - F(h, 2) for h in (0, 1) for n in range(2 * k + 1, 2 * k + 30, 4)]:
            check(1 - R(k, y) == sum((w(k, j) / 2 / (y - 2 * j + 1) for j in range(1, k + 1)), F(0)), "Lemma 6.1")
            expectation = sum((bb_pmf(k, x) / (2 * y - 2 - 4 * x) for x in range(k)), F(0))
            check(P(k, y) == F(1, 2 * k + 1) - F(2 * k * (2 * k * k + 1), 3 * (2 * k + 1)) * expectation, "eq:P-mixture")


def block_sum_mp(k, h, first=None, last=None):
    """Proposition 6.2 with a moving digamma window (two digamma values per block)."""
    first = t(k, h) if first is None else first
    last = t(k + 1, h) - 2 if last is None else last
    if last < first:
        return mp.mpf(0)
    half = mp.mpf(h) / 2
    harmonic = (mp.digamma(last - half) - mp.digamma(first - 1 - half)) / 2
    cs = [mp.mpf(1)]
    for j in range(1, k):
        cs.append(cs[-1] * mp.mpf(2 * j - 1) / (2 * j))
    weighted = mp.mpf(0)
    for j in range(1, k + 1):
        weighted += 2 * (2 * k - 2 * j + 1) * cs[j - 1] * cs[k - j] * harmonic
        if j < k:
            harmonic += (1 / mp.mpf(2 * first - h - 4 * j - 2) + 1 / mp.mpf(2 * first - h - 4 * j)
                         - 1 / mp.mpf(2 * last - h - 4 * j) - 1 / mp.mpf(2 * last - h - 4 * j + 2))
    return (last - first + 1 - mp.mpf(2 * k * k + 1) * weighted / 3) / (2 * k + 1)


def check_block_sums():
    Hm = lambda m: mp.digamma(m + 1) - mp.digamma(1)
    Hodd = lambda m: (mp.digamma(m + mp.mpf(1) / 2) - mp.digamma(mp.mpf(1) / 2)) / 2
    for m in range(0, 40):
        check(abs(Hodd(m) - mp.fsum(mp.mpf(1) / (2 * i - 1) for i in range(1, m + 1))) < mp.mpf(10) ** -45, "H^odd via digamma")
        check(abs(Hm(m) - mp.fsum(mp.mpf(1) / i for i in range(1, m + 1))) < mp.mpf(10) ** -45, "H via digamma")
    for h in (0, 1):
        for k in range(1, 15):
            exact = sum((P(k, F(n) - F(h, 2)) for n in range(t(k, h), t(k + 1, h) - 1)), F(0))
            ell = t(k + 1, h) - t(k, h) - 1
            br = lambda j: (mp.digamma(t(k + 1, h) - mp.mpf(h) / 2 - 2 * j) - mp.digamma(t(k, h) - mp.mpf(h) / 2 - 2 * j + 1))
            via = mp.mpf(ell) / (2 * k + 1) - mp.mpf(2 * k * k + 1) / (3 * (2 * k + 1)) * mp.fsum(
                mp.mpf(w(k, j).numerator) / w(k, j).denominator / 2 * br(j) for j in range(1, k + 1))
            check(abs(via - mp.mpf(exact.numerator) / exact.denominator) < mp.mpf(10) ** -40, "Proposition 6.2")
            for j in range(1, k + 1):
                alt = (Hm(t(k + 1, h) - 2 * j - 1) - Hm(t(k, h) - 2 * j)) if h == 0 else 2 * (Hodd(t(k + 1, h) - 2 * j - 1) - Hodd(t(k, h) - 2 * j))
                check(abs(br(j) - alt) < mp.mpf(10) ** -40, "bracket as harmonic numbers")
            check(abs(block_sum_mp(k, h) - via) < mp.mpf(10) ** -40, "moving window")
    check(sum((P(1, F(n)) for n in range(3, 7)), F(0)) == F(163, 180), "B_1(0)")
    S10 = sum((P(1, F(n)) for n in range(3, 7)), F(0)) + P(2, F(8)) + P(2, F(9)) + P(2, F(10))
    check(S10 == F(25493, 25200) and sum(equilibrium(n, F(0))[2].get(n, F(0)) for n in range(3, 11)) == S10, "S_0(10)")
    # Proposition 6.3
    for h in (0, 1):
        for k in range(1, 7):
            Wpoly = lambda x: mp.fsum(mp.mpf(w(k, j).numerator) / w(k, j).denominator * x ** (4 * (k - j)) for j in range(1, k + 1))
            for x in (mp.mpf("0.3"), mp.mpf("0.8"), mp.mpf(1)):
                integral = 4 * k / mp.pi * mp.quad(lambda th: mp.sqrt(th) / mp.sqrt(1 - th) * (1 - th + th * x**4) ** (k - 1), [0, 1])
                check(abs(integral - Wpoly(x)) < mp.mpf(10) ** -25, "W_k integral")
            ell = t(k + 1, h) - t(k, h) - 1
            integrand = lambda x: x ** (2 * t(k, h) - h - 4 * k + 1) * mp.fsum(x ** (2 * i) for i in range(ell)) * Wpoly(x)
            Bk = mp.mpf(ell) / (2 * k + 1) - mp.mpf(2 * k * k + 1) / (3 * (2 * k + 1)) * mp.quad(integrand, [0, 1])
            check(abs(Bk - block_sum_mp(k, h)) < mp.mpf(10) ** -30, "Proposition 6.3")
    for h in (0, 1):
        for k in range(1, 6):
            n = 10**8
            x = (P(k, F(n) - F(h, 2)) - F(1, 2 * k + 1)) * n
            check(abs(x + F(k * (2 * k * k + 1), 3 * (2 * k + 1))) < F(1, 10**6), "1/n term of P_k")


# ---------------------------------------------------------------- section 7
def moments_exact(k, order):
    pmf = [w(k, x + 1) / (2 * k) for x in range(k)]  # equals bb_pmf(k, x), checked for k <= 30
    return [sum((pmf[x] * x**i for x in range(k)), F(0)) for i in range(order + 1)]


def T_exact(k, h, order):
    eta = (-1) ** (k + h)
    ell = t(k + 1, h) - t(k, h) - 1
    aeta = F(5 * k, 3) + F(1 - eta, 2)
    mom = moments_exact(k, order)
    T = []
    for m in range(order + 1):
        total = F(0)
        for i in range(m + 1):
            shift_sum = sum(((aeta + 2 * r) ** (m - i) for r in range(ell)), F(0))
            total += comb(m, i) * (-4) ** i * mom[i] * shift_sum
        T.append(total)
    return T, ell, aeta


def check_asymptotic_algebra():
    for h in (0, 1):
        for k in range(1, 1001):
            eta = (-1) ** (k + h)
            check(2 * t(k, h) - h - 2 == F(4, 3) * k**3 + F(5, 3) * k + F(1 - eta, 2), "2 t_k - h - 2 = D + a_eta")
        for k in range(1, 13):
            D = F(4, 3) * k**3
            _, ell, aeta = T_exact(k, h, 0)
            for r in range(ell):
                e = sum((bb_pmf(k, x) / (1 + (aeta + 2 * r - 4 * x) / D) for x in range(k)), F(0))
                check(P(k, F(t(k, h) + r) - F(h, 2)) == F(1, 2 * k + 1) * (1 - (1 + F(1, 2 * k * k)) * e), "eq:P-geometric")
        for k in range(1, 41):
            T, ell, _ = T_exact(k, h, 2)
            mu = 2 * k * k + F(8, 3) * k + 1
            mom = moments_exact(k, 2)
            check(mom[1] == F(k - 1, 4) and 16 * (mom[2] - mom[1] ** 2) == k * k - 1, "E X and Var(4X)")
            check(T[0] == ell and T[1] == ell * mu and T[2] == ell * (mu * mu + F(ell * ell - 1, 3) + k * k - 1), "T_0, T_1, T_2")
        for k in (100, 200, 400):
            T, _, _ = T_exact(k, h, 3)
            check(abs(T[3] / F(k) ** 8 - 32) * k < 300, "T_3 = 32 k^8 + O(k^7)")
        for k in range(8, 300):
            eta = (-1) ** (k + h)
            ell = t(k + 1, h) - t(k, h) - 1
            aeta = F(5 * k, 3) + F(1 - eta, 2)
            xmax, xmin = aeta + 2 * (ell - 1), aeta - 4 * (k - 1)
            D = F(4, 3) * k**3
            check(xmax <= 4 * k * k + F(17, 3) * k and xmin >= -F(7, 3) * k + 4, "range of xi")
            check(max(abs(xmax), abs(xmin)) / D <= F(4, k) and ell <= 3 * k * k, "Lemma 7.1")


def coefficients(h):
    return [F(x) for x in CERT[h]["coefficients_even"]], [F(x) for x in CERT[h]["coefficients_odd"]]


def check_expansion():
    e0, o0 = coefficients(0)
    e1, o1 = coefficients(1)
    check(e0 == o1 and o0 == e1, "b_{j,eta} independent of h (eta = +1: even k for h = 0, odd k for h = 1)")
    table2 = [(F(3, 2), F(3, 2)), (F(-3, 4), F(-3, 4)), (F(3, 8), F(-3, 8)), (F(-9, 20), F(57, 40)),
              (F(291, 160), F(-219, 160)), (F(-2859, 560), F(-327, 280)), (F(6417, 560), F(17793, 2240)),
              (F(-203421, 8960), F(-168561, 8960)), (F(390009, 8960), F(274299, 8960))]
    check(all((e0[j], o0[j]) == table2[j] for j in range(9)), "Table 2 equals the certificate coefficients")
    zeta_coeff = {}
    for h, (E, O) in ((0, (e0, o0)), (1, (e1, o1))):
        u = lambda j: (E[j] + O[j]) / 2
        v = lambda j: (E[j] - O[j]) / 2
        zc = lambda j: u(j) - v(j) * (1 - F(2) ** (1 - j))   # coefficient of zeta(j)
        zeta_coeff[h] = [zc(2), zc(3), zc(4)]
        check(zc(2) * F(1, 6) == F((-1) ** (h + 1), 32), "pi^2/32 term, sign (-1)^(h+1)")
    check(zeta_coeff[0][1:] == [F(381, 320), F(-1497, 1280)] and zeta_coeff[1][1:] == [F(-69, 320), F(2073, 1280)],
          "zeta(3), zeta(4) contributions")
    REPORT["abs_b_j_eta_plus_minus"] = {j: [float(abs(e0[j])), float(abs(o0[j]))] for j in (10, 20, 30, 40, 50, 59)}
    for h in (0, 1):
        worst = 0
        for k in range(1, 1001):
            eta = (-1) ** (k + h)
            worst = max(worst, abs(block_sum_mp(k, h) - mp.mpf(3) / 2 + mp.mpf(3) / (4 * k) - 3 * eta * mp.mpf(1) / (8 * k * k)) * k**3)
        check(worst < 3, "Proposition 7.2(a)")
        REPORT["max_k3_remainder_7_2a_h%d" % h] = float(worst)
    for k in (400, 401):
        coeffs = e0 if k % 2 == 0 else o0
        approx = mp.fsum(mp.mpf(coeffs[j].numerator) / coeffs[j].denominator / mp.mpf(k) ** j for j in range(9))
        scaled = (block_sum_mp(k, 0) - approx) * mp.mpf(k) ** 9
        b9 = mp.mpf(coeffs[9].numerator) / coeffs[9].denominator
        check(abs(scaled - b9) < abs(b9) * mp.mpf("0.02"), "Table 2 checked against B_k(0) at k = 400, 401")
        REPORT["k9_residual_h0_k_%d" % k] = [float(scaled), float(b9)]
    for h in (0, 1):
        for k in (8, 9, 12, 17, 24, 33, 40):
            T, ell, _ = T_exact(k, h, 6)
            D = F(4, 3) * k**3
            Bk = block_sum_mp(k, h)
            for M in range(1, 7):
                bracket = ell - (1 + F(1, 2 * k * k)) * sum(((-1) ** m * T[m] / D**m for m in range(M + 1)), F(0))
                approx = mp.mpf(bracket.numerator) / bracket.denominator / (2 * k + 1)
                bound = mp.mpf(3) / 2 * k * (mp.mpf(4) / k) ** (M + 1) / (1 - mp.mpf(4) / k)
                check(abs(Bk - approx) <= bound, "Proposition 7.2(c)")


def c_tau(tau):
    return -(1 - 3 * tau * (1 - tau) * (4 * tau - 1)) / 4


def check_asymptotics():
    tex = TEX.read_text(encoding="utf-8")
    compact = re.sub(r"\\\\|[\s&]", "", tex)
    for h in (0, 1):
        cert = CERT[h]
        A = mp.mpf(cert["certified_80_decimals"]["A"])
        E = mp.mpf(cert["certified_80_decimals"]["sum_e"])
        check(abs(A - (mp.mpf(3) / 2 + 3 * mp.euler / 4 + mp.log(mp.mpf(3) / 2) / 4 - E)) < mp.mpf(10) ** -45, "series (12)")
        rem_text = re.search(r"0\.%s[0-9]{60,}" % ("5655" if h == 0 else "0794"), compact).group(0)
        rem = mp.mpf(cert["remainder_after_pi2_32"].strip("[").split()[0])
        check(abs(E - (-1) ** (h + 1) * mp.pi**2 / 32 - rem) < mp.mpf(10) ** -45, "remainder after the pi^2/32 term")
        check(abs(mp.mpf(rem_text) - rem) < mp.mpf(10) ** -(len(rem_text) - 3), "remainder digits in the text")
        S, worst_end = mp.mpf(0), 0
        for K in range(1, 1201):
            S += block_sum_mp(K, h)
            if K >= 10:
                N = t(K + 1, h) - 1
                d = S - (KAPPA * mp.cbrt(N) - mp.log(N) / 4 - A - mp.mpf(1) / (4 * K))
                worst_end = max(worst_end, abs(d) * K * K)
        check(worst_end < 10, "Theorem 7.3 at block ends")
        REPORT["max_K2_residual_block_ends_h%d" % h] = float(worst_end)
        worst_general, worst_p = 0, 0
        for K in (64, 128, 256, 512):
            S_prev = mp.fsum(block_sum_mp(k, h) for k in range(1, K))
            DK = t(K + 1, h) - t(K, h)
            for r in sorted(set([0, 1, DK // 7, DK // 3, DK // 2, (2 * DK) // 3, (6 * DK) // 7, DK - 2, DK - 1])):
                N = t(K, h) - 1 + r
                tau = mp.mpf(r) / DK
                S_N = S_prev + (block_sum_mp(K, h, t(K, h), N) if r >= 1 else 0)
                d = S_N - (KAPPA * mp.cbrt(N) - mp.log(N) / 4 - A - mp.mpf(3) / 2 * tau * (1 - tau) + c_tau(tau) / K)
                worst_general = max(worst_general, abs(d) * K * K)
                if 1 <= r <= DK - 1:
                    s = r
                    val = P(K, mp.mpf(N) - mp.mpf(h) / 2)
                    approx = mp.mpf(3 * s) / (4 * mp.mpf(K) ** 4) - mp.mpf(3 * s) / (8 * mp.mpf(K) ** 5) - mp.mpf(9 * s * s) / (8 * mp.mpf(K) ** 7)
                    worst_p = max(worst_p, abs(val - approx) * mp.mpf(K) ** 4)
        check(worst_general < 10 and worst_p < 30, "Theorem 7.3 for general N")
        REPORT["max_K2_residual_general_N_h%d" % h] = float(worst_general)
        REPORT["max_K4_residual_P_expansion_h%d" % h] = float(worst_p)


# ---------------------------------------------------------------- tables, digits and text
def check_tables_and_text():
    tex = TEX.read_text(encoding="utf-8")
    compact = re.sub(r"\\\\|[\s&]", "", tex)
    for h in (0, 1):
        check(CERT[h]["certified_80_decimals"]["A"] in compact, "certified digits of A(%d) in the text" % h)
    check(CERT[0]["certified_80_decimals"]["sum_e"] in compact, "certified digits of sum e_k(0)")
    check("1.7771799879" in tex and "1.6464176534" in tex, "abstract digits")
    # S_h(10^m) tables are published for h = 0 only (the values for h = 1 are a Project Euler answer)
    cpp = {}
    for line in CPP[0].read_text().split("\n"):
        if line.strip():
            e, value, _ = line.split()
            cpp[int(e)] = mp.mpf(value)
    arb30 = {int(m): mp.mpf(v) for m, v in CERT[0]["certified_30_decimals_S_10m"].items()}
    diff = max(abs(cpp[m] - arb30[m]) for m in range(1, 21))
    check(diff < mp.mpf("5e-31"), "C++ table agrees with Arb certificate (h = 0)")
    REPORT["max_diff_cpp_vs_arb_h0"] = mp.nstr(diff, 5)
    check("certified_30_decimals_S_10m" not in CERT[1], "no S_1 values are published")
    arb30 = {int(m): mp.mpf(v) for m, v in CERT[0]["certified_30_decimals_S_10m"].items()}
    A = mp.mpf(CERT[0]["certified_80_decimals"]["A"])
    rows = re.findall(r"^\$10\^\{(\d+)\}\$ & ([0-9.]+) & ([0-9.]+) & \$(-?[0-9.]+)\$ & \$(-?[0-9.]+)\$ \\\\$", tex, flags=re.M)
    check(len(rows) == 20, "Table 3 has 20 rows")
    diffs = {}
    for e, s_txt, ratio_txt, kd_txt, c_txt in rows:
        m = int(e)
        N = 10**m
        K = 1
        while t(K + 1, 0) - 1 <= N:
            K += 1
        tau = mp.mpf(N - t(K, 0) + 1) / (t(K + 1, 0) - t(K, 0))
        S = arb30[m]
        check(abs(mp.mpf(s_txt) - S) <= mp.mpf("5e-21") + mp.mpf("1e-30") and len(s_txt.split(".")[1]) == 20, "Table 3 S_0(N)")
        check(abs(mp.mpf(ratio_txt) - S / mp.cbrt(N)) <= mp.mpf("5e-13"), "Table 3 ratio")
        delta = S - (KAPPA * mp.cbrt(N) - mp.log(N) / 4 - A - mp.mpf(3) / 2 * tau * (1 - tau))
        check(abs(mp.mpf(kd_txt) - K * delta) <= mp.mpf("5e-7") and abs(mp.mpf(c_txt) - c_tau(tau)) <= mp.mpf("5e-7"), "Table 3 residuals")
        diffs[m] = float(abs(K * delta - c_tau(tau)))
        check(abs(K * delta - c_tau(tau)) * K < 3, "second-order residual is O(1/K)")
    REPORT["abs_K_delta_minus_c_tau_h0"] = diffs
    check(all(diffs[m] < 1e-5 for m in range(15, 21)), "text: below 1e-5 for m >= 15")
    check(mp.nstr(arb30[3] / 10, 4, strip_zeros=False) == "1.330" and mp.nstr(arb30[5] / mp.cbrt(10**5), 4, strip_zeros=False) == "1.613", "intro ratios")
    check(mp.nstr(KAPPA, 5) == "1.7171" and mp.nstr(KAPPA, 20).startswith("1.7170713"), "kappa")
    check("3,8,22,47,89,150,236,349" in compact and "3,9,22,48,89,151,236,350" in compact, "first thresholds in the text")


def main():
    start = time.time()
    steps = [check_game_basics, check_rows_lemma, check_kernel_and_even_intervals, check_interval_theorem, check_example,
             check_symbolic]
    for step in steps:
        t0 = time.time()
        step()
        print("%-34s ok  (%5.1f s, %d assertions so far)" % (step.__name__, time.time() - t0, COUNT["assertions"]), flush=True)
    t0 = time.time()
    theta = check_root_bound()
    print("%-34s ok  (%5.1f s, %d assertions so far)" % ("check_root_bound", time.time() - t0, COUNT["assertions"]), flush=True)
    t0 = time.time()
    check_structure_theorem(theta)
    print("%-34s ok  (%5.1f s, %d assertions so far)" % ("check_structure_theorem", time.time() - t0, COUNT["assertions"]), flush=True)
    for step in [check_block_table, check_partial_fractions, check_block_sums, check_asymptotic_algebra, check_expansion,
                 check_asymptotics, check_tables_and_text]:
        t0 = time.time()
        step()
        print("%-34s ok  (%5.1f s, %d assertions so far)" % (step.__name__, time.time() - t0, COUNT["assertions"]), flush=True)
    REPORT.update({"status": "PASS", "assertions": COUNT["assertions"], "seconds": round(time.time() - start, 1),
                   "python": platform.python_version(), "mpmath": mp.__version__})
    try:
        import sympy
        import scipy
        REPORT["sympy"], REPORT["scipy"] = sympy.__version__, scipy.__version__
    except ImportError:
        pass
    (HERE / "verify_paper_report.json").write_text(json.dumps(REPORT, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(REPORT, indent=2))


if __name__ == "__main__":
    sys.exit(main())
