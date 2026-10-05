"""Figure 1 and the LaTeX table rows of the paper, for the game G_n(0) (a win with m pays 2m).

Inputs: ../certificate/certificate_h0.json (certified constant, certified S_0(10^m), expansion
coefficients) and S_h0_40dp.txt (independent 100-digit C++ values, written by alperconstant.cpp).
"""
from decimal import Decimal, ROUND_HALF_EVEN
from fractions import Fraction as F
from pathlib import Path
import json
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
mp.mp.dps = 60
CERT = json.loads((HERE.parent / "certificate" / "certificate_h0.json").read_text(encoding="utf-8"))
KAPPA = mp.cbrt(mp.mpf(81) / 16)
ALPER = mp.mpf(CERT["certified_80_decimals"]["A"])


def t(k):
    """First n of the block k for h = 0: ceil((4k^3+5k)/6) + 1."""
    return (4 * k**3 + 5 * k + 5) // 6 + 1


def P_float(n):
    """Equilibrium probability of the largest number for h = 0 (Theorem A)."""
    k = 1
    while t(k + 1) <= n:
        k += 1
    if n == t(k + 1) - 1:
        return 0.0
    r = 1.0
    for j in range(1, k + 1):
        r *= (n - 2 * j) / (n - 2 * j + 1)
    return ((2 * k * k + 1) * r - 2 * (k * k - 1)) / (3 * (2 * k + 1))


def P_exact(k, n):
    r = F(1)
    for j in range(1, k + 1):
        r *= F(n - 2 * j, n - 2 * j + 1)
    return ((2 * k * k + 1) * r - 2 * (k * k - 1)) / (3 * (2 * k + 1))


def fixed(x, decimals):
    return str(Decimal(mp.nstr(x, 60, min_fixed=-80, max_fixed=80)).quantize(Decimal(10) ** -decimals, rounding=ROUND_HALF_EVEN))


S = {int(m): mp.mpf(v) for m, v in CERT["certified_30_decimals_S_10m"].items()}
S_cpp = {}
for line in (HERE / "S_h0_40dp.txt").read_text().split("\n"):
    if line.strip():
        e, v, _ = line.split()
        S_cpp[int(e)] = mp.mpf(v)
assert max(abs(S[m] - S_cpp[m]) for m in S) < mp.mpf("5e-31")

# Table 3: N, S_0(N) to 20 decimals, S_0(N)/N^(1/3) to 12 decimals, K*delta(N) and c(tau) to 6 decimals.
rows = []
for e in range(1, 21):
    N = 10**e
    K = 1
    while t(K + 1) - 1 <= N:
        K += 1
    tau = mp.mpf(N - t(K) + 1) / (t(K + 1) - t(K))
    delta = S[e] - (KAPPA * mp.cbrt(N) - mp.log(N) / 4 - ALPER - mp.mpf(3) / 2 * tau * (1 - tau))
    c = -(1 - 3 * tau * (1 - tau) * (4 * tau - 1)) / 4
    ratio = S[e] / mp.power(10, mp.mpf(e) / 3)
    rows.append(r"$10^{%d}$ & %s & %s & $%s$ & $%s$ \\" % (
        e, fixed(S[e], 20), fixed(ratio, 12), fixed(K * delta, 6), fixed(c, 6)))
(HERE / "table_S_rows.tex").write_text("\n".join(rows) + "\n", encoding="utf-8")

# Table 1: the first blocks
block_rows = []
for k in range(1, 7):
    B = sum((P_exact(k, n) for n in range(t(k), t(k + 1) - 1)), F(0))
    block_rows.append(r"%d & %d & %d--%d & %d & %d & %s\\" % (
        k, t(k), t(k), t(k + 1) - 2, t(k + 1) - 1, t(k + 1) - t(k) - 1, fixed(mp.mpf(B.numerator) / B.denominator, 10)))
(HERE / "table_blocks_rows.tex").write_text("\n".join(block_rows) + "\n", encoding="utf-8")


# Table 2: coefficients b_{j,eta}, eta = (-1)^(k+h); for h = 0 the even-k column is eta = +1
def tex_fraction(x):
    f = F(x)
    if f.denominator == 1:
        return "$%d$" % f.numerator
    return "$%s%d/%d$" % ("-" if f < 0 else "", abs(f.numerator), f.denominator)


coef_rows = ["%d & %s & %s \\\\" % (j, tex_fraction(CERT["coefficients_even"][j]), tex_fraction(CERT["coefficients_odd"][j]))
             for j in range(9)]
(HERE / "table_b_rows.tex").write_text("\n".join(coef_rows) + "\n", encoding="utf-8")

# Figure: (a) P(n) sawtooth, (b) S_0(10^m)/10^(m/3) against the limit and the asymptotic curve
plt.rcParams.update({"font.family": "serif", "font.size": 9})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.7))
ns = list(range(3, 501))
ax1.plot(ns, [P_float(n) for n in ns], lw=0.8, color="#1f4e8c")
for k in range(2, 8):
    ax1.axvline(t(k) - 1, color="#bbbbbb", lw=0.5, zorder=0)
ax1.set_xlabel(r"$n$")
ax1.set_ylabel(r"$P(n)$")
ax1.set_title(r"(a) probability of the largest number", fontsize=9)
ax1.set_xlim(0, 500)
ax1.set_ylim(0, 0.3)
es = list(range(3, 21))
ax2.plot(es, [float(S[e] / mp.power(10, mp.mpf(e) / 3)) for e in es], "o", ms=3, color="#1f4e8c",
         label=r"$S(N)/N^{1/3}$, $N=10^m$")
xs = [2.5 + i * 0.05 for i in range(351)]
ax2.plot(xs, [float(KAPPA - (mp.log(mp.power(10, x)) / 4 + ALPER) / mp.power(10, mp.mpf(x) / 3)) for x in xs],
         lw=0.8, color="#c0504d", label=r"$\kappa-(\frac{1}{4}\log N+\mathcal{A})/N^{1/3}$")
ax2.axhline(float(KAPPA), ls="--", lw=0.8, color="black", label=r"$\kappa=(3/2)^{4/3}$")
ax2.set_xlabel(r"$m=\log_{10}N$")
ax2.set_ylim(1.30, 1.75)
ax2.set_xticks([3, 5, 10, 15, 20])
ax2.set_xlim(2.5, 20.5)
ax2.set_title(r"(b) approach to the limit", fontsize=9)
ax2.legend(fontsize=7, loc="lower right", frameon=False)
fig.tight_layout()
fig.savefig(HERE.parent.parent / "fig_equilibrium.pdf")
print("wrote table_S_rows.tex, table_blocks_rows.tex, table_b_rows.tex and fig_equilibrium.pdf")
print("\n".join(block_rows))
print("\n".join(coef_rows))
print("\n".join(rows))
