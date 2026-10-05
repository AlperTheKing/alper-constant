"""Figure 2 of the paper: the Alper function on [0, 2), drawn separately on every interval between
consecutive points of the exceptional set E (first points: k <= 120).  Each value is an Arb ball
from alper_function.py; the script also checks that the sampled values increase and are convex on
every interval, and writes alper_function_grid.json and ../../fig_alper_function.pdf.

Run:  python make_alper_function_figure.py
"""
import json
from math import ceil
from multiprocessing import Pool
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
KMAX = 120
mp.mp.dps = 30


def theta(k):
    """beta_k - (4k^3+5k+3)/6 by bisection on the bracket of Lemma 5.2."""
    def P(y):
        r = mp.mpf(1)
        for j in range(1, k + 1):
            r *= (y - 2 * j) / (y - 2 * j + 1)
        return (2 * k * k + 1) * r - 2 * (k * k - 1)
    lo, hi = mp.mpf(4 * k**3 + 5 * k + 3) / 6, mp.mpf(4 * k**3 + 5 * k + 6) / 6
    assert P(lo) < 0 < P(hi)
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if P(mid) < 0 else (lo, mid)
    return (lo + hi) / 2 - mp.mpf(4 * k**3 + 5 * k + 3) / 6


def evaluate(h_pair):
    from flint import fmpq
    import alper_function as af
    num, den = h_pair
    A, _, _ = af.alper(fmpq(num, den))
    return float(A.mid()), float(A.rad())


def main():
    points = {k: (k % 2) + 1 - 2 * theta(k) for k in range(2, KMAX + 1)}
    E0 = [float(points[k]) for k in range(2, KMAX + 1, 2)]
    E1 = [float(points[k]) for k in range(3, KMAX + 1, 2)]
    assert all(x < y for x, y in zip(E0, E0[1:])) and all(x < y for x, y in zip(E1, E1[1:]))
    components = [(0.0, E0[0], True)] + [(x, y, False) for x, y in zip(E0, E0[1:])]
    components += [(1.0, E1[0], True)] + [(x, y, False) for x, y in zip(E1, E1[1:])]
    den = 10**12
    jobs, owner = [], []
    for c, (a, b, closed) in enumerate(components):
        n = max(3, ceil((b - a) * 1500))
        nudge = (b - a) * 1e-4
        for i in range(n):
            x = a + (b - a) * i / (n - 1)
            if i == 0 and not closed:
                x = a + nudge
            if i == n - 1:
                x = b - nudge
            jobs.append((round(x * den), den))
            owner.append(c)
    with Pool(60) as pool:
        values = pool.map(evaluate, jobs, chunksize=4)
    grid = [{"a": a, "b": b, "h": [], "A": [], "radius_max": 0.0} for a, b, _ in components]
    for (num, d), c, (val, rad) in zip(jobs, owner, values):
        grid[c]["h"].append(num / d)
        grid[c]["A"].append(val)
        grid[c]["radius_max"] = max(grid[c]["radius_max"], rad)
    # increasing and convex on every interval (sampled values)
    for g in grid:
        hs, As = g["h"], g["A"]
        assert all(y > x for x, y in zip(As, As[1:])), (g["a"], g["b"])
        slopes = [(As[i + 1] - As[i]) / (hs[i + 1] - hs[i]) for i in range(len(hs) - 1)]
        assert all(s2 > s1 - 1e-9 for s1, s2 in zip(slopes, slopes[1:])), (g["a"], g["b"])
    data = {"kmax": KMAX, "theta": {str(k): mp.nstr(theta(k), 15) for k in range(2, KMAX + 1)},
            "exceptional_points": {str(k): mp.nstr(points[k], 15) for k in range(2, KMAX + 1)},
            "components": grid, "evaluations": len(jobs)}
    (HERE / "alper_function_grid.json").write_text(json.dumps(data) + "\n", encoding="utf-8")
    print("evaluations:", len(jobs), " max radius:", max(g["radius_max"] for g in grid))
    print("min sampled value: %.10f at h = %.6f" % min((A, h) for g in grid for h, A in zip(g["h"], g["A"])))
    print("max sampled value: %.10f at h = %.6f" % max((A, h) for g in grid for h, A in zip(g["h"], g["A"])))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.size": 9})
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    for k in range(2, 12):
        ax.axvline(float(points[k]), color="#cccccc", lw=0.5, zorder=0)
    for g in grid:
        ax.plot(g["h"], g["A"], color="#1f4e8c", lw=0.9)
    A0, A1 = grid[0]["A"][0], grid[len(E0)]["A"][0]
    ax.plot([0, 1], [A0, A1], "o", color="#c0504d", ms=3.5, zorder=3)
    ax.plot([2], [A0], "o", mfc="white", mec="#c0504d", ms=3.5, zorder=3)
    ax.annotate(r"$\mathcal{A}=\mathcal{A}(0)$", (0, A0), textcoords="offset points", xytext=(5, -11), fontsize=8)
    ax.annotate(r"$\mathcal{A}(1)$", (1, A1), textcoords="offset points", xytext=(6, -3), fontsize=8)
    ax.set_xlabel(r"$h$ (a win with $m$ pays $2m-h$)")
    ax.set_ylabel(r"$\mathcal{A}(h)$")
    ax.set_xlim(0, 2.03)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2])
    fig.tight_layout()
    fig.savefig(HERE.parent.parent / "fig_alper_function.pdf")
    print("wrote alper_function_grid.json and fig_alper_function.pdf")


if __name__ == "__main__":
    main()
