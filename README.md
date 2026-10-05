# The Alper constant

Code, certificates and paper for

**Parity Rock–Paper–Scissors Games and the Alper Constant**
Alper Ferudun, 5 October 2026 ([paper.pdf](paper.pdf)).

## The games

Two players choose numbers from {1, …, n}. Equal numbers draw, an odd
difference is won by the smaller number, a nonzero even difference by the
larger number, and a win with the number m pays 2m − h (0 ≤ h < 2).
For h = 1 this is the game of
[Project Euler Problem 1012](https://projecteuler.net/problem=1012);
for h = 0 the payment is proportional to the winning number.

The paper proves that the equilibrium is unique and solves it in closed
form for every h outside a countable exceptional set, and it determines
the asymptotics of S(N) = Σ_{n=3}^{N} P(n), where P(n) is the equilibrium
probability of the largest number:

S(N) = (3/2)^{4/3} N^{1/3} − ¼ log N − 𝒜(h) − (3/2) τ(1 − τ) − (1 − 3τ(1 − τ)(4τ − 1))/(4K) + O(N^{−2/3}),

where τ is the relative position of N in its block (h = 0, 1).

## The constant

The Alper constant is the Euler-type constant of the game with h = 0:

    𝒜 = lim (κ N_K^{1/3} − ¼ log N_K − S(N_K)),   κ = (3/2)^{4/3},

along the block ends N_K. Its value correctly rounded to 80 decimal
places, certified with Arb ball arithmetic, is

    𝒜 ≈ 1.77717998790542318852801821627634156242557103207026437429288159344562194194456353

For h = 1 the same construction gives, also rounded to 80 decimal places,
𝒜(1) ≈ 1.64641765344056003385448411598581328630299683869991753007336025510008825020074857.

Review status: **Author reviewed**. This is an unrefereed preprint; no
independent peer review or proof-assistant verification is claimed.

## Reproduction

| Path | Purpose |
|---|---|
| `reproducibility/certificate/alper_constant_certificate.py` | exact expansion coefficients and certified values (`--h 0`, `--h 1`) |
| `reproducibility/figures/alperconstant.cpp` | independent 100-digit C++ evaluation of S(10^m) |
| `reproducibility/figures/make_figure_and_tables.py` | figure and tables of the paper |
| `reproducibility/checks/verify_paper.py` | exact and high-precision checks of every identity, example and table |

See [README_submission.md](README_submission.md) for the commands,
versions and timings, and [verification_report.md](verification_report.md)
for what was checked.

## License

CC BY 4.0, see [LICENSE.txt](LICENSE.txt).
