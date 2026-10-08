# The Alper constant

Code, certificates and paper for

**Parity Rock–Paper–Scissors Games and the Alper Constant**
Alper Ferudun, 5 October 2026 ([paper.pdf](paper.pdf); [paper page on eulersolve.org](https://eulersolve.org/papers/parity-rps-alper-constant/)).

## The games

Two players choose numbers from {1, …, n}. Equal numbers draw, an odd
difference is won by the smaller number, a nonzero even difference by the
larger number, and a win with the number m pays 2m − h (0 ≤ h < 2).
For h = 1 this is the game of
[Project Euler Problem 1012](https://projecteuler.net/problem=1012);
for h = 0 the payment is proportional to the winning number.

The paper proves that the equilibrium is unique and solves it in closed
form for every h outside a countable exceptional set E, and it determines
the asymptotics of S(N) = Σ_{n=3}^{N} P(n), where P(n) is the equilibrium
probability of the largest number:

S(N) = (3/2)^{4/3} N^{1/3} − ¼ log N − 𝒜(h) − (3/2) τ(1 − τ) − (1 − 3φ − 3τ(1 − τ)(4τ − 1))/(4K) + O_h(N^{−2/3}),

for each fixed h in [0,2) outside E, where τ is the relative position of
N in its block and φ = h − ⌊h⌋.

## The constant

The Alper constant is the Euler-type constant of the game with h = 0:

    𝒜 = lim (κ N_K^{1/3} − ¼ log N_K − S(N_K)),   κ = (3/2)^{4/3},

along the block ends N_K. Its value correctly rounded to 80 decimal
places, certified with Arb ball arithmetic, is

    𝒜 ≈ 1.77717998790542318852801821627634156242557103207026437429288159344562194194456353

For h = 1 the same construction gives, also rounded to 80 decimal places,
𝒜(1) ≈ 1.64641765344056003385448411598581328630299683869991753007336025510008825020074857.

## The Alper function

The same limit exists for every h outside E and defines the Alper
function h ↦ 𝒜(h), whose value at 0 is the Alper constant (Section 10 of
the paper). It is continuous, strictly increasing and strictly convex
between consecutive points of E, and jumps down at each of them; the
first point is 6 − 2√7 = 0.70849…, with a jump of (√7 − 1)/6. It tends to
𝒜(1) as h → 1⁻ and to 𝒜(0) as h → 2⁻. The extension satisfies
𝒜(h + 2) = 𝒜(h) for base parameters h in [0,2) outside E. This is not
global periodicity: for example, 𝒜(5) = 𝒜(1) − 1.

| h | 𝒜(h) | h | 𝒜(h) |
|---|---|---|---|
| 0.0 | 1.7771799879 | 1.0 | 1.6464176534 |
| 0.1 | 1.8116524222 | 1.1 | 1.6854092798 |
| 0.2 | 1.8465611039 | 1.2 | 1.7252229864 |
| 0.3 | 1.8819309085 | 1.3 | 1.7659309340 |
| 0.4 | 1.9177891039 | 1.4 | 1.8076154543 |
| 0.5 | 1.9541656688 | 1.5 | 1.8503710032 |
| 0.6 | 1.9910936650 | 1.6 | 1.8943065925 |
| 0.7 | 2.0286096750 | 1.7 | 1.9395488435 |
| 0.8 | 1.7913637669 | 1.8 | 1.8508040685 |
| 0.9 | 1.7105887967 | 1.9 | 1.8174408673 |

(Rounded to 10 decimals; the Arb enclosures have radius below 3·10⁻⁵⁴.)

Review status: **Author reviewed**. This is an unrefereed preprint; no
independent peer review or proof-assistant verification is claimed.

## Reproduction

| Path | Purpose |
|---|---|
| `reproducibility/certificate/alper_constant_certificate.py` | exact expansion coefficients and certified values (`--h 0`, `--h 1`) |
| `reproducibility/alper_function/alper_function.py` | certified values of the Alper function (Table 4) |
| `reproducibility/alper_function/make_alper_function_figure.py` | graph of the Alper function (Figure 2) |
| `reproducibility/figures/alperconstant.cpp` | independent 100-digit C++ evaluation of S(10^m) |
| `reproducibility/figures/make_figure_and_tables.py` | Figure 1 and Tables 1–3 of the paper |
| `reproducibility/checks/verify_paper.py` | exact and high-precision checks of every identity, example and table |

See [README_submission.md](README_submission.md) for the commands,
versions and timings, and [verification_report.md](verification_report.md)
for what was checked.

## License

CC BY 4.0, see [LICENSE.txt](LICENSE.txt).
