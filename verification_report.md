# Verification report

Manuscript: Parity Rock--Paper--Scissors Games and the Alper Constant
Author: Alper Ferudun, Mercury Software GmbH. Version 1.0, 5 October 2026.

## What is proved in the text

For the games G_n(h), in which a win with m pays 2m - h, the structure
theorem combines four steps: a row-comparison lemma that forces interval
supports; an explicit solution of the equilibrium equations on intervals
of odd length; a kernel lemma showing that these intervals have a
one-dimensional kernel with nonzero coordinate sum; and closed forms of
the two neighbouring payoffs. Together with a root bound proved by an
integral comparison, these exclude even supports for every h outside a
countable exceptional set and give the transition law. Positivity of the
closed-form probabilities is never assumed; it is inherited from the
equilibrium whose existence the minimax theorem guarantees. For h = 0 and
h = 1 the asymptotic theorem, including the second-order term
-(1-3tau(1-tau)(4tau-1))/(4K), is proved from a beta-binomial mixture
representation and a geometric expansion with explicit error bounds.

## Computer checks

`reproducibility/checks/verify_paper.py` passed 1,117,664 assertions in
139 s (Python 3.12.4, mpmath 1.3.0, SymPy 1.14.0, SciPy 1.18.1):

- Lemma 2.1 against the verbal payoff rule and the n = 3 solution for
  h in {0, 1, 1/3, 7/5, 19/10}; Lemma 3.1 for random rational strategies.
- Lemma 4.4: rank 2k of every odd interval matrix tested (h in
  {0, 1, 1/3, 7/5, 9/5}, k <= 5), spanned by the Theorem 4.3 vector.
  Even interval matrices are nonsingular for h in {0, 1} (k <= 5,
  a <= 40); their Pfaffians follow an explicit product formula for k <= 4.
- Theorem 4.3, Lemmas 4.1-4.2 and Proposition 5.1 in exact arithmetic:
  6,416 interval systems (h = 0 and h = 1 with k <= 10, b <= 300; h = 1/3
  and h = 7/5 with k <= 6, b <= 60), plus symbolic SymPy proofs of the
  value of Q, of p_b = P_k(y) and of the neighbour payoffs for general h.
- Lemma 5.2 for 2 <= k <= 400, the thresholds t_k(0) and t_k(1), the
  description of the exceptional set, and theta_k ~ 3/(8k).
- Theorem 5.4 in exact arithmetic for every 3 <= n <= 600 and h in
  {0, 1}: positive support probabilities, zero payoffs on the support and
  strictly negative payoffs off it.
- Independent floating-point linear programming (SciPy HiGHS) for
  3 <= n <= 150 and h in {0, 1, 0.2, 0.6, 1.3, 1.9}, without using any
  formula of the paper: maximum deviation from Theorem A 9.1e-12. At the
  exceptional value h = 1 - 2 theta_2 = 0.70850 the game G_7(h) has a
  segment of equilibria (P(7) ranges over [0, 0.2743]), as stated in
  Remark 5.5.
- Section 6 (partial fractions, beta-binomial law, block sums with
  ordinary and odd harmonic numbers, integral representation) and
  Section 7 (equation (10) exactly for k <= 12, T_0, T_1, T_2 for
  k <= 40, Lemma 7.1, Proposition 7.2(a) for k <= 1000 with k^3 times the
  remainder below 1.43, Table 2 against B_400 and B_401, Proposition
  7.2(c)) for h = 0 and h = 1.
- Theorem 7.3 for h = 0 and h = 1: K^2 times the residual stays below
  0.20 at all block ends 10 <= K <= 1200 and below 1.69 at 36 points of the
  blocks K = 64, 128, 256, 512.
- Every entry of Tables 1-3 and the numbers quoted in the text.

## Certified constants

`reproducibility/certificate/alper_constant_certificate.py` (python-flint
0.9.0, Arb balls at 160 digits) sums 191 exact head blocks and an
expansion tail through order 59 with proved bounds 3.481e-100 and
1.26e-109, for h = 0 and for h = 1. The enclosure radius is below
7.4e-100 for both constants, and the 80-place rounding is checked to be
constant on each ball. For h = 0 a second parameter set (cutoff 128,
order 48, 110 digits; radius below 6e-71) overlaps the first. Independent
130-digit mpmath summations with a separate digamma pair for every pole
differ from the certified midpoints by 7.7e-94 (h = 0) and 4.3e-94
(h = 1). S_0(10^m), m <= 20, is certified to 30 decimals; the independent
100-digit C++/Boost program `alperconstant.cpp` agrees to within 4.1e-31.
PSLQ in 120-digit arithmetic (tolerance 1e-85, coefficients up to 1000)
found no relation for A(h), A(h) - log(3/2)/4 or the deviation sums, for
h = 0 and h = 1, against two bases of classical constants.

These are independent computational routes, NOT an independent human or
model review, and the finite checks do not prove the general theorems.
The certification is a computer-assisted evaluation of proved error
bounds, not a formal proof. Values of S_1(N) are deliberately not
published, because they include the answer of Project Euler Problem 1012.

## Review status

The final PDF was compiled with pdfLaTeX (MiKTeX 25.12) without warnings,
overfull boxes or unresolved references, and all 16 pages were inspected
visually by the preparing AI agent. Author review is pending. There is no
independent refereeing or proof-assistant verification.

Scope: the case h = 1 is Project Euler Problem 1012, which is cited; its
solvers may have obtained parts of the h = 1 results, and no priority is
claimed for that case. A closed form of the Alper constant remains open.

No DOI or public release is asserted by this local verification report.
