# Verification report

Manuscript: Parity Rock--Paper--Scissors Games and the Alper Constant
Author: Alper Ferudun, Mercury Software GmbH. Version 1.1, 5 October 2026.

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
Remark 5.5 now also proves theta_k = 3/(8k) + O(k^-3), by comparing
second moments in the root-bound integrals and applying the mean value
theorem. Thus the assertions about accumulation at 1 and 2 and eventual
agreement of the thresholds no longer rely on numerical evidence.

Version 1.1 adds Section 10, the Alper function. Lemma 10.1 proves
theta_k < 1/(2k) for every k >= 2 by the integral comparison of Lemma 5.2
at x0 + 1/(2k), so t_k(h) = t_k(floor h) for k >= 1/(1 - phi). Theorem 10.3
transfers the proofs of Proposition 7.2 and Theorem 7.3 to every h in
[0,2) outside E: the regular blocks have the offset a_eta - phi, which changes B_k by
-3 phi/(4k^2) and the asymptotic formula by +3 phi/(4K). Proposition 10.4
proves continuity, strict monotonicity and strict convexity between the
points of E (each P_k(n - h/2) is decreasing and concave in h by the
beta-binomial mixture), the jump -sum P_{k-1}(beta_k - 1) at each point
of E, and the limits A(1) and A(0) at 1- and 2- (dominated convergence
with a uniform O(k^-2) bound). Proposition 10.5 relabels G_n(h+2j) as
G_{n-j}(h) plus strictly losing numbers, which gives A(h+2) = A(h) for
base parameters h in [0,2) outside E. This is not global periodicity:
A(5) = A(1) - 1. Remark 10.6 treats larger shifts and positive integers (odd payments give odd
Pfaffians and unique equilibria; at even h >= 4 a zero payment creates
two pure equilibria).

## Computer checks

`reproducibility/checks/verify_paper.py` was rerun in this second review
and passed **1,154,124 assertions in 198.7 s** (Python 3.12.4,
mpmath 1.3.0, SymPy 1.14.0, SciPy 1.18.1). This includes the Section 10
checks and 60 new checks on saved interval certificates:

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
- Lemma 5.2 for 2 <= k <= 400, the thresholds t_k(0) and t_k(1), and the
  description of the exceptional set. The new checks validate the exact
  second-moment difference and root-gap inequality used in Remark 5.5;
  the largest observed k^3 |theta_k - 3/(8k)| is 0.333990. The asymptotic
  assertion is established by the written proof, not this finite check.
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
- Section 10 (version 1.1): theta_k < 1/(2k) and theta_k decreasing for
  2 <= k <= 400; the inequalities of the proof of Lemma 10.1 in exact
  arithmetic for 2 <= k <= 3000 and G(x0 + 1/(2k)) < lambda_k for k <= 60;
  the threshold rule of Lemma 10.1 by exact sign tests for h = j/20,
  0 <= j < 40, and k <= 60. Theorem 10.3 for h = 0.3, 0.8, 1.3, 1.9
  (h = 0.8 has an irregular threshold at k = 2; h = 1.9 at k = 3, 5, 7): k^3 times
  the remainder of the B_k expansion stays below 3.11 for k <= 800; K^2
  times the residual stays below 1.36 at the block ends 10 <= K <= 800
  and below 1.30 at 28 points of the blocks K = 64, 128, 256, 512,
  whereas without the term 3 phi/(4K) it grows to 180-541; Richardson
  extrapolation of the exact partial sums gives A(h) independently, within
  9e-12 of the certified values. Proposition 10.4 on the 3207 certified
  values of Figure 2: increasing and convex on each of the 119 intervals,
  the points of E for k <= 120, and the first jump against (sqrt 7 - 1)/6
  (difference 4e-13 at distance 1e-12 from 6 - 2 sqrt 7), with SymPy for
  beta_2 = 5 + sqrt 7. Proposition 10.5 in exact arithmetic: the raised
  equilibria of G_n(h + 2j) for h in {0, 1, 1/3, 7/5, 3/4, 19/10},
  j <= 3 and j + 3 <= n <= j + 30, with strictly negative payoffs off the
  support, P(3) = 0 for G_3(h + 2) and c_2(h) = 1. Remark 10.6: every
  even principal minor is odd for h in {1, 3, 5, 7, 9}, n <= 8;
  c_j(1) = 0, 1, 14/9, 91/45 for j <= 4 by exact support enumeration; two
  pure equilibria of G_{j+1}(2j) for j = 2, ..., 5. Table 4 and the
  numbers of Section 10 against the certified values. The 20 saved Arb
  balls are read back at 100 digits: each radius is below 3e-54, and
  rounding is constant on the whole ball at both 10 and 20 decimal places.

## Certified constants and Alper function values

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

`reproducibility/alper_function/alper_function.py` (version 1.1) evaluates
A(h) in 60-digit Arb arithmetic with 191 exact head blocks, thresholds
from exact sign tests for k <= 193, Lemma 10.1 beyond, and the shifted
expansion through order 39 with the same proved truncation bounds. The
20 values of Table 4 have saved enclosure radius below 3e-54; the values
at h = 0 and h = 1 overlap the 80-digit certificates. All 20 values and
all 3207 Figure 2 grid values were recomputed in this second review.
The 10- and 20-place roundings and the threshold lists of the 20 table
values agree exactly with the pre-review version.

The second review also added
`reproducibility/checks/review_function_high_precision.py`. It uses
110-digit mpmath arithmetic and a separate digamma pair for every
partial-fraction pole, with cutoff 128 and expansion order 56 (the table
uses cutoff 192 and order 40). For h = 0.3, 0.8, 1.3, 1.9 the differences
from the saved Arb midpoints are at most 2.14e-55. The proved algebraic
tail bound is below 1.061e-83; this bound excludes mpmath roundoff, so the
extra digits printed by this cross-check are not presented as certified
digits. Exact threshold, moment and expansion generators are shared
with the certificate code; the numerical summation route and parameters
are different. The check passed in 1.8 s.

The original version 1.0 review had rerun both 80-place certificates,
the second h = 0 parameter set, the mpmath and PSLQ checks, and the
100-digit C++/Boost computation (g++ 16.1.0); the packaged outputs
matched. Its final re-check passed 1,118,462 assertions and inspected
16 PDF pages. These older certificate generators/results and C++
source/output are byte-for-byte unchanged in the present review. They
were compared with the pre-review backup, but were not freshly rerun
as part of this Section 10 review.

These are cross-checks with the shared components specified above; the general theorems rest
on the written proofs. They do not constitute independent peer review.
The certification is a computer-assisted evaluation of proved error
bounds, not a formal proof. Values of S_1(N) are deliberately not
published, because they include the answer of Project Euler Problem 1012.

## Corrections from the original review

- Replaced the numerical inference about exceptional-set accumulation
  with a proof of the root-offset asymptotic. Also showed explicitly
  that the two limiting exceptional equilibria are distinct.
- Qualified odd-support and uniqueness statements by h not in E, and
  separated the general structure theorem from the asymptotic theorem
  proved for h = 0 and h = 1. The abstract now includes P(n) = 0 when the
  support ends at n - 1 and specifies the parity dependence of thresholds.
- Corrected the sign in the narrative describing discovery of A(1).
  Corrected the figure caption: P(n) vanishes one step before the
  support widens.
- Distinguished values rounded to 80 decimal places from literal
  decimal prefixes followed by an ellipsis, in the manuscript and README.
- Removed an unverified novelty claim and clarified that failed PSLQ
  searches do not establish nonexistence of a closed form. Extension
  of the constant to all nonexceptional h was then an open task; version
  1.1 carries it out in Section 10 and replaces the open question by the
  monotonicity of theta_k and the behaviour near 1 and 2.
- Kept numerical tables and figures ahead of the concluding section and
  bibliography, and added PDF title, author and subject metadata. Table 3
  and Figure 1 later received more flexible float placement, so that they
  share one page.

## Version 1.1

Section 10 (the Alper function), Theorem C, the related sentences of the
abstract, introduction and concluding remarks, Table 4 and Figure 2 were
added after the original review. The second review checked their written
arguments and reran the computations described above. It made these
corrections and clarifications:

- Made the uniform O(k^-2) estimate in Proposition 10.4(c) explicit,
  including the extra block endpoint r = ell_k, which lies outside the
  earlier lemma's range. This supplies a common summable bound for
  dominated convergence at h = 1 and h = 2.
- Restricted A(h+2) = A(h) explicitly to base parameters in [0,2)
  outside E throughout the abstract and supporting documentation.
  Remark 10.6 now states the range N >= j+2 for the finite initial-game
  correction, qualifies the odd-integer statement as positive, and
  explicitly rules out a global-periodicity interpretation.
- Gave the short proof that 6 - 2 sqrt(7) is the first exceptional
  point, and made the fixed-h dependence and uniformity in tau of the
  asymptotic error explicit.
- Fixed a certificate serialization defect: `A.str(30)` discarded
  precision and widened the saved interval to about 1e-30, while a
  separate field reported an internal radius of order 1e-54. The code
  now saves the full ball, reads it back, checks containment and the
  radius, and computes the printed rounding from that saved ball.
  The displayed table values did not change.
- Added the separate 110-digit numerical cross-check and the 60
  saved-ball checks; corrected the report's irregular-parameter list
  and its Section 10 block-check cutoff from 1000 to 800.
- Kept Table 3 and Figure 1 before Section 10 so that they no longer
  interrupt the new lemma's proof. Rebuilt the manuscript with the
  installed pdfLaTeX, without warnings, box problems, missing glyphs or
  unresolved references, and visually inspected all **21 pages**.

## Review status

**Author reviewed**. This label records the author-directed review
and corrections in this revision; it is not a claim of external human
refereeing or a separately documented manual read-through by the author.

The final version 1.1 PDF was compiled with pdfLaTeX (pdfTeX 1.40.28,
MiKTeX 25.12) without warnings, overfull or underfull boxes, missing glyphs
or unresolved references; all 21 pages were inspected visually during
the AI-assisted review. The built-in editor's compiler
returned an environment error ("Unable to find standard directories for
platform"); the PDF was successfully exported with the installed
MiKTeX compiler instead. There is no independent refereeing or
proof-assistant verification.

Scope: the case h = 1 is Project Euler Problem 1012, which is cited; its
solvers may have obtained parts of the h = 1 results, and no priority is
claimed for that case. A closed form of the Alper constant remains open.

The earlier package records publication of version 1.0 at
https://github.com/AlperTheKing/alper-constant on 5 October 2026.
The remote release and DOI status were not rechecked in this local
review. This further-reviewed version 1.1 was prepared locally; this
review did not upload or publish it. The source archive and build
receipt correspond to the corrected 21-page PDF.
