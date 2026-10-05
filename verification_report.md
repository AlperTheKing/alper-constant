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
Remark 5.5 now also proves theta_k = 3/(8k) + O(k^-3), by comparing
second moments in the root-bound integrals and applying the mean value
theorem. Thus the assertions about accumulation at 1 and 2 and eventual
agreement of the thresholds no longer rely on numerical evidence.

## Computer checks

`reproducibility/checks/verify_paper.py` was rerun during review and passed
1,118,462 assertions in 122.3 s (Python 3.12.4, mpmath 1.3.0,
SymPy 1.14.0, SciPy 1.18.1):

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

Both main certificates, the second h = 0 parameter set, and the mpmath
and PSLQ checks were rerun for this review. The 80-place results and
exact coefficients agree with the packaged certificates. The C++
program was recompiled with g++ 16.1.0 and run for h = 0; all 40 printed
decimal places at the 20 test inputs agree with the packaged output.

A further re-check after the review rebuilt the PDF with latexmk (no
warnings, 16 pages), reran verify_paper.py on the final source (1,118,462
assertions passed, 124.5 s), reran both main certificates and the second
h = 0 parameter set from scratch (identical to the packaged certificates
apart from run times), recompiled and reran `alperconstant.cpp` for h = 0
(identical output), and inspected all 16 pages again.

These are independent computational routes; the general theorems rest
on the written proofs. They do not constitute independent peer review.
The certification is a computer-assisted evaluation of proved error
bounds, not a formal proof. Values of S_1(N) are deliberately not
published, because they include the answer of Project Euler Problem 1012.

## Corrections made in the review

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
  of the constant to all nonexceptional h is explicitly an open task.
- Kept numerical tables and figures ahead of the concluding section and
  bibliography, and added PDF title, author and subject metadata. Table 3
  and Figure 1 later received more flexible float placement, so that they
  share one page.

## Review status

**Author reviewed**. This label records the author-directed review
and corrections in this revision; it is not a claim of external human
refereeing or a separately documented manual read-through by the author.

The final PDF was compiled with pdfLaTeX (MiKTeX 25.12) without warnings,
overfull boxes or unresolved references, and all 16 pages were inspected
visually during the AI-assisted review. The built-in editor's compiler
returned an environment error ("Unable to find standard directories for
platform"); the PDF was successfully exported with the installed
MiKTeX compiler instead. There is no independent refereeing or
proof-assistant verification.

Scope: the case h = 1 is Project Euler Problem 1012, which is cited; its
solvers may have obtained parts of the h = 1 results, and no priority is
claimed for that case. A closed form of the Alper constant remains open.

The reviewed revision was uploaded to the public GitHub repository
https://github.com/AlperTheKing/alper-constant on 5 October 2026. No DOI has been
assigned. The source archive and build receipt correspond to this PDF.
