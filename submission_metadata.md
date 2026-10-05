# Submission metadata

- Title: Parity Rock--Paper--Scissors Games and the Alper Constant
- Author: Alper Ferudun
- Affiliation: Mercury Software GmbH
- Email: alper@mercurycodelab.com
- GitHub: https://github.com/AlperTheKing
- Code and data: https://github.com/AlperTheKing/alper-constant
- Primary category: math.OC
- Secondary category: math.CA
- Language: English
- Resource type: Publication / Preprint
- Version: 1.1 (adds Section 10, the Alper function)
- Manuscript date: 2026-10-05
- License: Creative Commons Attribution 4.0 International (CC BY 4.0)
- Keywords: rock-paper-scissors; zero-sum game; Nash equilibrium; skew-symmetric matrix; beta-binomial distribution; digamma function; asymptotic expansion; interval arithmetic; Alper constant; Alper function
- Related problem: Project Euler Problem 1012 (the case h = 1), cited in the paper
- Suggested site path: /papers/parity-rps-alper-constant/
- DOI: Not assigned
- Publication status: Unrefereed preprint; reviewed revision prepared locally
- Review status: Author reviewed
- Latest review: Section 10 reviewed and corrected locally; 21-page PDF, 1,154,124 checks passed, table/grid regenerated, separate 110-digit cross-check added
- Release scope: Version 1.0 (paper, code, certificates) is public at https://github.com/AlperTheKing/alper-constant; version 1.1 prepared locally; site publication pending

## Abstract

Two players choose numbers from {1,...,n}. Equal numbers draw, an odd
difference is won by the smaller number, a nonzero even difference is won
by the larger number, and a win with the number m pays 2m-h, where
0 <= h < 2 is fixed. For h = 1 this is the game of Problem 1012 of Project
Euler; for h = 0 the payment is proportional to the winning number. For
every h outside an explicit countable set E we prove that the game has a
unique equilibrium, supported on 2k+1 consecutive numbers ending at n or
n-1, and we solve it in closed form. When the support ends at n, the
probability P(n) of the largest number is
((2k^2+1)R_k - 2(k^2-1))/(3(2k+1)) for an explicit product R_k
that depends only on the largest payment; otherwise P(n) = 0.
The support widens at
thresholds given by the zeros of these functions; for h in {0,1} the
thresholds are cubic polynomials in k for each parity. For S_h(N) = sum_{n=3}^N P(n)
we prove S_h(N) = (3/2)^{4/3} N^{1/3} - (1/4) log N - A(h)
- (3/2) tau(1-tau) + O_h(N^{-1/3}), where tau in [0,1) is the relative
position of N in its block, and we determine the next term. The constant
A(h) plays the role that Euler's constant plays for the harmonic series.
We call A = A(0) = 1.7771799879... the Alper constant, give a series and an
integral representation, and certify 80 decimal places of A and of
A(1) = 1.6464176534... with interval arithmetic. No expression of either
constant in classical constants was found. The Alper function h -> A(h) is
strictly increasing and convex between consecutive points of E and jumps
down at each of them; it tends to A(1) as h -> 1- and to A(0) as h -> 2-,
and A(h+2) = A(h) for base parameters 0 <= h < 2 outside E when the
games are extended to larger parameters.

## Public-notice scope

Complete theorems for the family G_n(h): uniqueness, closed form and
transition law for every h outside a countable exceptional set, with
explicit cubic thresholds for h = 0 and h = 1. Rigorous asymptotics of
S_h(N) for every h not in E with an explicit second-order term. Definition
of the Alper constant A = A(0) as an Euler-type limit, with 80 certified
decimals (computer-assisted interval evaluation of proved error bounds).
The Alper function A(h): proved shape (monotone and convex pieces, jumps at
E, limits at 1 and 2, and A(h+2) = A(h) for 0 <= h < 2 outside E) and
certified values at h = 0.1 steps. Global periodicity is not claimed.
The case h = 1 is Project Euler Problem 1012, which is cited; its solvers
may have obtained parts of the h = 1 results, and no priority is claimed
for that case. NOT claimed: a closed form, irrationality or transcendence
of the constants; the negative PSLQ searches do not exclude a closed form. Values of
S_1(N) are deliberately not published. Unrefereed; no proof-assistant
verification.

The release-scope entry records the earlier package's publication status;
the remote state and DOI were not rechecked in this second local review.
No upload or publication was performed by this review.
Post only after actual public record and site links have been verified.
