# Publication package

Title: Parity Rock--Paper--Scissors Games and the Alper Constant

Author: Alper Ferudun. Affiliation: Mercury Software GmbH.
English unrefereed preprint; version 1.1; manuscript date 5 October 2026.
License: CC BY 4.0. Primary category math.OC; secondary math.CA.
Review status: Author reviewed; no independent peer review is claimed.
Suggested site path: /papers/parity-rps-alper-constant/ (PDF as paper.pdf).
Code and data: https://github.com/AlperTheKing/alper-constant

The author is shown alone at the top; affiliation, email and GitHub are
in a footnote. AI assistance is disclosed in running prose, not in a
separate section. No ORCID or other identity is invented.

main.tex is standalone and has an embedded bibliography.
references.bib contains matching reusable citation metadata and is not
required to compile this standalone source. Build with
`latexmk -pdf main.tex` (fig_equilibrium.pdf and fig_alper_function.pdf
must be next to main.tex). If latexmk's Perl runtime is unavailable, run
`pdflatex -interaction=nonstopmode -halt-on-error main.tex` until the
references and labels are stable (normally two or three passes). Copy the
resulting main.pdf to paper.pdf.

The paper studies the games G_n(h) in which a win with the number m pays
2m - h, 0 <= h < 2. The claims are: uniqueness, closed form and transition
law for every h outside a countable exceptional set E (Theorems 4.3 and
5.4, Theorem A), with explicit cubic thresholds for h = 0 and h = 1; the
asymptotic formula for S_h(N), h in {0,1}, with an explicit second-order
term (Theorem 7.3, Theorem B); 80 certified decimals of the Alper
constant A = A(0) and of A(1) (Theorem 8.2); and, new in version 1.1, the
Alper function (Section 10, Theorem C): the bound theta_k < 1/(2k) and the
resulting description of the thresholds for every h (Lemma 10.1), the
asymptotic formula for every h not in E with the extra term 3 phi/(4K),
phi = h - floor(h) (Theorem 10.3), continuity, strict monotonicity and
strict convexity between the points of E, the jumps at E and the limits at
1 and 2 (Proposition 10.4), the shift A(h+2) = A(h) for base parameters
0 <= h < 2 outside E (Proposition 10.5; not global periodicity),
the values at integers (Remark 10.6), and certified values of A(h) for
h = 0, 0.1, ..., 1.9 (Table 4, Figure 2). The case h = 1 is Project Euler
Problem 1012, which is cited together with the EulerSolve notes; no
priority is claimed for that case. A closed form of the constants is NOT
claimed. Values of S_1(N), which include the answer of the Project Euler
problem, are deliberately not published: the certificate computes the S
table only for h = 0 unless `--powers yes` is given.

## Files

- `paper.pdf`: compiled manuscript (21 pages).
- `main.tex`, `fig_equilibrium.pdf`, `fig_alper_function.pdf`,
  `references.bib`: source.
- `submission_metadata.md`, `verification_report.md`, `LICENSE.txt`.
- `release_build_receipt.json`: build and visual-inspection receipt.
- `source.zip`: everything above except `paper.pdf`, plus `reproducibility/`.
- `README.md`: landing page of the GitHub repository.
- `reproducibility/checks/verify_paper.py`: exact and high-precision
  checks of every displayed identity, example, table and asymptotic
  claim, including LP cross-checks for several h and the Section 10
  checks; writes `verify_paper_report.json`.
- `reproducibility/checks/review_function_high_precision.py`: a second
  numerical route for h = 0.3, 0.8, 1.3, 1.9, with 110-digit mpmath,
  separate digamma pairs, cutoff 128 and order 56. Shares the exact
  threshold and expansion generators; its roundoff is not interval
  certified. Writes the matching `.json` report.
- `reproducibility/certificate/alper_constant_certificate.py`: exact
  expansion coefficients and the Arb-certified evaluation of A(h) for
  h = 0 and h = 1, and of S_0(10^m), m <= 20; outputs `certificate_h0.json`,
  `certificate_h1.json`, `certificate_h0_cutoff128.json` (second parameter
  set) and `research_checks_h0.json`, `research_checks_h1.json`
  (independent mpmath sums, symbolic identities, Gosper and PSLQ logs).
- `reproducibility/alper_function/alper_function.py`: Arb-certified values
  of the Alper function (60 digits, exact thresholds for k <= 193, shifted
  expansion in the tail) and the first jump; outputs
  `alper_function_values.json` and `table_alper_rows.tex` (Table 4).
  The saved Arb balls retain the certified precision, and are read back
  to verify their error radii and the printed rounding.
  `make_alper_function_figure.py`: Figure 2 from 3207 certified
  evaluations on the intervals between the points of E with k <= 120;
  outputs `alper_function_grid.json` and `fig_alper_function.pdf`.
- `reproducibility/figures/alperconstant.cpp`: independent 100-digit C++
  evaluation of S_h(10^m); its output for h = 0 is `S_h0_40dp.txt`.
  `make_figure_and_tables.py`: Figure 1 and the rows of Tables 1-3.

## Reproduction

Tested with Python 3.12.4, python-flint 0.9.0, mpmath 1.3.0, SymPy 1.14.0,
SciPy 1.18.1 (optional; enables the independent LP checks), NumPy 2.2.6,
matplotlib 3.9.0, g++ 16.1.0 with Boost 1.91, MiKTeX 25.12
(pdfTeX 1.40.28). The current PDF was built by direct pdfLaTeX calls.

    cd reproducibility/certificate
    python alper_constant_certificate.py --h 0
    python alper_constant_certificate.py --h 1
    python alper_constant_certificate.py --h 0 --cutoff 128 --order 48 --digits 110 --output certificate_h0_cutoff128.json
    python alper_constant_certificate.py --research-checks certificate_h0.json --output research_checks_h0.json
    python alper_constant_certificate.py --research-checks certificate_h1.json --output research_checks_h1.json

    cd ../alper_function
    python alper_function.py
    python make_alper_function_figure.py

    cd ../figures
    g++ -O2 -std=c++17 -pthread alperconstant.cpp -o alperconstant
    ./alperconstant 0 > S_h0_40dp.txt
    python make_figure_and_tables.py

    cd ../checks
    python verify_paper.py
    python review_function_high_precision.py

Measured times on the preparation machine: each certificate about 3 s;
research checks under one minute each; alper_function.py 5 s;
make_alper_function_figure.py 39 s with 60 worker processes; C++ table
11 s with all 128 logical processors in the preparation run;
verify_paper.py 198.7 s (1,154,124 assertions) and
review_function_high_precision.py 1.8 s in the second review. The
original certificate/C++ timings are historical; these unchanged files
were compared with the backup in this review. Assertions must remain
enabled. Finite checks do not replace the written proofs.

The earlier package records version 1.0 as published at
https://github.com/AlperTheKing/alper-constant on 5 October 2026;
version 1.1 adds Section 10 and has now received a second local review.
This review did not upload or publish the revised files. The remote
release and DOI status were not rechecked; verify the actual public
record before announcing it.
