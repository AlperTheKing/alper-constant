# Publication package

Title: Parity Rock--Paper--Scissors Games and the Alper Constant

Author: Alper Ferudun. Affiliation: Mercury Software GmbH.
English unrefereed preprint; version 1.0; manuscript date 5 October 2026.
License: CC BY 4.0. Primary category math.OC; secondary math.CA.
Suggested site path: /papers/parity-rps-alper-constant/ (PDF as paper.pdf).
Code and data: https://github.com/AlperTheKing/alper-constant

The author is shown alone at the top; affiliation, email and GitHub are
in a footnote. AI assistance is disclosed in running prose, not in a
separate section. No ORCID or other identity is invented.

main.tex is standalone and has an embedded bibliography.
references.bib contains matching reusable citation metadata and is not
required to compile this standalone source. Build with
`latexmk -pdf main.tex` (fig_equilibrium.pdf must be next to main.tex).

The paper studies the games G_n(h) in which a win with the number m pays
2m - h, 0 <= h < 2. The claims are: uniqueness, closed form and transition
law for every h outside a countable exceptional set (Theorems 4.3 and
5.4, Theorem A), with explicit cubic thresholds for h = 0 and h = 1; the
asymptotic formula for S_h(N), h in {0,1}, with an explicit second-order
term (Theorem 7.3, Theorem B); and 80 certified decimals of the Alper
constant A = A(0) and of A(1) (Theorem 8.2). The case h = 1 is Project
Euler Problem 1012, which is cited together with the EulerSolve notes;
no priority is claimed for that case. A closed form of the constants is
NOT claimed. Values of S_1(N), which include the answer of the Project
Euler problem, are deliberately not published: the certificate computes
the S table only for h = 0 unless `--powers yes` is given.

## Files

- `paper.pdf`: compiled manuscript (16 pages).
- `main.tex`, `fig_equilibrium.pdf`, `references.bib`: source.
- `submission_metadata.md`, `verification_report.md`, `LICENSE.txt`.
- `release_build_receipt.json`: build and visual-inspection receipt.
- `source.zip`: everything above except `paper.pdf`, plus `reproducibility/`.
- `README.md`: landing page of the GitHub repository.
- `reproducibility/checks/verify_paper.py`: exact and high-precision
  checks of every displayed identity, example, table and asymptotic
  claim, including LP cross-checks for several h; writes
  `verify_paper_report.json`.
- `reproducibility/certificate/alper_constant_certificate.py`: exact
  expansion coefficients and the Arb-certified evaluation of A(h) for
  h = 0 and h = 1, and of S_0(10^m), m <= 20; outputs `certificate_h0.json`,
  `certificate_h1.json`, `certificate_h0_cutoff128.json` (second parameter
  set) and `research_checks_h0.json`, `research_checks_h1.json`
  (independent mpmath sums, symbolic identities, Gosper and PSLQ logs).
- `reproducibility/figures/alperconstant.cpp`: independent 100-digit C++
  evaluation of S_h(10^m); its output for h = 0 is `S_h0_40dp.txt`.
  `make_figure_and_tables.py`: Figure 1 and the rows of Tables 1-3.

## Reproduction

Tested with Python 3.12.4, python-flint 0.9.0, mpmath 1.3.0, SymPy 1.14.0,
SciPy 1.18.1 (optional; enables the independent LP checks), NumPy 2.2.6,
matplotlib 3.9.0, g++ 16.1.0 with Boost 1.91, MiKTeX 25.12 (pdfTeX 4.23,
latexmk 4.88).

    cd reproducibility/certificate
    python alper_constant_certificate.py --h 0
    python alper_constant_certificate.py --h 1
    python alper_constant_certificate.py --h 0 --cutoff 128 --order 48 --digits 110 --output certificate_h0_cutoff128.json
    python alper_constant_certificate.py --research-checks certificate_h0.json --output research_checks_h0.json
    python alper_constant_certificate.py --research-checks certificate_h1.json --output research_checks_h1.json

    cd ../figures
    g++ -O2 -std=c++17 -pthread alperconstant.cpp -o alperconstant
    ./alperconstant 0 > S_h0_40dp.txt
    python make_figure_and_tables.py

    cd ../checks
    python verify_paper.py

Measured times on the preparation machine: each certificate about 3 s;
research checks under one minute each; C++ table 11 s with all 128
logical processors; verify_paper.py 139 s. Assertions must remain
enabled. Finite checks do not replace the written proofs.

No DOI, Zenodo record or site page exists for this paper at local
preparation. Review the manuscript before publishing; never announce a
reserved DOI as a public release.
