# Validation status

Results obtained during development on 2026-09-15 with RDKit 2026.03.6 and
Python 3.11; rerun before initial publication on 2026-09-16.

## Reproducible synthetic testing

- 15 regression tests: correct and mirrored coordinates, single-centre inversion,
  pseudoasymmetric inversion, symmetry mapping, explicit/implicit hydrogen
  handling, atom-order changes, planar geometry and chemical-state mismatch.
- 560 stress cases: 20 ligand structures, four random seeds and seven variations.
- Expected results: 408 PASS and 152 CHIRALITY_MISMATCH; no inconclusive outcomes
  after correcting the symmetric-polyol case.

The suite includes a meso polyol whose full mirror is chemically equivalent. Its
mirror should pass; changing only the pseudoasymmetric centre should fail.
These cases guard against accepting all symmetric molecules indiscriminately.

## Real AF3 examples (aggregate observations)

| User-run complex | Models | Ligand instances | Centres examined | Result |
| --- | ---: | ---: | ---: | --- |
| E. coli DHFR / methotrexate | 5 | 5 | 5 | All matched reference |
| Streptavidin / biotin | 5 | 20 | 60 | All matched reference |

Each complex used one seed; models and equivalent binding sites are correlated.
No natural AF3 inversion was found. A separate intentionally reflected copy of
the AF3 methotrexate ligand changed S to R and was correctly flagged while all
interatomic distances remained unchanged. This is an injected negative control,
not a failed AF3 prediction. Raw AF3 files are not included here, so these real
examples are reported observations rather than bundled reproducible fixtures.

## Interpretation

The prototype detects controlled tetrahedral chirality errors and handles these
two real examples. These results do not establish general sensitivity,
specificity, an AF3 error rate, or reliability for untested stereochemical classes.
Larger, independently labelled datasets remain needed.
