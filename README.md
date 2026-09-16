# AF3 ligand chirality checker

A lightweight Python/RDKit command-line tool that compares a ligand's observed
3D tetrahedral stereochemistry with an authoritative isomeric SMILES.

**Status: research prototype, version 0.1.0.** This is a standalone validator,
not an upstream AlphaFold patch, a structure repair tool, or a predictor of binding
affinity. It does not modify the input coordinates. It is not affiliated with
Google DeepMind or the RDKit project.

## Install

Tested with Python 3.11 and RDKit 2026.03.6. From this repository directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

On Windows, activate with `.venv\Scripts\activate`. Alternatively, create a
Conda environment with Python 3.11 and RDKit 2026.03.6, then run the Python script
directly. No GPU or AlphaFold installation is required for checking coordinates.

## Run

Supply a **ligand-only PDB**, not the full protein complex:

```bash
af3-chirality-check \
  --ligand-pdb predicted_ligand.pdb \
  --reference-smiles 'C[C@H](N)C(=O)O' \
  --af3-smiles 'C[C@H](N)C(=O)O' \
  --report chirality_report.json
```

`python af3_ligand_chirality_check.py` accepts the same arguments.
`--af3-smiles` is optional: provide the chemical definition actually supplied to
AF3, not a guessed copy of the reference. A coordinate-only PDB cannot independently
establish protonation and tautomer state. Standard AF3 mmCIF output must first be
converted and filtered to one ligand instance; general mmCIF extraction is not
included in this release.

The checker restores bond orders using the reference, enumerates graph mappings,
infers stereochemistry from observed coordinates, and assigns accurate CIP labels
(including pseudoasymmetric r/s). It screens local geometry and writes JSON.
One complete, reliable mapping must match all reference centres to produce PASS.
Matches from different mappings are never combined. The mapping search is capped.

## Interpret the report

| Verdict | Interpretation |
| --- | --- |
| `PASS` | A complete mapping matches every defined centre and passes the local geometry checks. |
| `CHIRALITY_MISMATCH` | All enumerated mappings disagree with the defined reference chirality, with assessable geometry. |
| `LOW_GEOMETRY_QUALITY` | The local geometry or stereochemical assignment is insufficient for a reliable comparison. |
| `AMBIGUOUS_MAPPING` | The mapping cap was exceeded, or mixed inconclusive outcomes remain without a complete match. |
| `CHEMICAL_STATE_MISMATCH` | The supplied chemical definitions or mapped graphs are incompatible. |
| `NO_DEFINED_STEREOCENTERS` | No assigned tetrahedral R/S/r/s centres were found in the reference. |
| `INPUT_ERROR` | Input parsing or coordinate-graph mapping failed. |

Exit status 0 means a report was produced, **not necessarily PASS**. Inspect its
`verdict` field. Input errors normally return 1; argument errors return 2.
Reported atom indices are zero-based RDKit indices. `selected_mapping_index`
identifies a matching mapping, not unique atom-name correspondence.

## Reproduce the tests without AF3 data

```bash
python -m unittest discover -s tests -v
python tests/stress_chirality.py --output validation/stress.json
```

The tests generate synthetic coordinates themselves. The stress runner returns
nonzero on an incorrect result or failed fixture generation.

Validation to date: 15 regression tests, 560 synthetic cases, and two user-run
AF3 example complexes. See [validation summary](VALIDATION.md) for counts and
limitations. No naturally occurring AF3 chirality error has yet been demonstrated
in these two examples. The synthetic mirrored negative controls are explicitly
artificial.

## Limitations

- Initial scope: single, complete, ordinary organic ligands with defined
  tetrahedral stereochemistry. E/Z, atropisomerism, metal stereochemistry and
  general peptide-chain validation are outside the claimed scope.
- Unspecified reference stereocentres are not validated. PASS applies only to
  the centres that were defined and assessed.
- PDB distance-based connectivity can fail for distorted coordinates. Mapping
  failure alone is not evidence of a chirality inversion.
- Template restoration cannot independently verify the ligand's chemical state.
- Local geometry thresholds are heuristic and uncalibrated against a large
  experimental benchmark. They do not certify global geometry or pocket fit.
- This release has no batch interface or general full-complex/mmCIF importer.
- Synthetic tests use RDKit for both generation and assessment. They are not an
  independent accuracy estimate; the 560 cases are correlated transformations.

## Contributing

Useful contributions include full-complex/mmCIF extraction, independent labelled
benchmarks, geometry-threshold calibration, and reproducible edge cases. In an
issue, include the RDKit version, command, reference chemical definition, expected
result and a minimal input you have permission to share. Label deliberately
altered structures as synthetic. Do not include private paths or credentials.

## References and license

- [RDKit](https://www.rdkit.org/)
- [AlphaFold 3](https://github.com/google-deepmind/alphafold3)
- [MTX reference](https://www.rcsb.org/ligand/MTX)
- [BTN reference](https://www.rcsb.org/ligand/BTN)

Original project code is released under the [MIT license](LICENSE).
Dependencies retain their own licenses. AF3 outputs, MSAs, templates and the local
development environment are not distributed in this repository.
