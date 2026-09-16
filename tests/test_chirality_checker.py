"""Regression tests for the AF3 ligand chirality checker.

Run with an RDKit-enabled interpreter:
    python -m unittest tests/test_chirality_checker.py
"""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


SCRIPT = Path(__file__).resolve().parents[1] / "af3_ligand_chirality_check.py"
SMILES = "C[C@H](F)Cl"


def write_ligand_pdb(path: Path, mirror: bool = False, smiles: str = SMILES,
                     mode: str = "normal") -> None:
    """Create a known stereochemical structure, optionally reflected in x."""
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    result = AllChem.EmbedMolecule(molecule, randomSeed=7, enforceChirality=True)
    if result != 0:
        raise RuntimeError("Could not embed synthetic test ligand")
    AllChem.MMFFOptimizeMolecule(molecule)
    if mirror:
        conformer = molecule.GetConformer()
        for atom_index in range(molecule.GetNumAtoms()):
            position = conformer.GetAtomPosition(atom_index)
            position.x = -position.x
            conformer.SetAtomPosition(atom_index, position)
    if mode == "planar":
        conformer = molecule.GetConformer()
        for atom_index in range(molecule.GetNumAtoms()):
            position = conformer.GetAtomPosition(atom_index)
            position.z = 0
            conformer.SetAtomPosition(atom_index, position)
    if mode in {"no_h", "planar"}:
        molecule = Chem.RemoveHs(molecule)
    if mode == "reorder":
        molecule = Chem.RenumberAtoms(molecule, list(reversed(range(molecule.GetNumAtoms()))))
    Chem.MolToPDBFile(molecule, str(path))


def run_checker(pdb_path: Path, report_path: Path, smiles: str = SMILES,
                extra_args: tuple = ()) -> dict:
    completed = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--ligand-pdb", str(pdb_path),
            "--reference-smiles", smiles,
            "--af3-smiles", smiles,
            "--report", str(report_path),
            *extra_args,
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise AssertionError(f"Checker failed:\n{completed.stderr}")
    return json.loads(report_path.read_text(encoding="utf-8"))


class ChiralityCheckerTests(unittest.TestCase):
    def check_case(self, expected, smiles=SMILES, mode="normal", mirror=False, extra_args=()):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            pdb_path = base / "ligand.pdb"
            write_ligand_pdb(pdb_path, mirror=mirror, smiles=smiles, mode=mode)
            report = run_checker(pdb_path, base / "report.json", smiles, extra_args)
            self.assertEqual(report["verdict"], expected, report)
            return report

    def test_unique_mapping_at_limit_one(self):
        report = self.check_case("PASS", extra_args=("--max-mappings", "1"))
        self.assertEqual(report["mapping_count"], 1)

    def test_symmetry_exactly_at_limit(self):
        smiles = "C[C@H](O)c1ccccc1"
        baseline = self.check_case("PASS", smiles=smiles)
        count = baseline["mapping_count"]
        self.assertGreater(count, 1)
        self.check_case("PASS", smiles=smiles, extra_args=("--max-mappings", str(count)))

    def test_symmetry_exceeds_limit(self):
        self.check_case("AMBIGUOUS_MAPPING", smiles="C[C@H](O)c1ccccc1",
                        extra_args=("--max-mappings", "1"))

    def test_missing_hydrogen_coordinates(self):
        self.check_case("PASS", mode="no_h")

    def test_reordered_atoms(self):
        self.check_case("PASS", mode="reorder")

    def test_planar_coordinates(self):
        self.check_case("LOW_GEOMETRY_QUALITY", mode="planar")

    def test_protonation_mismatch(self):
        self.check_case("CHEMICAL_STATE_MISMATCH", smiles="C[C@H](N)C(=O)O",
                        extra_args=("--af3-smiles", "C[C@H](N)C(=O)[O-]"))

    def test_multiple_centres(self):
        self.check_case("PASS", smiles="C[C@H](F)[C@H](O)Cl")

    def test_mirrored_multiple_centres(self):
        self.check_case("CHIRALITY_MISMATCH", smiles="C[C@H](F)[C@H](O)Cl", mirror=True)

    def test_unspecified_stereochemistry(self):
        self.check_case("NO_DEFINED_STEREOCENTERS", smiles="CC(F)Cl")

    def test_meso_polyol_and_its_mirror_match(self):
        smiles = "OC[C@H](O)[C@@H](O)[C@H](O)CO"
        for mirror in (False, True):
            with self.subTest(mirror=mirror):
                report = self.check_case("PASS", smiles=smiles, mirror=mirror)
                self.assertIn("r", report["defined_reference_centres"].values())
                selected = report["mapping_assessments"][report["selected_mapping_index"]]
                self.assertTrue(all(c["status"] == "MATCH" for c in selected["centres"]))

    def test_pseudoasymmetric_inversion_is_detected(self):
        smiles = "OC[C@H](O)[C@@H](O)[C@H](O)CO"
        opposite = "OC[C@H](O)[C@H](O)[C@H](O)CO"
        # Only the middle centre changes: global mirroring is not this case.
        self.check_case("CHIRALITY_MISMATCH", smiles=opposite,
                        extra_args=("--reference-smiles", smiles))

    def test_one_of_two_centres_inverted(self):
        self.check_case("CHIRALITY_MISMATCH", smiles="C[C@H](F)[C@@H](O)Cl",
                        extra_args=("--reference-smiles", "C[C@H](F)[C@H](O)Cl"))

    def test_correct_coordinates_pass(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            pdb_path = base / "correct.pdb"
            report_path = base / "correct.json"
            write_ligand_pdb(pdb_path)
            report = run_checker(pdb_path, report_path)
            self.assertEqual(report["verdict"], "PASS")

    def test_mirrored_coordinates_are_flagged(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            pdb_path = base / "mirrored.pdb"
            report_path = base / "mirrored.json"
            write_ligand_pdb(pdb_path, mirror=True)
            report = run_checker(pdb_path, report_path)
            self.assertEqual(report["verdict"], "CHIRALITY_MISMATCH")


if __name__ == "__main__":
    unittest.main()
