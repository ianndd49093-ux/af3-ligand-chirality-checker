"""Reproducible synthetic stress benchmark; no AF3 predictions are implied.

Run: python tests/stress_chirality.py --output validation/stress.json
"""
import argparse
import collections
import contextlib
import io
import json
from pathlib import Path
import random
import sys
import tempfile
import time
from unittest.mock import patch

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import af3_ligand_chirality_check as checker

# Explicit structures, not inferred protein sequences or asserted drug identities.
MOLECULES = {
    "halogenated": "C[C@H](F)Cl",
    "amino_acid": "C[C@H](N)C(=O)O",
    "amino_acid_zwitterion": "C[C@H]([NH3+])C(=O)[O-]",
    "hydroxy_acid": "C[C@H](O)C(=O)O",
    "thiol_sidechain": "N[C@@H](CS)C(=O)O",
    "branched_sidechain": "N[C@@H](C(C)C)C(=O)O",
    "aromatic_sidechain": "N[C@@H](Cc1ccccc1)C(=O)O",
    "cyclic_amino_acid": "O=C(O)[C@@H]1CCCN1",
    "aryl_alcohol": "C[C@H](O)c1ccccc1",
    "aryl_hydroxy_acid": "O=C(O)[C@H](O)c1ccccc1",
    "heteroaryl_alcohol": "C[C@H](O)c1ccncc1",
    "flexible_alcohol": "CCCC[C@H](O)CCC",
    "ether": "CO[C@H](C)CC",
    "amide": "CC(=O)N[C@H](C)CC",
    "quaternary_carbon": "C[C@](F)(Cl)Br",
    "two_centres": "C[C@H](F)[C@H](O)Cl",
    "three_centres": "C[C@H](F)[C@H](O)[C@H](Cl)Br",
    "dipeptide": "N[C@@H](C)C(=O)N[C@@H](CO)C(=O)O",
    "substituted_ring": "C[C@H]1CCC[C@H](O)C1",
    "polyol": "OC[C@H](O)[C@@H](O)[C@H](O)CO",
}
SEEDS = (7, 23, 101, 2026)
MODES = ("original", "mirror", "no_h", "mirror_no_h", "reordered", "rigid_motion", "no_conect")


def prepare(smiles, seed):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.enforceChirality = True
    if AllChem.EmbedMolecule(mol, p) != 0:
        raise RuntimeError("fixture embedding failed")
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise RuntimeError("fixture MMFF parameters unavailable")
    if AllChem.MMFFOptimizeMolecule(mol, maxIters=2000) != 0:
        raise RuntimeError("fixture minimization did not converge")
    return mol


def transformed(mol, mode, seed):
    result = Chem.Mol(mol)
    conf = result.GetConformer()
    for i in range(result.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        if mode.startswith("mirror"):
            p.x = -p.x
        if mode == "rigid_motion":
            p.x, p.y, p.z = -p.y + 12, p.x - 8, p.z + 3
        conf.SetAtomPosition(i, p)
    if mode in ("no_h", "mirror_no_h", "no_conect"):
        result = Chem.RemoveHs(result)
    if mode == "reordered":
        order = list(range(result.GetNumAtoms()))
        random.Random(seed).shuffle(order)
        result = Chem.RenumberAtoms(result, order)
    return result


def run_case(mol, smiles, path, mode):
    pdb, report = path / "ligand.pdb", path / "report.json"
    block = Chem.MolToPDBBlock(mol)
    if mode == "no_conect":
        block = "\n".join(line for line in block.splitlines() if not line.startswith("CONECT")) + "\n"
    pdb.write_text(block)
    argv = ["checker", "--pdb", str(pdb), "--reference-smiles", smiles,
            "--af3-smiles", smiles, "--report", str(report)]
    start = time.perf_counter()
    with patch.object(sys, "argv", argv), contextlib.redirect_stderr(io.StringIO()):
        code = checker.main()
    data = json.loads(report.read_text())
    return code, data, time.perf_counter() - start


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows, failures = [], []
    start = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="chirality-stress-") as directory:
        path = Path(directory)
        for name, smiles in MOLECULES.items():
            for seed in SEEDS:
                try:
                    mol = prepare(smiles, seed)
                except Exception as exc:
                    failures.append({"molecule": name, "seed": seed, "error": str(exc)})
                    continue
                for mode in MODES:
                    # A meso molecule can be identical to its mirror image.
                    # Derive that relationship from the reference graph, not
                    # from the checker's verdict or generated coordinates.
                    reference = Chem.MolFromSmiles(smiles)
                    inverted = Chem.Mol(reference)
                    for atom in inverted.GetAtoms():
                        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                            atom.InvertChirality()
                    mirror_equivalent = Chem.MolToSmiles(reference) == Chem.MolToSmiles(inverted)
                    expected = ("CHIRALITY_MISMATCH"
                                if mode.startswith("mirror") and not mirror_equivalent else "PASS")
                    try:
                        with rdBase.BlockLogs():
                            code, result, elapsed = run_case(transformed(mol, mode, seed), smiles, path, mode)
                        row = dict(molecule=name, smiles=smiles, seed=seed, mode=mode,
                                   expected=expected, actual=result["verdict"], exit_code=code,
                                   seconds=round(elapsed, 6))
                        row["ok"] = code == 0 and row["actual"] == expected
                        if not row["ok"]:
                            row["report"] = result
                    except Exception as exc:
                        row = dict(molecule=name, seed=seed, mode=mode, expected=expected,
                                   actual="EXCEPTION", ok=False, error=str(exc))
                    rows.append(row)
            print(f"{name}: {sum(r['ok'] for r in rows if r['molecule'] == name)}/{sum(r['molecule'] == name for r in rows)}", flush=True)
    output = dict(rdkit_version=rdBase.rdkitVersion, molecule_count=len(MOLECULES),
                  seeds=list(SEEDS), modes=list(MODES), planned_cases=len(MOLECULES)*len(SEEDS)*len(MODES),
                  tested=len(rows), passed=sum(r["ok"] for r in rows), fixture_failures=failures,
                  seconds=round(time.perf_counter()-start, 3),
                  verdict_counts=dict(collections.Counter(r["actual"] for r in rows)), cases=rows,
                  scope="Synthetic RDKit-generated ligands; not independent or real AF3 validation.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps({k:v for k,v in output.items() if k != "cases"}, indent=2))
    return int(output["passed"] != output["planned_cases"])


if __name__ == "__main__":
    raise SystemExit(main())
