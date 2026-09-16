#!/usr/bin/env python3
"""Check whether ligand chirality in AF3 coordinates agrees with a reference SMILES.

This is a validation tool: it never changes, minimises, embeds, or docks the
input structure. The ligand PDB must contain the ligand only (not the protein).
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, rdCIPLabeler

LOG = logging.getLogger("af3_chirality_check")


def defined_cip_labels(mol: Chem.Mol) -> Dict[int, str]:
    """Assign accurate CIP labels, including pseudoasymmetric r/s centres."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            atom.ClearProp("_CIPCode")
    rdCIPLabeler.AssignCIPLabels(mol, maxRecursiveIterations=1000000)
    return {atom.GetIdx(): atom.GetProp("_CIPCode") for atom in mol.GetAtoms()
            if atom.HasProp("_CIPCode") and atom.GetProp("_CIPCode") in {"R", "S", "r", "s"}}


def nonstereo_canonical_smiles(smiles: str) -> str:
    """Canonicalise a SMILES while deliberately ignoring stereochemistry."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES could not be parsed or sanitized")
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def vector(a: Any, b: Any) -> Tuple[float, float, float]:
    return (b.x - a.x, b.y - a.y, b.z - a.z)


def norm(v: Sequence[float]) -> float:
    return math.sqrt(sum(x * x for x in v))


def distance(a: Any, b: Any) -> float:
    return norm(vector(a, b))


def determinant(a: Sequence[float], b: Sequence[float], c: Sequence[float]) -> float:
    """Scalar triple product of three vectors."""
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def local_geometry_check(
    mol: Chem.Mol,
    center_idx: int,
    min_bond_ratio: float,
    max_bond_ratio: float,
    clash_ratio: float,
    min_tetrahedral_volume: float,
) -> Dict[str, Any]:
    """Apply conservative, local coordinate-quality checks at one stereocentre.

    The normalised scalar triple product is scale-independent. It remains valid
    when the fourth substituent is an implicit hydrogen, so made-up hydrogen
    coordinates are never introduced as AF3 evidence.
    """
    conf = mol.GetConformer()
    periodic_table = Chem.GetPeriodicTable()
    center = mol.GetAtomWithIdx(center_idx)
    center_pos = conf.GetAtomPosition(center_idx)
    neighbours = [atom for atom in center.GetNeighbors() if atom.GetAtomicNum() > 1]
    issues: List[str] = []
    bond_ratios: List[Dict[str, float]] = []

    if len(neighbours) < 3:
        issues.append("fewer than three heavy-atom substituents")

    unit_vectors: List[Tuple[float, float, float]] = []
    for atom in neighbours:
        atom_pos = conf.GetAtomPosition(atom.GetIdx())
        observed = distance(center_pos, atom_pos)
        expected = (periodic_table.GetRcovalent(center.GetAtomicNum())
                    + periodic_table.GetRcovalent(atom.GetAtomicNum()))
        ratio = observed / expected if expected else float("inf")
        bond_ratios.append({"atom_index": atom.GetIdx(), "ratio": round(ratio, 4)})
        if not min_bond_ratio <= ratio <= max_bond_ratio:
            issues.append(
                f"implausible bond length to atom {atom.GetIdx()} (ratio {ratio:.2f})"
            )
        if observed > 1e-8:
            v = vector(center_pos, atom_pos)
            unit_vectors.append((v[0] / observed, v[1] / observed, v[2] / observed))

    for first, second in itertools.combinations(neighbours, 2):
        # Directly bonded neighbour pairs are expected in small rings and are
        # not treated as a non-bonded steric clash.
        if mol.GetBondBetweenAtoms(first.GetIdx(), second.GetIdx()) is not None:
            continue
        first_pos = conf.GetAtomPosition(first.GetIdx())
        second_pos = conf.GetAtomPosition(second.GetIdx())
        observed = distance(first_pos, second_pos)
        expected = (periodic_table.GetRcovalent(first.GetAtomicNum())
                    + periodic_table.GetRcovalent(second.GetAtomicNum()))
        if expected and observed / expected < clash_ratio:
            issues.append(
                f"neighbour clash between atoms {first.GetIdx()} and {second.GetIdx()}"
            )

    volumes = [abs(determinant(*triple))
               for triple in itertools.combinations(unit_vectors, 3)]
    min_volume = min(volumes) if volumes else 0.0
    if min_volume < min_tetrahedral_volume:
        issues.append(
            f"near-planar or degenerate tetrahedron (normalised volume {min_volume:.3f})"
        )

    return {
        "ok": not issues,
        "issues": issues,
        "bond_ratios": bond_ratios,
        "min_normalised_tetrahedral_volume": round(min_volume, 4),
    }


def read_and_map_ligand(template: Chem.Mol, ligand_pdb: Path) -> Chem.Mol:
    """Read ligand coordinates and assign SMILES bond orders to their graph."""
    coordinate_mol = Chem.MolFromPDBFile(str(ligand_pdb), sanitize=False, removeHs=False)
    if coordinate_mol is None:
        raise ValueError("RDKit could not parse the ligand PDB")
    coordinate_heavy = Chem.RemoveHs(coordinate_mol, sanitize=False)
    if coordinate_heavy.GetNumAtoms() != template.GetNumAtoms():
        raise ValueError(
            "heavy-atom count differs between PDB coordinates and reference SMILES"
        )
    try:
        mapped = AllChem.AssignBondOrdersFromTemplate(template, coordinate_heavy)
        Chem.SanitizeMol(mapped)
    except Exception as exc:
        raise ValueError(f"SMILES-to-coordinate bond-order mapping failed: {exc}") from exc
    return mapped


def mapping_assessment(
    mapping: Sequence[int],
    reference_centers: Dict[int, str],
    observed_centers: Dict[int, str],
    mapped_mol: Chem.Mol,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """Compare all reference chiral centres for one possible atom mapping."""
    checks: List[Dict[str, Any]] = []
    statuses: List[str] = []
    for reference_idx, expected_cip in sorted(reference_centers.items()):
        target_idx = mapping[reference_idx]
        geometry = local_geometry_check(
            mapped_mol, target_idx, args.min_bond_ratio, args.max_bond_ratio,
            args.clash_ratio, args.min_tetrahedral_volume,
        )
        observed_cip = observed_centers.get(target_idx)
        if not geometry["ok"] or observed_cip not in {"R", "S", "r", "s"}:
            centre_status = "LOW_GEOMETRY_QUALITY"
        elif observed_cip == expected_cip:
            centre_status = "MATCH"
        else:
            centre_status = "MISMATCH"
        statuses.append(centre_status)
        checks.append({
            "reference_atom_index": reference_idx,
            "structure_atom_index": target_idx,
            "expected_cip": expected_cip,
            "observed_cip": observed_cip or "unassigned",
            "status": centre_status,
            "geometry": geometry,
        })

    if "LOW_GEOMETRY_QUALITY" in statuses:
        verdict = "LOW_GEOMETRY_QUALITY"
    elif "MISMATCH" in statuses:
        verdict = "CHIRALITY_MISMATCH"
    else:
        verdict = "PASS"
    return {"verdict": verdict, "centres": checks}


def write_report(report_path: Path, report: Dict[str, Any]) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligand-pdb", "--pdb", dest="ligand_pdb", required=True,
                        type=Path, help="Ligand-only PDB extracted from an AF3 structure")
    parser.add_argument("--reference-smiles", required=True,
                        help="Known authoritative isomeric SMILES")
    parser.add_argument("--af3-smiles",
                        help="Optional SMILES supplied to AF3; enables chemical-state checking")
    parser.add_argument("--report", required=True, type=Path, help="Output JSON report")
    parser.add_argument("--max-mappings", type=int, default=100,
                        help="Maximum symmetry-equivalent matches to inspect (default: 100)")
    parser.add_argument("--min-bond-ratio", type=float, default=0.60,
                        help="Minimum observed/covalent-radius-sum bond ratio")
    parser.add_argument("--max-bond-ratio", type=float, default=1.50,
                        help="Maximum observed/covalent-radius-sum bond ratio")
    parser.add_argument("--clash-ratio", type=float, default=0.65,
                        help="Non-bonded neighbour clash ratio")
    parser.add_argument("--min-tetrahedral-volume", type=float, default=0.10,
                        help="Minimum normalised scalar-triple-product magnitude")
    args = parser.parse_args()

    if not args.ligand_pdb.is_file():
        parser.error(f"ligand PDB does not exist: {args.ligand_pdb}")
    if args.max_mappings < 1:
        parser.error("--max-mappings must be at least one")
    if not 0 < args.min_bond_ratio < args.max_bond_ratio:
        parser.error("bond-ratio limits must be positive and ordered")

    report: Dict[str, Any] = {
        "tool": "af3_ligand_chirality_check",
        "ligand_pdb": str(args.ligand_pdb),
        "reference_smiles": args.reference_smiles,
        "thresholds": {
            "min_bond_ratio": args.min_bond_ratio,
            "max_bond_ratio": args.max_bond_ratio,
            "clash_ratio": args.clash_ratio,
            "min_tetrahedral_volume": args.min_tetrahedral_volume,
        },
    }

    try:
        template = Chem.MolFromSmiles(args.reference_smiles)
        if template is None:
            raise ValueError("reference SMILES could not be parsed or sanitized")
        template = Chem.RemoveHs(template)
        reference_centers = defined_cip_labels(Chem.Mol(template))
        report["defined_reference_centres"] = reference_centers

        if args.af3_smiles:
            reference_graph = nonstereo_canonical_smiles(args.reference_smiles)
            af3_graph = nonstereo_canonical_smiles(args.af3_smiles)
            report["chemical_state_check"] = {
                "reference_nonstereo_smiles": reference_graph,
                "af3_nonstereo_smiles": af3_graph,
                "compatible": reference_graph == af3_graph,
            }
            if reference_graph != af3_graph:
                report["verdict"] = "CHEMICAL_STATE_MISMATCH"
                report["message"] = "Reference and AF3 SMILES differ before stereochemistry is considered."
                write_report(args.report, report)
                return 0
        else:
            report["chemical_state_check"] = {
                "compatible": None,
                "note": "Not independently verifiable from a coordinate-only PDB; provide --af3-smiles when available.",
            }

        if not reference_centers:
            report["verdict"] = "NO_DEFINED_STEREOCENTERS"
            report["message"] = "The reference SMILES has no assigned tetrahedral R/S centres to check."
            write_report(args.report, report)
            return 0

        mapped_mol = read_and_map_ligand(template, args.ligand_pdb)
        # A match tuple is indexed by template atom index and stores the
        # corresponding coordinate-molecule atom index.
        matches = mapped_mol.GetSubstructMatches(
            # Request one extra match to distinguish a complete enumeration
            # exactly at the limit from a genuinely truncated search.
            template, uniquify=False, useChirality=False, maxMatches=args.max_mappings + 1
        )
        if not matches:
            report["verdict"] = "CHEMICAL_STATE_MISMATCH"
            report["message"] = "No exact SMILES-to-coordinate graph mapping was found."
            write_report(args.report, report)
            return 0
        if len(matches) > args.max_mappings:
            report["verdict"] = "AMBIGUOUS_MAPPING"
            report["message"] = "Mapping limit exceeded; increase --max-mappings only if needed."
            report["mapping_count_returned"] = len(matches)
            write_report(args.report, report)
            return 0

        observed_mol = Chem.Mol(mapped_mol)
        Chem.AssignStereochemistryFrom3D(observed_mol, replaceExistingTags=True)
        observed_centers = defined_cip_labels(observed_mol)
        assessments = [
            mapping_assessment(match, reference_centers, observed_centers, mapped_mol, args)
            for match in matches
        ]
        verdicts = {assessment["verdict"] for assessment in assessments}
        report["mapping_count"] = len(matches)
        report["mapping_assessments"] = assessments
        # Molecular identity requires one complete stereo-consistent mapping,
        # not agreement across every graph automorphism. Never combine matches
        # from different mappings to fabricate a whole-molecule match.
        if "PASS" in verdicts:
            report["verdict"] = "PASS"
            report["selected_mapping_index"] = next(
                i for i, assessment in enumerate(assessments) if assessment["verdict"] == "PASS")
            report["message"] = "A complete mapping matches every defined stereocentre with reliable geometry."
        elif len(verdicts) == 1:
            report["verdict"] = verdicts.pop()
        else:
            report["verdict"] = "AMBIGUOUS_MAPPING"
            report["message"] = "Valid atom mappings produced different chirality outcomes."
        write_report(args.report, report)
        LOG.info("%s: %s", report["verdict"], args.report)
        return 0
    except (OSError, ValueError) as exc:
        report["verdict"] = "INPUT_ERROR"
        report["message"] = str(exc)
        write_report(args.report, report)
        LOG.error("INPUT_ERROR: %s", exc)
        return 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    raise SystemExit(main())
