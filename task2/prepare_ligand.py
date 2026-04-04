"""Pipeline подготовки лигандов: RDKit → 3D (ETKDGv3) → Meeko → PDBQT."""
import logging
from pathlib import Path
from typing import List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from meeko import MoleculePreparation, PDBQTWriterLegacy
from utils import load_sdf_molecules, load_smi_molecules


def normalize_mol(mol: Chem.Mol, random_seed: int = 42) -> Optional[Chem.Mol]:
    """Добавляет водороды, генерирует 3D-конформер через ETKDGv3
    и оптимизирует геометрию."""
    try:
        mol_h = Chem.AddHs(mol)
        Chem.SanitizeMol(mol_h)

        params = rdDistGeom.ETKDGv3()
        params.randomSeed = random_seed
        params.enforceChirality = True
        params.useMacrocycleTorsions = True

        AllChem.EmbedMolecule(mol_h, params)
        AllChem.MMFFOptimizeMolecule(mol_h)

        return mol_h
    except Exception as e:
        logging.error("Ошибка нормализации молекулы: %s", e)
        return None


def meeko_setup(
    mol: Chem.Mol, flexible_rings: bool = False, min_ring_size: int = 6
) -> Optional[List]:
    """Подготавливает молекулу через Meeko
    и возвращает список MoleculeSetup."""
    try:
        preparator = MoleculePreparation(
            rigid_macrocycles=not flexible_rings,
            min_ring_size=min_ring_size,
        )
        molsetups = preparator.prepare(mol)
        return molsetups
    except Exception as e:
        logging.error("Ошибка подготовки Meeko: %s", e)
        return None


def mol_to_pdbqt(molsetup, output_path: Path) -> bool:
    """Сохраняет MoleculeSetup в PDBQT-файл."""
    try:
        pdbqt_string, success, err_msg = PDBQTWriterLegacy.write_string(molsetup)
        if not success:
            raise RuntimeError(err_msg or "Неизвестная ошибка записи PDBQT")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(pdbqt_string)

        logging.info("Сохранён: %s", output_path.name)
        return True
    except Exception as e:
        logging.error("Ошибка записи PDBQT: %s", e)
        return False


def convert_single_molecule(
    mol: Chem.Mol,
    name: str,
    output_dir: Path,
    flexible_rings: bool = False,
    min_ring_size: int = 6,
) -> Tuple[bool, str]:
    """Выполняет полную конвертацию одной молекулы."""
    normalized = normalize_mol(mol)
    if normalized is None:
        return False, "Ошибка нормализации"

    molsetups = meeko_setup(normalized, flexible_rings, min_ring_size)
    if not molsetups:
        return False, "Ошибка подготовки Meeko"

    output_file = output_dir / f"{name}.pdbqt"
    success = mol_to_pdbqt(molsetups[0], output_file)
    return success, output_file.name if success else "Ошибка записи"


def convert_smi_file(
    smi_file: Path,
    output_dir: Path,
    flexible_rings: bool = False,
    min_ring_size: int = 6,
) -> List[Tuple[str, str, bool]]:
    """Конвертирует все молекулы из SMILES-файла в PDBQT."""

    molecules = load_smi_molecules(smi_file)
    logging.info("Загружено %d молекул из %s", len(molecules), smi_file.name)

    results = []
    for mol, name in molecules:
        ok, msg = convert_single_molecule(mol, name, output_dir,
                                          flexible_rings, min_ring_size)
        results.append((name, msg, ok))
    return results


def convert_sdf_file(
    sdf_file: Path,
    output_dir: Path,
    flexible_rings: bool = False,
    min_ring_size: int = 6,
) -> List[Tuple[str, str, bool]]:
    """Конвертирует все молекулы из SDF-файла в PDBQT."""

    molecules = load_sdf_molecules(sdf_file)
    logging.info("Загружено %d молекул из %s", len(molecules), sdf_file.name)

    results = []
    for mol, name in molecules:
        ok, msg = convert_single_molecule(mol, name, output_dir,
                                          flexible_rings, min_ring_size)
        results.append((name, msg, ok))
    return results
