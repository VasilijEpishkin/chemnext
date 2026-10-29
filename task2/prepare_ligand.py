"""Подготовка лигандов для AutoDock Vina: RDKit → Meeko → PDBQT."""
import logging
import math
from pathlib import Path
from typing import List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from meeko import MoleculePreparation, PDBQTWriterLegacy
from utils import load_sdf_molecules, load_smi_molecules


def normalize_mol(mol: Chem.Mol, random_seed: int = 42) -> Optional[Chem.Mol]:
    """Добавляет H, генерирует ETKDGv3-конформер и оптимизирует MMFF94.

    Возвращает None, если молекулу нельзя корректно подготовить этим методом.
    """
    try:
        if mol is None:
            raise ValueError("Пустая молекула")
        mol_h = Chem.AddHs(Chem.Mol(mol))
        Chem.SanitizeMol(mol_h)
        if not AllChem.MMFFHasAllMoleculeParams(mol_h):
            raise ValueError("Для молекулы отсутствуют параметры MMFF94")

        params = rdDistGeom.ETKDGv3()
        params.randomSeed = random_seed
        params.enforceChirality = True
        params.useMacrocycleTorsions = True

        conformer_id = AllChem.EmbedMolecule(mol_h, params)
        if conformer_id < 0:
            raise ValueError("ETKDGv3 не смог построить 3D-конформер")

        status = AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
        if status != 0:
            reason = "нет сходимости после 500 итераций" if status == 1 else "сбой MMFF94"
            raise ValueError(f"Оптимизация геометрии: {reason} (код {status})")

        conformer = mol_h.GetConformer()
        for atom_index in range(mol_h.GetNumAtoms()):
            pos = conformer.GetAtomPosition(atom_index)
            if not all(math.isfinite(value) for value in (pos.x, pos.y, pos.z)):
                raise ValueError("Нечисловые координаты после оптимизации")

        properties = AllChem.MMFFGetMoleculeProperties(mol_h)
        forcefield = AllChem.MMFFGetMoleculeForceField(mol_h, properties) if properties else None
        if forcefield is None:
            raise ValueError("Не удалось вычислить энергию MMFF94")
        energy = float(forcefield.CalcEnergy())
        if not math.isfinite(energy):
            raise ValueError("Неконечная энергия MMFF94")
        logging.debug("Энергия подготовленного конформера MMFF94: %.4f ккал/моль", energy)
        return mol_h
    except Exception as exc:
        logging.error("Ошибка подготовки молекулы: %s", exc)
        return None


def meeko_setup(
    mol: Chem.Mol,
    flexible_rings: bool = False,
    min_ring_size: int = 6,
) -> Optional[List]:
    """Подготавливает атомные типы, заряды и торсии для PDBQT."""
    try:
        preparator = MoleculePreparation(
            rigid_macrocycles=not flexible_rings,
            min_ring_size=min_ring_size,
        )
        return preparator.prepare(mol)
    except Exception as exc:
        logging.error("Ошибка подготовки Meeko: %s", exc)
        return None


def validate_pdbqt_text(pdbqt_string: str) -> Tuple[bool, str]:
    """Проверяет наличие атомов и конечность координат в PDBQT."""
    atoms = [
        line for line in pdbqt_string.splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    if not atoms:
        return False, "PDBQT не содержит атомных записей"
    try:
        for line in atoms:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if not all(math.isfinite(value) for value in xyz):
                return False, "PDBQT содержит нечисловые координаты"
    except ValueError as exc:
        return False, f"Не удалось прочитать координаты PDBQT: {exc}"
    return True, "Успешно"


def mol_to_pdbqt(molsetup, output_path: Path) -> bool:
    """Сохраняет проверенный PDBQT; некорректный файл не создаётся."""
    try:
        pdbqt_string, success, err_msg = PDBQTWriterLegacy.write_string(molsetup)
        if not success:
            raise ValueError(err_msg or "Ошибка Meeko при записи PDBQT")

        valid, message = validate_pdbqt_text(pdbqt_string)
        if not valid:
            raise ValueError(message)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(pdbqt_string, encoding="utf-8")
        logging.info("Сохранён: %s", output_path.name)
        return True
    except Exception as exc:
        logging.error("Ошибка записи PDBQT: %s", exc)
        output_path.unlink(missing_ok=True)
        return False


def convert_single_molecule(
    mol: Chem.Mol,
    name: str,
    output_dir: Path,
    flexible_rings: bool = False,
    min_ring_size: int = 6,
) -> Tuple[bool, str]:
    """Готовит и конвертирует одну молекулу."""
    normalized = normalize_mol(mol)
    if normalized is None:
        return False, "Ошибка генерации или оптимизации 3D-конформера"

    molsetups = meeko_setup(normalized, flexible_rings, min_ring_size)
    if not molsetups:
        return False, "Ошибка параметризации Meeko"

    output_file = output_dir / f"{name}.pdbqt"
    success = mol_to_pdbqt(molsetups[0], output_file)
    return success, output_file.name if success else "Ошибка записи PDBQT"


def convert_smi_file(
    smi_file: Path,
    output_dir: Path,
    flexible_rings: bool = False,
    min_ring_size: int = 6,
) -> List[Tuple[str, str, bool]]:
    """Конвертирует все SMILES из входного файла."""
    molecules = load_smi_molecules(smi_file)
    logging.info("Загружено %d молекул из %s", len(molecules), smi_file.name)
    results = []
    for mol, name in molecules:
        ok, msg = convert_single_molecule(mol, name, output_dir, flexible_rings, min_ring_size)
        results.append((name, msg, ok))
    return results


def convert_sdf_file(
    sdf_file: Path,
    output_dir: Path,
    flexible_rings: bool = False,
    min_ring_size: int = 6,
) -> List[Tuple[str, str, bool]]:
    """Конвертирует все молекулы из SDF."""
    molecules = load_sdf_molecules(sdf_file)
    logging.info("Загружено %d молекул из %s", len(molecules), sdf_file.name)
    results = []
    for mol, name in molecules:
        ok, msg = convert_single_molecule(mol, name, output_dir, flexible_rings, min_ring_size)
        results.append((name, msg, ok))
    return results
