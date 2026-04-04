"""Core-логика: параметризация OpenFF и экспорт в Gromacs."""
import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit
from rdkit import Chem


def rdkit_to_openff_mol(rdkit_mol: Chem.Mol, name: str) -> Optional[Molecule]:
    """Конвертирует RDKit молекулу в OpenFF Molecule,
    сохраняя 3D-координаты."""
    try:
        off_mol = Molecule.from_rdkit(
            rdkit_mol,
            allow_undefined_stereo=True,
        )
        off_mol.name = name
        return off_mol
    except Exception as e:
        logging.error("Ошибка конвертации '%s' в OpenFF: %s", name, e)
        return None


def verify_coordinates_preserved(
    rdkit_mol: Chem.Mol,
    off_mol: Molecule,
    tolerance: float = 1e-3,
) -> bool:
    """Проверяет, что координаты RDKit и OpenFF молекул совпадают."""
    rdkit_coords = rdkit_mol.GetConformer().GetPositions()
    off_coords = off_mol.conformers[0].m_as(unit.angstrom)

    if len(rdkit_coords) != len(off_coords):
        logging.warning(
            "Несоответствие числа атомов: RDKit=%d, OpenFF=%d",
            len(rdkit_coords), len(off_coords),
        )
        return False

    max_diff = np.max(np.abs(rdkit_coords - off_coords))
    if max_diff > tolerance:
        logging.warning(
            "Координаты изменились (max_diff=%.6f A)", max_diff
        )
        return False

    return True


def create_interchange(
    off_mol: Molecule,
    forcefield: str = "openff-2.2.0.offxml",
) -> Optional[object]:
    """Создаёт Interchange объект для одной молекулы."""
    try:
        ff = ForceField(forcefield)
        topology = Topology.from_molecules([off_mol])
        interchange = ff.create_interchange(topology)
        return interchange
    except Exception as e:
        logging.error("Ошибка создания Interchange: %s", e)
        return None


def export_to_gromacs(
    interchange,
    output_dir: Path,
    prefix: str = "system",
) -> Tuple[bool, str]:
    """Экспортирует Interchange в файлы Gromacs (.gro + .top).

    Returns:
        (success, message)
    """
    try:
        output_dir.mkdir(parents=True, exist_ok=True)

        interchange.to_gromacs(
            prefix=str(output_dir / prefix),
        )

        gro_file = output_dir / f"{prefix}.gro"
        top_file = output_dir / f"{prefix}.top"

        if not gro_file.exists():
            return False, "Файл .gro не создан"
        if not top_file.exists():
            return False, "Файл .top не создан"

        logging.info(
            "Gromacs файлы: %s, %s", gro_file.name, top_file.name
        )
        return True, "Успешно"

    except Exception as e:
        logging.error("Ошибка экспорта в Gromacs: %s", e)
        return False, str(e)


def parameterize_molecule(
    rdkit_mol: Chem.Mol,
    name: str,
    output_dir: Path,
    forcefield: str = "openff-2.0.0.offxml",
    check_coords: bool = True,
) -> Tuple[bool, str]:
    """Полный pipeline параметризации одной молекулы.

    1. Конвертация RDKit → OpenFF (с сохранением координат)
    2. Создание Interchange
    3. Экспорт в Gromacs (.gro + topol.top)
    """
    off_mol = rdkit_to_openff_mol(rdkit_mol, name)
    if off_mol is None:
        return False, "Ошибка конвертации в OpenFF"

    if check_coords:
        if not verify_coordinates_preserved(rdkit_mol, off_mol):
            logging.warning(
                "Координаты молекулы '%s' могли измениться", name
            )

    interchange = create_interchange(off_mol, forcefield)
    if interchange is None:
        return False, "Ошибка создания Interchange"

    success, msg = export_to_gromacs(interchange, output_dir, "topol")
    return success, msg
