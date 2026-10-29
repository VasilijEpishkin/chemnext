"""Параметризация малых молекул OpenFF и экспорт для GROMACS."""
import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit
from rdkit import Chem


def rdkit_to_openff_mol(rdkit_mol: Chem.Mol, name: str) -> Optional[Molecule]:
    """Преобразует структуру RDKit в OpenFF с исходными конформерами."""
    try:
        off_mol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        off_mol.name = name
        return off_mol
    except Exception as exc:
        logging.error("Ошибка конвертации '%s' в OpenFF: %s", name, exc)
        return None


def verify_coordinates_preserved(
    rdkit_mol: Chem.Mol,
    off_mol: Molecule,
    tolerance: float = 1e-3,
) -> bool:
    """Проверяет координаты в ангстремах (порог — максимальное отличие по оси)."""
    try:
        rdkit_coords = np.asarray(rdkit_mol.GetConformer().GetPositions(), dtype=float)
        off_coords = np.asarray(off_mol.conformers[0].m_as(unit.angstrom), dtype=float)
    except (AttributeError, IndexError, ValueError) as exc:
        logging.error("Нет доступных 3D-координат: %s", exc)
        return False

    if rdkit_coords.shape != off_coords.shape or rdkit_coords.ndim != 2:
        logging.error("Несоответствие размеров координат RDKit и OpenFF")
        return False
    if rdkit_coords.shape[1] != 3 or not np.isfinite(rdkit_coords).all() or not np.isfinite(off_coords).all():
        logging.error("Некорректные или нечисловые координаты молекулы")
        return False

    max_diff = float(np.max(np.abs(rdkit_coords - off_coords)))
    if max_diff > tolerance:
        logging.error("Координаты изменились: max_diff=%.6f Å (порог %.6f Å)", max_diff, tolerance)
        return False
    return True


def create_interchange(
    off_mol: Molecule,
    forcefield: str = "openff-2.2.0.offxml",
) -> Optional[object]:
    """Формирует топологию и параметры OpenFF для одной молекулы."""
    try:
        ff = ForceField(forcefield)
        topology = Topology.from_molecules([off_mol])
        interchange = ff.create_interchange(topology)

        # Большая коробка для одной молекулы; это не настройка MD-эксперимента.
        interchange.box = [10, 10, 10] * unit.nanometer
        return interchange
    except Exception as exc:
        logging.error("Ошибка создания Interchange: %s", exc)
        return None


def validate_gromacs_files(gro_file: Path, top_file: Path, expected_atoms: int) -> Tuple[bool, str]:
    """Базовый контроль полноты экспорта: число атомов и конечность координат."""
    try:
        if not gro_file.is_file() or not top_file.is_file() or top_file.stat().st_size == 0:
            return False, "Отсутствует или пустой файл .gro/.top"

        lines = gro_file.read_text(encoding="utf-8").splitlines()
        if len(lines) < 3:
            return False, "Повреждённый .gro"
        atom_count = int(lines[1].strip())
        if atom_count != expected_atoms or len(lines) != atom_count + 3:
            return False, f"Некорректное число атомов в .gro: {atom_count} вместо {expected_atoms}"

        for line in lines[2:2 + atom_count]:
            # GRO хранит координаты в нм, обычно по три поля шириной 8 символов.
            xyz_nm = [float(line[i:i + 8]) for i in (20, 28, 36)]
            if not np.isfinite(xyz_nm).all():
                return False, "Нечисловые координаты в .gro"

        box = [float(value) for value in lines[-1].split()]
        if len(box) not in (3, 9) or not np.isfinite(box).all() or any(v <= 0 for v in box[:3]):
            return False, "Некорректная периодическая коробка в .gro"
        return True, "Успешно"
    except (OSError, ValueError, IndexError) as exc:
        return False, f"Ошибка проверки GROMACS-файлов: {exc}"


def export_to_gromacs(
    interchange,
    output_dir: Path,
    prefix: str = "system",
) -> Tuple[bool, str]:
    """Экспортирует .gro и .top и проверяет структуру результата."""
    gro_file = output_dir / f"{prefix}.gro"
    top_file = output_dir / f"{prefix}.top"

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        interchange.to_gromacs(prefix=str(output_dir / prefix))
        ok, message = validate_gromacs_files(
            gro_file, top_file, expected_atoms=interchange.topology.n_atoms
        )
        if not ok:
            for path in (gro_file, top_file):
                path.unlink(missing_ok=True)
            return False, message
        logging.info("GROMACS-файлы: %s, %s", gro_file.name, top_file.name)
        return True, "Успешно"
    except Exception as exc:
        logging.error("Ошибка экспорта в GROMACS: %s", exc)
        for path in (gro_file, top_file):
            path.unlink(missing_ok=True)
        return False, str(exc)


def parameterize_molecule(
    rdkit_mol: Chem.Mol,
    name: str,
    output_dir: Path,
    forcefield: str = "openff-2.2.0.offxml",
    check_coords: bool = True,
) -> Tuple[bool, str]:
    """Один лиганд: RDKit → OpenFF → Interchange → GROMACS."""
    off_mol = rdkit_to_openff_mol(rdkit_mol, name)
    if off_mol is None:
        return False, "Ошибка конвертации в OpenFF"

    if check_coords and not verify_coordinates_preserved(rdkit_mol, off_mol):
        # Недопустимо создавать топологию с незаметно изменившейся геометрией.
        return False, "Проверка сохранения координат не пройдена"

    interchange = create_interchange(off_mol, forcefield)
    if interchange is None:
        return False, "Ошибка параметризации OpenFF"
    return export_to_gromacs(interchange, output_dir, "topol")
