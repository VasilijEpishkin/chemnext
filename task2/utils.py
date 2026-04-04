"""Вспомогательные функции: загрузка молекул и настройка логирования."""
import logging
from pathlib import Path
from typing import List, Tuple
from rdkit import Chem


def setup_logging(log_file: Path, level: int = logging.INFO) -> None:
    """Настраивает логирование в файл и консоль."""
    log_file.parent.mkdir(parents=True, exist_ok=True)

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(formatter)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)

    logging.basicConfig(level=level, handlers=[file_handler, console_handler])


def load_smi_molecules(smi_file: Path) -> List[Tuple[Chem.Mol, str]]:
    """Загружает молекулы из SMILES-файла (SMILES + имя)."""
    molecules: List[Tuple[Chem.Mol, str]] = []

    with open(smi_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    start_idx = (
        1 if lines and lines[0].strip().lower().startswith("smiles") else 0
    )

    for i, line in enumerate(lines[start_idx:], start=start_idx):
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        smiles = parts[0]
        name = parts[1] if len(parts) > 1 else f"mol_{i}"

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            molecules.append((mol, name))
        else:
            logging.warning("Не удалось распарсить SMILES: %s", smiles[:80])

    return molecules


def load_sdf_molecules(sdf_file: Path) -> List[Tuple[Chem.Mol, str]]:
    """Загружает молекулы из SDF-файла."""
    molecules: List[Tuple[Chem.Mol, str]] = []
    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)

    for i, mol in enumerate(supplier):
        if mol is None:
            logging.warning("Не удалось считать молекулу #%d из %s", i, sdf_file.name)
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i}"
        molecules.append((mol, name))

    return molecules
