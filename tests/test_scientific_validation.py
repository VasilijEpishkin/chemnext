"""Химические smoke-тесты подготовки структур.

Запуск: conda activate chemnext && python -m pytest -q
Без OpenFF или Meeko соответствующие тесты пропускаются, а не имитируют расчёт.
"""
import importlib.util
import math
from pathlib import Path
import sys

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]


def ethanol_3d():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) >= 0
    assert AllChem.MMFFOptimizeMolecule(mol) == 0
    return mol


def load_gromacs_core():
    pytest.importorskip("openff.toolkit")
    spec = importlib.util.spec_from_file_location("chemnext_gromacs_core", ROOT / "task1/core.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_docking_preparation(monkeypatch):
    pytest.importorskip("meeko")
    # Исходный CLI импортирует utils из собственного каталога.
    monkeypatch.syspath_prepend(str(ROOT / "task2"))
    sys.modules.pop("utils", None)
    spec = importlib.util.spec_from_file_location("chemnext_prepare_ligand", ROOT / "task2/prepare_ligand.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_gromacs_export_for_known_ethanol(tmp_path):
    core = load_gromacs_core()
    mol = ethanol_3d()
    ok, message = core.parameterize_molecule(mol, "ethanol", tmp_path)
    assert ok, message

    gro = tmp_path / "topol.gro"
    top = tmp_path / "topol.top"
    assert gro.is_file() and top.is_file()
    assert core.validate_gromacs_files(gro, top, mol.GetNumAtoms())[0]
    off = core.rdkit_to_openff_mol(mol, "ethanol")
    assert core.verify_coordinates_preserved(mol, off)


def test_coordinate_mismatch_prevents_parameterization(monkeypatch, tmp_path):
    core = load_gromacs_core()
    monkeypatch.setattr(core, "verify_coordinates_preserved", lambda *args: False)
    monkeypatch.setattr(core, "create_interchange", lambda *args: pytest.fail("Should stop before FF"))
    ok, message = core.parameterize_molecule(ethanol_3d(), "ethanol", tmp_path)
    assert not ok
    assert "координат" in message.lower()
    assert not (tmp_path / "topol.gro").exists()


def test_export_check_rejects_invalid_gro(tmp_path):
    core = load_gromacs_core()
    (tmp_path / "test.top").write_text("[ molecules ]\nETH 1\n")
    (tmp_path / "test.gro").write_text("ethanol\n2\ninvalid\n10 10 10\n")
    ok, _ = core.validate_gromacs_files(tmp_path / "test.gro", tmp_path / "test.top", 9)
    assert not ok


@pytest.mark.parametrize("smiles", ["CCO", "c1ccccc1"])
def test_docking_known_molecules_have_finite_mmff_energy_and_geometry(monkeypatch, smiles):
    prep = load_docking_preparation(monkeypatch)
    mol = Chem.MolFromSmiles(smiles)
    out = prep.normalize_mol(mol)
    assert out is not None
    conf = out.GetConformer()
    for bond in out.GetBonds():
        length = conf.GetAtomPosition(bond.GetBeginAtomIdx()).Distance(
            conf.GetAtomPosition(bond.GetEndAtomIdx())
        )
        assert 0.8 < length < 2.0, f"Неправдоподобная длина связи: {length:.3f} Å"
    props = AllChem.MMFFGetMoleculeProperties(out)
    ff = AllChem.MMFFGetMoleculeForceField(out, props)
    assert math.isfinite(ff.CalcEnergy())


def test_docking_pdbqt_for_ethanol(monkeypatch, tmp_path):
    prep = load_docking_preparation(monkeypatch)
    ok, message = prep.convert_single_molecule(Chem.MolFromSmiles("CCO"), "ethanol", tmp_path)
    assert ok, message
    pdbqt = (tmp_path / "ethanol.pdbqt").read_text()
    assert prep.validate_pdbqt_text(pdbqt) == (True, "Успешно")
    atoms = [line for line in pdbqt.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    assert len(atoms) >= 3  # C, C, O (неполярные H Meeko может объединять)


def test_docking_rejects_embedding_failure(monkeypatch):
    prep = load_docking_preparation(monkeypatch)
    monkeypatch.setattr(AllChem, "EmbedMolecule", lambda *args, **kwargs: -1)
    assert prep.normalize_mol(Chem.MolFromSmiles("CCO")) is None


def test_docking_rejects_mmff_nonconvergence(monkeypatch):
    prep = load_docking_preparation(monkeypatch)
    monkeypatch.setattr(AllChem, "MMFFOptimizeMolecule", lambda *args, **kwargs: 1)
    assert prep.normalize_mol(Chem.MolFromSmiles("CCO")) is None


def test_docking_rejects_missing_mmff_params(monkeypatch):
    prep = load_docking_preparation(monkeypatch)
    monkeypatch.setattr(AllChem, "MMFFHasAllMoleculeParams", lambda mol: False)
    assert prep.normalize_mol(Chem.MolFromSmiles("CCO")) is None


def test_docking_rejects_empty_pdbqt(monkeypatch):
    prep = load_docking_preparation(monkeypatch)
    assert prep.validate_pdbqt_text("REMARK no atom records\n")[0] is False
