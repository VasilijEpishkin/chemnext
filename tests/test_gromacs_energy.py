"""Необязательный сквозной sanity-check GROMACS на этаноле.

Нужны GROMACS (gmx в PATH), OpenFF и AmberTools из environment.yml.
Результат проверяет расчёт энергии, но не валидирует физическую модель MD.
"""
import math
from pathlib import Path
import shutil
import subprocess

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from test_scientific_validation import load_gromacs_core

ROOT = Path(__file__).resolve().parents[1]


def _run(cmd, cwd, input_text=None):
    result = subprocess.run(
        cmd, input=input_text, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        cwd=cwd, check=False, timeout=120,
    )
    assert result.returncode == 0, result.stdout[-3000:]


def test_gromacs_can_minimize_ethanol_and_report_finite_energy(tmp_path):
    if not shutil.which("gmx"):
        pytest.skip("GROMACS (gmx) не установлен")
    core = load_gromacs_core()

    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) >= 0
    assert AllChem.MMFFOptimizeMolecule(mol) == 0
    ok, msg = core.parameterize_molecule(mol, "ethanol", tmp_path)
    assert ok, msg

    tpr = tmp_path / "minimization.tpr"
    _run([
        "gmx", "grompp", "-f", str(ROOT / "validation/minimization.mdp"),
        "-c", "topol.gro", "-p", "topol.top", "-o", str(tpr)
    ], tmp_path)
    _run(["gmx", "mdrun", "-s", str(tpr), "-deffnm", "minimization", "-nt", "1"], tmp_path)

    xvg = tmp_path / "potential.xvg"
    _run(["gmx", "energy", "-f", "minimization.edr", "-o", str(xvg)], tmp_path, "Potential\n0\n")
    samples = [
        float(line.split()[1])
        for line in xvg.read_text().splitlines()
        if line.strip() and not line.startswith(("#", "@"))
    ]
    assert samples and all(math.isfinite(energy) for energy in samples)
