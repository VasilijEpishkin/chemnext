# ChemNEXT — Тестовое задание

## Задание 1. Параметризация малых молекул полем OpenFF для Gromacs

Скрипт подготавливает лиганды из SDF-файла для молекулярной динамики в Gromacs с использованием force field OpenFF Sage. Оригинальные 3D-координаты из SDF сохраняются без изменений.

**Запуск:**

```bash
python task1/convert_to_gromacs.py --input input/example.sdf --output output/task1
```

**Параметры:**

| Флаг | По умолчанию | Описание |
|------|-------------|----------|
| `--input`, `-i` | — | SDF-файл с лигандами |
| `--output`, `-o` | `output/ligands` | Директория результатов |
| `--forcefield`, `-ff` | `openff-2.2.0.offxml` | Force field OpenFF |
| `--no-coord-check` | `False` | Отключить проверку координат |
| `--log-file` | `logs/parametrization_log.txt` | Путь к лог-файлу |
| `--verbose`, `-v` | `False` | Подробный DEBUG-вывод |
| `--limit`, `-n` | `None` | Ограничить количество молекул |

**Выходные данные:** для каждой молекулы создаётся отдельная папка с именем молекулы (`_Name` из SDF), содержащая:
- `topol.gro` — структура
- `topol.top` — топология

**Pipeline:**
1. Загрузка SDF через RDKit (с сохранением координат и водородов)
2. Конвертация в OpenFF Molecule (`Molecule.from_rdkit`, `allow_undefined_stereo=True`)
3. Проверка сохранения координат (сравнение RDKit/OpenFF, допуск 0.001 Å)
4. Создание Interchange (`ForceField.create_interchange`)
5. Экспорт в Gromacs (`interchange.to_gromacs`)

---

## Задание 2. Подготовка лигандов для докинга в AutoDock Vina

Скрипт конвертирует лиганды из SMILES/SDF в формат PDBQT. Поддерживает генерацию 3D-конформаций (ETKDGv3), добавление водородов, оптимизацию геометрии (MMFF) и гибкие алифатические макроциклы.

**Запуск:**

```bash
python task2/convert_to_pdbqt.py --input input/example.sdf --output output/task2
python task2/convert_to_pdbqt.py --input input/example.smi --output output/task2
python task2/convert_to_pdbqt.py --input input/ --output output/task2
```

**Параметры:**

| Флаг | По умолчанию | Описание |
|------|-------------|----------|
| `--input`, `-i` | — | `.smi`/`.sdf` файл или директория |
| `--output`, `-o` | `output/ligands` | Директория результатов |
| `--flex-rings`, `-f` | `False` | Гибкие алифатические циклы |
| `--min-ring-size` | `6` | Минимальный размер цикла для гибкости |
| `--log-file` | `logs/conversion_log.txt` | Путь к лог-файлу |
| `--verbose`, `-v` | `False` | Подробный DEBUG-вывод |

**Выходные данные:** по одному `.pdbqt` файлу на каждый лиганд.

**Pipeline:**
1. Загрузка SMILES/SDF через RDKit
2. Добавление водородов (`Chem.AddHs`)
3. Генерация 3D-конформера (ETKDGv3) + оптимизация (MMFF)
4. Подготовка через Meeko (атомные типы AD4, заряды Gasteiger)
5. Запись PDBQT

---

## Установка

### Основной способ (conda)

```bash
conda env create -f environment.yml
conda activate chemnext
```

### Альтернативный способ (pip)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Структура проекта

```
chemnext/
├── environment.yml            # Conda-окружение
├── requirements.txt           # Pip-зависимости (альтернатива)
├── .gitignore
├── input/                     # Тестовые данные
│   ├── example.sdf            
│   └── example.smi            
├── task1/
│   ├── convert_to_gromacs.py  # CLI: параметризация OpenFF → Gromacs
│   ├── core.py                # Pipeline: RDKit → OpenFF → Interchange
│   └── utils.py               # Загрузка SDF, логирование
├── task2/
│   ├── convert_to_pdbqt.py    # CLI: конвертация в PDBQT
│   ├── prepare_ligand.py      # Pipeline: RDKit → ETKDGv3 → Meeko
│   └── utils.py               # Загрузка SDF/SMI, логирование
└── output/                    # Результаты
```
