"""CLI-скрипт для параметризации малых молекул полем OpenFF для Gromacs."""
import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Tuple
from rdkit import Chem
from utils import setup_logging, load_sdf_molecules
from core import parameterize_molecule


def parse_args() -> argparse.Namespace:
    """Парсит аргументы командной строки."""
    parser = argparse.ArgumentParser(
        description="Параметризация лигандов полем OpenFF для Gromacs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Путь к SDF-файлу с лигандами",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("output/ligands"),
        help="Директория для сохранения результатов",
    )
    parser.add_argument(
        "--forcefield",
        "-ff",
        type=str,
        default="openff-2.2.0.offxml",
        help="Force field OpenFF (по умолчанию: openff-2.2.0.offxml)",
    )
    parser.add_argument(
        "--no-coord-check",
        action="store_true",
        help="Отключить проверку сохранения координат",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/parametrization_log.txt"),
        help="Путь к лог-файлу",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Включить подробный вывод (DEBUG)",
    )
    parser.add_argument(
        "--limit",
        "-n",
        type=int,
        default=None,
        help="Ограничить количество обрабатываемых молекул",
    )

    return parser.parse_args()


def main() -> None:
    """Основная функция CLI."""
    args = parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    setup_logging(args.log_file, log_level)

    logging.info("=" * 60)
    logging.info("НАЧАЛО ПАРАМЕТРИЗАЦИИ OPENFF → GROMACS")
    logging.info("Время: %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logging.info("Force field: %s", args.forcefield)
    logging.info("=" * 60)

    if not args.input.is_file():
        logging.error("Файл не найден: %s", args.input)
        sys.exit(1)

    if args.input.suffix.lower() != ".sdf":
        logging.error("Ожидается .sdf файл, получен: %s", args.input.suffix)
        sys.exit(1)

    molecules = load_sdf_molecules(args.input)
    if not molecules:
        logging.error("Не найдено молекул в файле %s", args.input)
        sys.exit(1)

    if args.limit:
        molecules = molecules[:args.limit]
        logging.info("Ограничение: обрабатываем %d молекул", args.limit)

    logging.info("Загружено %d молекул из %s", len(molecules), args.input.name)

    results: List[Tuple[str, str, bool]] = []
    name_counts: dict = {}

    for rdkit_mol, name in molecules:
        # Обработка дубликатов имён
        if name in name_counts:
            name_counts[name] += 1
            unique_name = name + "_" + str(name_counts[name])
        else:
            name_counts[name] = 1
            unique_name = name

        logging.info("Обработка молекулы: %s", unique_name)

        mol_dir = args.output / unique_name
        success, msg = parameterize_molecule(
            rdkit_mol=rdkit_mol,
            name=unique_name,
            output_dir=mol_dir,
            forcefield=args.forcefield,
            check_coords=not args.no_coord_check,
        )

        if success:
            logging.info("Молекула '%s' успешно параметризована", name)
        else:
            logging.error("Ошибка параметризации '%s': %s", name, msg)

        results.append((name, msg, success))

    # Итоговая статистика
    if results:
        success = sum(1 for _, _, ok in results if ok)
        total = len(results)

        logging.info("=" * 60)
        logging.info("ПАРАМЕТРИЗАЦИЯ ЗАВЕРШЕНА")
        logging.info("Всего молекул: %d | Успешно: %d | Ошибок: %d",
                     total, success, total - success)
        logging.info("Результаты сохранены в: %s", args.output)
        logging.info("=" * 60)

        if total - success > 0:
            logging.warning("Некоторые молекулы не удалось параметризовать."
                            " Подробности в логе.")
    else:
        logging.error("Нет результатов для обработки")
        sys.exit(1)


if __name__ == "__main__":
    main()
