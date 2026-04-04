"""CLI-скрипт для массовой конвертации SMILES/SDF → PDBQT для AutoDock Vina."""
import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Tuple
from utils import setup_logging
from prepare_ligand import convert_smi_file, convert_sdf_file


def parse_args() -> argparse.Namespace:
    """Парсит аргументы командной строки."""
    parser = argparse.ArgumentParser(
        description="Конвертация лигандов в PDBQT для AutoDock Vina",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Путь к .smi/.sdf файлу или директории",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("output/ligands"),
        help="Директория для сохранения PDBQT-файлов",
    )
    parser.add_argument(
        "--flex-rings",
        "-f",
        action="store_true",
        help="Делать алифатические циклы гибкими",
    )
    parser.add_argument(
        "--min-ring-size",
        type=int,
        default=6,
        help="Минимальный размер цикла для гибкости (по умолчанию: 6)",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("logs/conversion_log.txt"),
        help="Путь к лог-файлу (по умолчанию: logs/conversion_log.txt)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Включить подробный вывод (DEBUG)",
    )

    return parser.parse_args()


def main() -> None:
    """Основная функция CLI."""
    args = parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    setup_logging(args.log_file, log_level)

    logging.info("=" * 60)
    logging.info("НАЧАЛО КОНВЕРТАЦИИ PDBQT")
    logging.info("Время: %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logging.info("=" * 60)

    args.output.mkdir(parents=True, exist_ok=True)

    all_results: List[Tuple[str, str, bool]] = []

    if args.input.is_file():
        suffix = args.input.suffix.lower()

        if suffix in (".smi", ".smiles"):
            logging.info("Обработка SMILES файла: %s", args.input.name)
            results = convert_smi_file(
                args.input, args.output, args.flex_rings, args.min_ring_size
            )
            all_results.extend(results)

        elif suffix == ".sdf":
            logging.info("Обработка SDF файла: %s", args.input.name)
            results = convert_sdf_file(
                args.input, args.output, args.flex_rings, args.min_ring_size
            )
            all_results.extend(results)

        else:
            logging.error("Неподдерживаемый формат файла: %s", suffix)
            sys.exit(1)

    elif args.input.is_dir():
        logging.info("Обработка директории: %s", args.input)
        for file_path in args.input.glob("*.smi"):
            results = convert_smi_file(
                file_path, args.output, args.flex_rings, args.min_ring_size
            )
            all_results.extend(results)
        for file_path in args.input.glob("*.sdf"):
            results = convert_sdf_file(
                file_path, args.output, args.flex_rings, args.min_ring_size
            )
            all_results.extend(results)

    else:
        logging.error("Указанный путь не существует: %s", args.input)
        sys.exit(1)

    # Итоговая статистика
    if all_results:
        success = sum(1 for _, _, ok in all_results if ok)
        total = len(all_results)

        logging.info("=" * 60)
        logging.info("КОНВЕРТАЦИЯ ЗАВЕРШЕНА")
        logging.info("Всего молекул: %d | Успешно: %d | Ошибок: %d",
                     total, success, total - success)
        logging.info("Результаты сохранены в: %s", args.output)
        logging.info("=" * 60)

        if total - success > 0:
            logging.warning("Некоторые молекулы не удалось конвертировать."
                            " Подробности в логе.")
    else:
        logging.error("Не найдено молекул для обработки")
        sys.exit(1)


if __name__ == "__main__":
    main()
