"""Process ranked AlphaFold 3 structure predictions.

The input directory must be an AF3 prediction directory containing
``ranking_scores.csv`` and the ``seed-*_sample-*`` model directories. Processed
structures can be written as either PDB or CIF.

Examples:
    python process_AF3_preds.py --input-dir /path/to/af3/run \\
        --output-dir /path/to/results --output-format pdb
    python process_AF3_preds.py --input-dir /path/to/af3/run \\
        --output-format cif

If ``--output-dir`` is omitted, generated files are saved in ``--input-dir``.
"""
import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any

from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO


def save_structure(
    cif_file: Path, output_file: Path, output_format: str
) -> None:
    """Reads an AF3 CIF model and saves it in the requested format.

    Args:
        cif_file: Path to the source CIF structure produced by AF3.
        output_file: Path where the processed structure will be written.
        output_format: Structure format to write. Supported values are ``pdb``
            and ``cif``.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", str(cif_file))
    io = PDBIO() if output_format == "pdb" else MMCIFIO()
    io.set_structure(structure)
    io.save(str(output_file))


def create_parser() -> argparse.ArgumentParser:
    """Creates the command-line argument parser.

    Returns:
        A parser configured with input, output, and structure-format arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Process AF3 predictions, save ranked structures as PDB or CIF, "
            "and write model metadata."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        required=True,
        help="AF3 prediction directory containing ranking_scores.csv",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for generated structures and metadata (default: input directory)",
    )
    parser.add_argument(
        "-f",
        "--output-format",
        choices=("pdb", "cif"),
        default="pdb",
        help="Format for processed structure files",
    )
    return parser


def clean_path(path: Path) -> Path:
    """Expands the user-home component of a filesystem path.

    Args:
        path: Filesystem path that may begin with a tilde.

    Returns:
        The path with its user-home component expanded.
    """
    return Path(path).expanduser()


def resolve_paths(args: argparse.Namespace) -> tuple[Path, Path]:
    """Resolves input and output directories from command-line arguments.

    When no output directory is supplied, the input directory is also used as
    the output directory.

    Args:
        args: Parsed arguments containing ``input_dir`` and ``output_dir``.

    Returns:
        A pair containing the resolved input and output directory paths.
    """
    input_dir = clean_path(args.input_dir)
    output_dir = clean_path(args.output_dir) if args.output_dir else input_dir

    return input_dir, output_dir


def load_json(json_file: Path) -> Any:
    """Loads and decodes a JSON file.

    Args:
        json_file: Path to the JSON file to load.

    Returns:
        The Python object decoded from the JSON document.
    """
    with json_file.open("r", encoding="utf-8") as conf_file:
        return json.load(conf_file)


def process_predictions(
    input_dir: Path, output_dir: Path, output_format: str
) -> None:
    """Processes ranked AF3 models and writes each complete structure.

    The function reads model rankings and confidence values from the AF3 output
    directory. It writes ranked structure files and a ``processed_models.dat``
    metadata summary to the output directory.

    Args:
        input_dir: AF3 directory containing ``ranking_scores.csv`` and model
            subdirectories.
        output_dir: Directory where structures and metadata will be written.
        output_format: Structure format to write. Supported values are ``pdb``
            and ``cif``.
    """
    print("Processing AF3 models...")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Structure format: {output_format.upper()}")

    ranking_score_file = input_dir / "ranking_scores.csv"
    outputinfofile = output_dir / "processed_models.dat"
    output_dir.mkdir(parents=True, exist_ok=True)

    with ranking_score_file.open("r", encoding="utf-8") as infile, outputinfofile.open("w", encoding="utf-8") as fout:
        csv_reader = csv.reader(infile)
        next(csv_reader, None)  # skip the headers
        for i, row in enumerate(csv_reader):
            seed = row[0]
            sample = row[1]
            score = str(round(float(row[2]), 2))
            modeldir = f"seed-{seed}_sample-{sample}"

            summary_conf_file = input_dir / modeldir / "summary_confidences.json"
            conf_data = load_json(summary_conf_file)
            ptm = round(conf_data.get("ptm", 0), 3)
            iptm = round(conf_data.get("iptm", 0), 3)

            conf_file_path = input_dir / modeldir / "confidences.json"
            conf_data = load_json(conf_file_path)
            plddt_scores = conf_data.get("atom_plddts", [])
            avg_plddt = round(statistics.mean(plddt_scores), 2) if plddt_scores else 0

            fout.write(f"rank {i}: seed {seed}, sample {sample}, score {score}, ptm {ptm}, iptm {iptm}, avg_plddt {avg_plddt}\n")
            modelfilepath = input_dir / modeldir / "model.cif"
            outputfilepath = output_dir / f"rank_{i}.{output_format}"
            save_structure(modelfilepath, outputfilepath, output_format)


def main() -> None:
    """Parses command-line arguments and processes the AF3 predictions."""
    parser = create_parser()
    args = parser.parse_args()
    input_dir, output_dir = resolve_paths(args)
    if not input_dir.is_dir():
        parser.error(f"input directory does not exist: {input_dir}")
    if not (input_dir / "ranking_scores.csv").is_file():
        parser.error(
            f"input directory does not contain ranking_scores.csv: {input_dir}"
        )

    process_predictions(input_dir, output_dir, args.output_format)


if __name__ == "__main__":
    main()
