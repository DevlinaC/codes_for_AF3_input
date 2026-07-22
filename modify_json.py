"""Modify MSA and template settings in an AlphaFold 3 input JSON file.

The script accepts an exact input JSON path instead of assuming a PDB-based
filename or a particular directory layout. By default, the modified file is
written beside the input file with ``_modified`` appended to its stem.

Typical usage examples:

    python modify_json.py input.json --single-seq --no-templates
    python modify_json.py input.json --output output.json \
        --custom-msa-file alignment.a3m
    python modify_json.py input.json \
        --custom-templates templates.json

Custom templates can be supplied as a JSON list or as a path to a JSON file
containing a list.

Reference:
    AlphaFold 3 input documentation:
    https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def remove_invalid_msa_keys(chain: dict[str, Any]) -> None:
    """Removes incorrectly capitalized MSA keys from a chain payload.

    Args:
        chain: AF3 protein or RNA payload to modify in place.
    """
    chain.pop("unpairedMSA", None)
    chain.pop("pairedMSA", None)


def set_protein_msa(
    protein: dict[str, Any], unpaired_msa: str, paired_msa: str
) -> None:
    """Sets inline unpaired and paired MSA values for a protein chain.

    Existing MSA path fields and incorrectly capitalized MSA keys are removed
    so the resulting payload contains only the supported inline MSA fields.

    Args:
        protein: AF3 protein payload to modify in place.
        unpaired_msa: A3M-formatted unpaired MSA, or an empty string.
        paired_msa: A3M-formatted paired MSA, or an empty string.
    """
    remove_invalid_msa_keys(protein)
    protein.pop("unpairedMsaPath", None)
    protein.pop("pairedMsaPath", None)
    protein["unpairedMsa"] = unpaired_msa
    protein["pairedMsa"] = paired_msa


def set_rna_msa(rna: dict[str, Any], unpaired_msa: str) -> None:
    """Sets an inline unpaired MSA value for an RNA chain.

    Args:
        rna: AF3 RNA payload to modify in place.
        unpaired_msa: A3M-formatted unpaired MSA, or an empty string.
    """
    remove_invalid_msa_keys(rna)
    rna.pop("unpairedMsaPath", None)
    rna["unpairedMsa"] = unpaired_msa


def modify_json(
    json_data: dict[str, Any],
    single_seq: bool = False,
    no_templates: bool = False,
    custom_msa: str | None = None,
    custom_templates: list[Any] | None = None,
) -> dict[str, Any]:
    """Modifies MSA and template settings in AF3 JSON data.

    Protein and RNA payloads receive MSA changes. Template changes apply only
    to protein payloads, as required by the AF3 input schema. The input mapping
    is modified in place and returned for convenience.

    Args:
        json_data: Decoded AF3 input JSON to modify.
        single_seq: Whether to disable MSA input by setting MSA fields to empty
            strings.
        no_templates: Whether to disable protein templates with an empty list.
        custom_msa: Optional inline A3M content to assign to supported MSA
            fields. Ignored when ``single_seq`` is true.
        custom_templates: Optional template objects to assign to each protein.
            Ignored when ``no_templates`` is true.

    Returns:
        The same decoded JSON mapping after modification.
    """
    for sequence_entry in json_data.get("sequences", []):
        if "protein" in sequence_entry:
            protein = sequence_entry["protein"]
            remove_invalid_msa_keys(protein)
            if single_seq:
                set_protein_msa(protein, unpaired_msa="", paired_msa="")
            elif custom_msa is not None:
                set_protein_msa(
                    protein,
                    unpaired_msa=custom_msa,
                    paired_msa=custom_msa,
                )

            if no_templates:
                protein["templates"] = []
            elif custom_templates is not None:
                protein["templates"] = custom_templates

        elif "rna" in sequence_entry:
            rna = sequence_entry["rna"]
            remove_invalid_msa_keys(rna)
            if single_seq:
                set_rna_msa(rna, unpaired_msa="")
            elif custom_msa is not None:
                set_rna_msa(rna, unpaired_msa=custom_msa)

    return json_data


def load_json(json_file: Path) -> dict[str, Any]:
    """Loads an AF3 JSON document from disk.

    Args:
        json_file: Path to the input JSON document.

    Returns:
        The decoded top-level JSON object.

    Raises:
        ValueError: If the JSON document does not contain a top-level object.
    """
    with json_file.open("r", encoding="utf-8") as input_handle:
        data = json.load(input_handle)
    if not isinstance(data, dict):
        raise ValueError("AF3 input JSON must contain a top-level object.")
    return data


def load_text(text_file: Path) -> str:
    """Loads text content from a file.

    Args:
        text_file: Path to the text file to read.

    Returns:
        The complete text content of the file.
    """
    return text_file.read_text(encoding="utf-8")


def parse_template_source(value: str) -> list[Any]:
    """Parses templates from a JSON string or JSON file path.

    Args:
        value: JSON list text or a path to a JSON file containing a list.

    Returns:
        The decoded list of AF3 template objects.

    Raises:
        argparse.ArgumentTypeError: If the value is invalid JSON or does not
            decode to a list.
    """
    candidate_path = Path(value).expanduser()
    try:
        if candidate_path.is_file():
            parsed = json.loads(candidate_path.read_text(encoding="utf-8"))
        else:
            parsed = json.loads(value)
    except (OSError, json.JSONDecodeError) as error:
        raise argparse.ArgumentTypeError(
            "custom templates must be a JSON list or a path to one"
        ) from error

    if not isinstance(parsed, list):
        raise argparse.ArgumentTypeError("custom templates must decode to a list")
    return parsed


def default_output_path(input_file: Path) -> Path:
    """Builds the default output path for an input JSON file.

    Args:
        input_file: Path to the source JSON document.

    Returns:
        A sibling path whose stem ends in ``_modified``.
    """
    return input_file.with_name(f"{input_file.stem}_modified{input_file.suffix}")


def write_json(json_data: dict[str, Any], output_file: Path) -> None:
    """Writes AF3 JSON data to disk.

    Parent directories are created when needed.

    Args:
        json_data: Decoded JSON object to serialize.
        output_file: Destination path for the modified JSON document.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8") as output_handle:
        json.dump(json_data, output_handle, indent=1)
        output_handle.write("\n")


def create_parser() -> argparse.ArgumentParser:
    """Creates the command-line argument parser.

    Returns:
        A parser configured for AF3 JSON paths and modification settings.
    """
    parser = argparse.ArgumentParser(
        description="Modify MSA and template settings in an AF3 input JSON file."
    )
    parser.add_argument(
        "input_json",
        type=Path,
        help="Exact path to the AF3 input JSON file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output JSON path (default: INPUT_JSON with _modified appended).",
    )

    msa_group = parser.add_mutually_exclusive_group()
    msa_group.add_argument(
        "--single-seq",
        "--single_seq",
        action="store_true",
        help="Disable MSA input for protein and RNA chains.",
    )
    msa_group.add_argument(
        "--custom-msa-file",
        "--custom_msa_file",
        type=Path,
        help="A3M file whose content will be embedded in protein/RNA MSA fields.",
    )

    template_group = parser.add_mutually_exclusive_group()
    template_group.add_argument(
        "--no-templates",
        "--no_templates",
        action="store_true",
        help="Disable templates for protein chains.",
    )
    template_group.add_argument(
        "--custom-templates",
        "--custom_templates",
        type=parse_template_source,
        help="Template JSON list or path to a JSON file containing that list.",
    )
    return parser


def main() -> None:
    """Parses arguments, modifies an AF3 JSON file, and writes the result."""
    parser = create_parser()
    args = parser.parse_args()

    input_file = args.input_json.expanduser()
    if not input_file.is_file():
        parser.error(f"input JSON file does not exist: {input_file}")

    output_file = (
        args.output.expanduser()
        if args.output is not None
        else default_output_path(input_file)
    )
    custom_msa = None
    if args.custom_msa_file is not None:
        msa_file = args.custom_msa_file.expanduser()
        if not msa_file.is_file():
            parser.error(f"custom MSA file does not exist: {msa_file}")
        custom_msa = load_text(msa_file)

    json_data = load_json(input_file)
    modified_data = modify_json(
        json_data,
        single_seq=args.single_seq,
        no_templates=args.no_templates,
        custom_msa=custom_msa,
        custom_templates=args.custom_templates,
    )
    write_json(modified_data, output_file)
    print(f"Wrote {output_file}")


if __name__ == "__main__":
    main()
