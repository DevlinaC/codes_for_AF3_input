"""
Create AlphaFold 3 input JSON from a FASTA file.

The script accepts simplified headers with chain IDs after the final underscore:

    >6BZ1_C
    TTCTTATAAATAGTT
    >6BZ1_D
    AACTATTTATAAGA
    >6BZ1_AB
    MGRKKIQITRIMDERNRQVTFTKRKFGLMKKAYELSVLCDCEIALIIFNSSNKLFQYASTDMDKVLLKYTEYSEPHESRTNTVILETLKRREHRG

The chain label is used as the AlphaFold 3 entity ID. Multi-letter labels such
as "AB" are interpreted as homomeric copies and become ["A", "B"].

RCSB FASTA headers are also detected automatically:

    >6BYY_1|Chains A, B, C, D|MEF2 CHIMERA|Homo sapiens (9606)

Use ``--include-chains`` to select a subset of chains from either header style.

Usage:
    python create_input_jsonfile_AF3.py /path/to/sequences.fasta
    python create_input_jsonfile_AF3.py /path/to/sequences.fasta my_job \
        --output-dir /path/to/output
    python create_input_jsonfile_AF3.py /path/to/rcsb.fasta \
        --include-chains A B E F

Reference:
    AlphaFold 3 input documentation:
    https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from Bio import SeqIO


DNA_ALPHABET = set("ACGT")
RNA_ALPHABET = set("ACGU")
NUCLEIC_ACID_ALPHABET = set("ACGTU")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")
RCSB_CHAINS_PATTERN = re.compile(r"(?:^|\|)Chains?\s+([^|]+)", re.IGNORECASE)
RCSB_AUTH_CHAIN_PATTERN = re.compile(r"\s*\[auth\s+[^\]]+\]", re.IGNORECASE)


def annotate_seq(seq_string: str) -> str:
    """Infers whether a sequence is protein, DNA, or RNA.

    Args:
        seq_string: Sequence containing standard one-letter residue codes.

    Returns:
        One of ``protein``, ``dna``, or ``rna``.

    Raises:
        ValueError: If the sequence is empty, mixes thymine and uracil, or
            contains unsupported residue codes.
    """
    sequence = seq_string.upper()
    letters = set(sequence)

    if not sequence:
        raise ValueError("Found an empty sequence.")

    if letters <= NUCLEIC_ACID_ALPHABET:
        if "T" in letters and "U" in letters:
            raise ValueError(
                f"Cannot infer DNA/RNA type for sequence with both T and U: {sequence}"
            )
        if letters <= RNA_ALPHABET and "U" in letters:
            return "rna"
        if letters <= DNA_ALPHABET:
            return "dna"

    invalid = sorted(letters - PROTEIN_ALPHABET)
    if invalid:
        invalid_text = ", ".join(invalid)
        raise ValueError(
            "Unsupported sequence character(s) for AF3 input: "
            f"{invalid_text}. Use standard protein, DNA, or RNA one-letter codes."
        )

    return "protein"


def parse_chain_ids(chain_label: str) -> str | list[str]:
    """Converts a FASTA chain label into an AF3 entity ID field.

    Multi-character labels are treated as multiple homomeric chain IDs. Comma,
    slash, plus, and colon can also be used as explicit ID separators.

    Args:
        chain_label: Chain label extracted from a FASTA record identifier.

    Returns:
        A single chain ID or a list of chain IDs for homomeric copies.

    Raises:
        ValueError: If the label is empty or contains invalid chain IDs.
    """
    clean_label = chain_label.strip().upper()
    if not clean_label:
        raise ValueError("FASTA record has an empty chain label.")

    if re.search(r"[,/+:]", clean_label):
        chain_ids = [
            part.strip()
            for part in re.split(r"[,/+:]+", clean_label)
            if part.strip()
        ]
    elif len(clean_label) > 1:
        chain_ids = list(clean_label)
    else:
        chain_ids = [clean_label]

    invalid = [
        chain_id
        for chain_id in chain_ids
        if not re.fullmatch(r"[A-Z0-9]+", chain_id)
    ]
    if invalid:
        raise ValueError(f"Invalid AF3 chain id(s): {', '.join(invalid)}")

    return chain_ids[0] if len(chain_ids) == 1 else chain_ids


def parse_record_id(record_id: str) -> tuple[str | None, str | list[str]]:
    """Extracts an optional prefix and chain IDs from a FASTA record ID.

    The text after the final underscore is interpreted as the chain label. If
    no underscore exists, the complete record ID is treated as the label.

    Args:
        record_id: Identifier from a FASTA record.

    Returns:
        A pair containing the optional identifier prefix and parsed chain IDs.
    """
    if "_" in record_id:
        prefix, chain_label = record_id.rsplit("_", 1)
    else:
        prefix, chain_label = None, record_id

    return prefix, parse_chain_ids(chain_label)


def parse_rcsb_chain_ids(description: str) -> str | list[str] | None:
    """Extracts chain IDs from an RCSB FASTA description.

    Author-chain annotations such as ``[auth K]`` are removed in favor of the
    RCSB label chain IDs used before those annotations.

    Args:
        description: Complete FASTA header text without the leading ``>``.

    Returns:
        A single chain ID, a list of chain IDs, or ``None`` when the description
        does not contain an RCSB ``Chains`` field.
    """
    match = RCSB_CHAINS_PATTERN.search(description)
    if match is None:
        return None

    chain_field = RCSB_AUTH_CHAIN_PATTERN.sub("", match.group(1))
    chain_ids = [
        part.strip().upper()
        for part in chain_field.split(",")
        if part.strip()
    ]
    invalid = [
        chain_id
        for chain_id in chain_ids
        if not re.fullmatch(r"[A-Z0-9]+", chain_id)
    ]
    if invalid:
        raise ValueError(f"Invalid RCSB chain id(s): {', '.join(invalid)}")
    if not chain_ids:
        raise ValueError("RCSB FASTA header has an empty Chains field.")
    return chain_ids[0] if len(chain_ids) == 1 else chain_ids


def parse_fasta_record_ids(
    record_id: str, description: str
) -> tuple[str | None, str | list[str]]:
    """Extracts entity and chain IDs from an RCSB or simplified FASTA header.

    Args:
        record_id: Identifier parsed from the first token of the FASTA header.
        description: Complete FASTA header text without the leading ``>``.

    Returns:
        A pair containing the optional entity prefix and parsed chain IDs.
    """
    rcsb_chain_ids = parse_rcsb_chain_ids(description)
    if rcsb_chain_ids is not None:
        entity_id = record_id.split("|", 1)[0]
        return entity_id, rcsb_chain_ids
    return parse_record_id(record_id)


def as_id_list(chain_id: str | list[str]) -> list[str]:
    """Normalizes an AF3 entity ID field to a list.

    Args:
        chain_id: A single chain ID or list of chain IDs.

    Returns:
        A list containing all chain IDs.
    """
    if isinstance(chain_id, list):
        return chain_id
    return [chain_id]


def normalize_chain_filter(include_chains: list[str] | None) -> set[str] | None:
    """Normalizes optional chain-selection values from the API or CLI.

    Comma-separated values are accepted in addition to separate list items.

    Args:
        include_chains: Chain IDs to retain, or ``None`` to retain all chains.

    Returns:
        An uppercase set of selected chain IDs, or ``None`` when no filter was
        supplied.

    Raises:
        ValueError: If the supplied filter contains no valid chain IDs.
    """
    if include_chains is None:
        return None

    normalized = {
        chain_id.strip().upper()
        for value in include_chains
        for chain_id in value.split(",")
        if chain_id.strip()
    }
    if not normalized:
        raise ValueError("At least one chain ID must be supplied.")

    invalid = sorted(
        chain_id
        for chain_id in normalized
        if not re.fullmatch(r"[A-Z0-9]+", chain_id)
    )
    if invalid:
        raise ValueError(f"Invalid chain ID(s): {', '.join(invalid)}")
    return normalized


def validate_unique_ids(sequences: list[dict[str, dict[str, Any]]]) -> None:
    """Ensures every AF3 entity ID is unique across all sequences.

    Args:
        sequences: AF3 sequence entries to validate.

    Raises:
        ValueError: If a chain ID occurs in more than one sequence entry.
    """
    seen: dict[str, str] = {}
    for sequence_entry in sequences:
        molecule_type, payload = next(iter(sequence_entry.items()))
        for chain_id in as_id_list(payload["id"]):
            if chain_id in seen:
                raise ValueError(
                    f"Duplicate AF3 chain id '{chain_id}' found in "
                    f"{seen[chain_id]} and {molecule_type} entries."
                )
            seen[chain_id] = molecule_type


def sequence_payload(
    molecule_type: str,
    chain_id: str | list[str],
    sequence: str,
    single_seq: bool = False,
    no_templates: bool = False,
) -> dict[str, Any]:
    """Builds an AF3 payload for one molecular entity.

    Args:
        molecule_type: AF3 entity type: ``protein``, ``dna``, or ``rna``.
        chain_id: A single chain ID or IDs for identical homomeric copies.
        sequence: Polymer sequence in one-letter notation.
        single_seq: Whether to disable MSA input for proteins and RNA.
        no_templates: Whether to disable template input for proteins.

    Returns:
        A JSON-compatible AF3 molecular entity payload.
    """
    payload: dict[str, Any] = {
        "id": chain_id,
        "sequence": sequence,
        "modifications": [],
    }

    if molecule_type == "protein":
        if single_seq:
            payload["unpairedMsa"] = ""
            payload["pairedMsa"] = ""
        if no_templates:
            payload["templates"] = []
    elif molecule_type == "rna" and single_seq:
        payload["unpairedMsa"] = ""

    return payload


def read_fasta(
    inputfile: str | Path,
    single_seq: bool = False,
    no_templates: bool = False,
    include_chains: list[str] | None = None,
) -> list[dict[str, dict[str, Any]]]:
    """Reads FASTA records and creates AF3 sequence entries.

    Args:
        inputfile: Path to the input FASTA file.
        single_seq: Whether to disable MSA input for proteins and RNA.
        no_templates: Whether to disable template input for proteins.
        include_chains: Optional chain IDs to retain from the FASTA records.

    Returns:
        JSON-compatible AF3 sequence entries in FASTA record order.

    Raises:
        ValueError: If the FASTA is empty, contains unsupported sequence
            characters, or produces duplicate chain IDs.
    """
    sequences: list[dict[str, dict[str, Any]]] = []
    input_path = Path(inputfile)
    selected_chains = normalize_chain_filter(include_chains)
    matched_chains: set[str] = set()

    with input_path.open(encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            _, chain_id = parse_fasta_record_ids(record.id, record.description)
            record_chain_ids = as_id_list(chain_id)
            if selected_chains is not None:
                retained_ids = [
                    item for item in record_chain_ids if item in selected_chains
                ]
                matched_chains.update(retained_ids)
                if not retained_ids:
                    continue
                chain_id = retained_ids[0] if len(retained_ids) == 1 else retained_ids

            sequence = str(record.seq).upper()
            molecule_type = annotate_seq(sequence)
            payload = sequence_payload(
                molecule_type=molecule_type,
                chain_id=chain_id,
                sequence=sequence,
                single_seq=single_seq,
                no_templates=no_templates,
            )
            sequences.append({molecule_type: payload})

    if selected_chains is not None:
        missing_chains = sorted(selected_chains - matched_chains)
        if missing_chains:
            raise ValueError(
                "Requested chain ID(s) not found in FASTA: "
                f"{', '.join(missing_chains)}"
            )

    if not sequences:
        raise ValueError(f"No FASTA records found in {input_path}")

    validate_unique_ids(sequences)
    return sequences


def build_af3_json(
    job_name: str,
    fastafile: str | Path,
    model_seeds: list[int],
    version: int = 1,
    single_seq: bool = False,
    no_templates: bool = False,
    include_chains: list[str] | None = None,
) -> dict[str, Any]:
    """Creates a complete AF3 JSON-compatible dictionary.

    Args:
        job_name: Name used to identify the AF3 prediction job.
        fastafile: Path to the source FASTA file.
        model_seeds: Random seeds AF3 should use for model inference.
        version: AF3 input JSON schema version.
        single_seq: Whether to disable MSA input for proteins and RNA.
        no_templates: Whether to disable template input for proteins.
        include_chains: Optional chain IDs to retain from the FASTA records.

    Returns:
        A JSON-compatible mapping containing the AF3 job configuration.
    """
    return {
        "name": job_name,
        "modelSeeds": model_seeds,
        "sequences": read_fasta(
            fastafile,
            single_seq=single_seq,
            no_templates=no_templates,
            include_chains=include_chains,
        ),
        "bondedAtomPairs": None,
        "userCCD": None,
        "dialect": "alphafold3",
        "version": version,
    }


def write_json(
    job_name: str,
    fastafile: str | Path,
    output_path: str | Path | None = None,
    output_dir: str | Path | None = None,
    model_seeds: list[int] | None = None,
    version: int = 1,
    single_seq: bool = False,
    no_templates: bool = False,
    include_chains: list[str] | None = None,
) -> Path:
    """Writes an AF3 input JSON file.

    Args:
        job_name: Name used to identify the AF3 prediction job.
        fastafile: Path to the source FASTA file.
        output_path: Exact destination path. Overrides ``output_dir``.
        output_dir: Destination directory used when ``output_path`` is absent.
        model_seeds: Random seeds AF3 should use for model inference. Defaults
            to seeds 1 through 5 when omitted.
        version: AF3 input JSON schema version.
        single_seq: Whether to disable MSA input for proteins and RNA.
        no_templates: Whether to disable template input for proteins.
        include_chains: Optional chain IDs to retain from the FASTA records.

    Returns:
        Path to the JSON file that was written.
    """
    seeds = model_seeds if model_seeds is not None else [1, 2, 3, 4, 5]

    if output_path is None:
        output_directory = Path(output_dir) if output_dir else Path(".")
        output_path = output_directory / f"{job_name}.json"
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    af3_json = build_af3_json(
        job_name=job_name,
        fastafile=fastafile,
        model_seeds=seeds,
        version=version,
        single_seq=single_seq,
        no_templates=no_templates,
        include_chains=include_chains,
    )

    with output_file.open("w", encoding="utf-8") as handle:
        json.dump(af3_json, handle, indent=1)
        handle.write("\n")

    return output_file


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments.

    Returns:
        The parsed FASTA, output, seed, and AF3 configuration arguments.
    """
    parser = argparse.ArgumentParser(
        description="Create an AlphaFold 3 JSON input file from FASTA."
    )
    parser.add_argument("fastafile", help="Input FASTA file.")
    parser.add_argument(
        "job_name",
        nargs="?",
        help="AF3 job name and output basename (default: FASTA filename stem).",
    )
    parser.add_argument(
        "--output",
        help="Exact output JSON path. Overrides --output_dir.",
    )
    parser.add_argument(
        "--output-dir",
        "--output_dir",
        help="Directory for the output JSON when --output is not set.",
    )
    parser.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=[1, 2, 3, 4, 5],
        help="One or more AF3 model seeds. Default: 1 2 3 4 5.",
    )
    parser.add_argument(
        "--version",
        type=int,
        default=1,
        help="AF3 JSON version. Default: 1.",
    )
    parser.add_argument(
        "--single-seq",
        "--single_seq",
        action="store_true",
        help="Set protein/RNA MSA fields to empty strings for MSA-free input.",
    )
    parser.add_argument(
        "--no-templates",
        "--no_templates",
        action="store_true",
        help="Set protein templates to an empty list.",
    )
    parser.add_argument(
        "--include-chains",
        "--include_chains",
        nargs="+",
        metavar="CHAIN",
        help="Only include these chain IDs; accepts spaces or commas.",
    )
    return parser.parse_args()


def main() -> None:
    """Creates an AF3 input JSON file from command-line arguments."""
    args = parse_args()
    fasta_path = Path(args.fastafile).expanduser()
    if not fasta_path.is_file():
        raise SystemExit(f"Input FASTA file does not exist: {fasta_path}")
    job_name = args.job_name or fasta_path.stem
    output_file = write_json(
        job_name=job_name,
        fastafile=fasta_path,
        output_path=args.output,
        output_dir=args.output_dir,
        model_seeds=args.seeds,
        version=args.version,
        single_seq=args.single_seq,
        no_templates=args.no_templates,
        include_chains=args.include_chains,
    )
    print(f"Wrote {output_file}")


if __name__ == "__main__":
    main()
