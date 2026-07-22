"""Tests FASTA header parsing and chain filtering for AF3 JSON creation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from create_input_jsonfile_AF3 import parse_fasta_record_ids
from create_input_jsonfile_AF3 import parse_rcsb_chain_ids
from create_input_jsonfile_AF3 import read_fasta


RCSB_FASTA = """>6BYY_1|Chains A, B, C, D|MEF2 CHIMERA|Homo sapiens (9606)
MGRKKIQITRIMDERNRQVTFTKRKFGLMKKAYELSVLCDCEIALIIFNSSNKL
>6BYY_2|Chains E, G[auth K]|DNA
AACTATTTATAAGA
>6BYY_3|Chains F, H[auth L]|DNA
TTCTTATAAATAGTT
"""


class FastaHeaderParsingTest(unittest.TestCase):
    """Tests RCSB and simplified FASTA header parsing."""

    def test_parses_rcsb_chain_list(self) -> None:
        """Parses comma-separated label chain IDs from an RCSB header."""
        description = "6BYY_1|Chains A, B, C, D|MEF2 CHIMERA"
        self.assertEqual(parse_rcsb_chain_ids(description), ["A", "B", "C", "D"])

    def test_removes_rcsb_author_chain_annotations(self) -> None:
        """Uses label chain IDs rather than bracketed author-chain IDs."""
        description = "6BYY_2|Chains E, G[auth K]|DNA"
        self.assertEqual(parse_rcsb_chain_ids(description), ["E", "G"])

    def test_preserves_multiletter_rcsb_chain_id(self) -> None:
        """Treats a multi-letter RCSB label as one chain ID."""
        description = "EXAMPLE_1|Chains AA|PROTEIN"
        self.assertEqual(parse_rcsb_chain_ids(description), "AA")

    def test_falls_back_to_simplified_header(self) -> None:
        """Retains support for the existing underscore-based header format."""
        self.assertEqual(
            parse_fasta_record_ids("6BZ1_AB", "6BZ1_AB"),
            ("6BZ1", ["A", "B"]),
        )


class ChainFilteringTest(unittest.TestCase):
    """Tests chain selection while reading RCSB FASTA records."""

    def setUp(self) -> None:
        """Creates a temporary RCSB-style FASTA file."""
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.fasta_path = Path(self.temporary_directory.name) / "6byy.fasta"
        self.fasta_path.write_text(RCSB_FASTA, encoding="utf-8")

    def tearDown(self) -> None:
        """Removes the temporary FASTA directory."""
        self.temporary_directory.cleanup()

    def test_reads_all_rcsb_chains(self) -> None:
        """Uses every label chain ID when no selection is supplied."""
        sequences = read_fasta(self.fasta_path)
        ids = [entry[next(iter(entry))]["id"] for entry in sequences]
        self.assertEqual(ids, [["A", "B", "C", "D"], ["E", "G"], ["F", "H"]])

    def test_selects_rcsb_chain_subset(self) -> None:
        """Filters copies within entities and omits unselected copies."""
        sequences = read_fasta(
            self.fasta_path,
            include_chains=["A", "B", "E", "F"],
        )
        ids = [entry[next(iter(entry))]["id"] for entry in sequences]
        self.assertEqual(ids, [["A", "B"], "E", "F"])

    def test_accepts_comma_separated_chain_selection(self) -> None:
        """Accepts comma-separated chain IDs in a selection item."""
        sequences = read_fasta(self.fasta_path, include_chains=["A,B,E,F"])
        ids = [entry[next(iter(entry))]["id"] for entry in sequences]
        self.assertEqual(ids, [["A", "B"], "E", "F"])

    def test_rejects_missing_requested_chain(self) -> None:
        """Reports requested chains that do not occur in the FASTA."""
        with self.assertRaisesRegex(ValueError, "Z"):
            read_fasta(self.fasta_path, include_chains=["A", "Z"])


if __name__ == "__main__":
    unittest.main()
