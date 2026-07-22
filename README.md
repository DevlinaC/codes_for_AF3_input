# AlphaFold 3 Input and Output Utilities

Small command-line utilities for preparing AlphaFold 3 inputs and processing
prediction outputs. Run any script with `--help` for its complete options; the
source files also contain detailed docstrings.

## Installation

Python 3.10 or newer is required. Install the Python dependency with:

```bash
python -m pip install -r requirements.txt
```

## Scripts

### `create_input_jsonfile_AF3.py`

Creates an AlphaFold 3 input JSON file from a FASTA file. It supports both RCSB
FASTA headers and simplified headers such as `>6BZ1_AB`, can select particular
chains, and can configure model seeds, MSA-free input, and template-free input.

```bash
python create_input_jsonfile_AF3.py input.fasta --output input.json
```

### `modify_json.py`

Modifies MSA and template settings in an existing AlphaFold 3 input JSON file.
It supports MSA-free input, custom A3M files, template-free input, and custom
template JSON.

```bash
python modify_json.py input.json --single-seq --no-templates
```

### `process_AF3_preds.py`

Processes ranked AlphaFold 3 predictions, saves each complete structure as PDB
or CIF, and writes confidence values to `processed_models.dat`.

```bash
python process_AF3_preds.py --input-dir /path/to/predictions \
    --output-dir /path/to/results --output-format cif
```

## Tests

Run the FASTA parsing and chain-selection tests with:

```bash
python -m unittest -v test_create_input_jsonfile_AF3.py
```

## Development disclosure

The initial scripts in this repository were written by the project maintainer.
OpenAI Codex was subsequently used as a development assistant for code review,
refactoring, command-line interface design, documentation, and test generation.
All AI-assisted changes were reviewed, tested, and accepted by the maintainer,
who remains responsible for the published code.

## License

This project is available under the [MIT License](LICENSE).

AlphaFold 3 input format reference:
[google-deepmind/alphafold3 input documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)
