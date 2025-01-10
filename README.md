# Corse Grain Converter for PDB Files

The coarse grain format simplifies the representation of protein structures by reducing the number of atoms to a smaller set of pseudo-atoms. This allows for faster computations and easier analysis of large biomolecular systems. In this format, groups of atoms are represented by a single pseudo-atom, preserving the overall structure and essential features of the protein.

## Requirements

- biopython==1.84
- numpy==2.1.3

## Installation

```
pip install -r requirements.txt
```

## Usage

**Reconstruct**
```
py coarse_grain_format.py --coarse_grain_file <COARSE_GRAIN_PDB_FILE_PATH>
```

**Generate Coarse Grain and Reconstruct**
```
py coarse_grain_format.py --pdb_file_input <PDB_FILE_PATH>
```
