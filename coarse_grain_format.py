import pathlib
import argparse

from Bio.PDB import (
    PDBParser,
    Structure,
    Select,
    Atom,
    PDBIO,
    Superimposer,
    Chain,
    Model,
    Residue,
)

TEMPLATES_PATH = pathlib.Path(__file__).parent / "templates"

PURINE_ATOMS = {"N9", "C2", "C6"}
PYRIMIDINE_ATOMS = {"N1", "C2", "C4"}
BACKBONE_ATOMS = {"P", "C4'"}

ATOMS_SUBSET = set().union(PURINE_ATOMS, PYRIMIDINE_ATOMS, BACKBONE_ATOMS)

RESIDUES = {"A", "G", "C", "U"}


class CoarseGrainStructure(Select):
    def accept_atom(self, atom: Atom.Atom) -> int:
        if (
            atom.get_name() in ATOMS_SUBSET
            and atom.get_parent().get_resname() in RESIDUES
        ):
            return 1
        return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse_grain_input", type=pathlib.Path, default=None)
    parser.add_argument("--pdb_file_input", type=pathlib.Path, default=None)
    return parser.parse_args()


def read_pdb(file_path: pathlib.Path, file_name: str) -> Structure.Structure:
    parser = PDBParser()
    structure = parser.get_structure(file_name, file_path)
    return structure


def save_structure(
    structure: Structure.Structure,
    file_path: pathlib.Path,
    as_coarse_grain: bool = False,
) -> None:
    io = PDBIO()
    io.set_structure(structure)
    if as_coarse_grain:
        io.save(file=str(file_path), select=CoarseGrainStructure())
    else:
        io.save(file=str(file_path))


def reconstruct_pdb(
    structure: Structure.Structure, templates: dict[str, Residue.Residue]
) -> Structure.Structure:
    superimposer = Superimposer()
    reconstructed_structure = Structure.Structure("reconstructed")

    i = 1
    for model in structure:
        if not any(m.id == model.id for m in reconstructed_structure):
            reconstructed_structure.add(Model.Model(model.id))
        for chain in model:
            if not any(c.id == chain.id for c in reconstructed_structure[model.id]):
                reconstructed_structure[model.id].add(Chain.Chain(chain.id))
            for residue in chain:
                template_residue = templates[residue.get_resname()].copy()

                template_residue.id = (" ", i, " ")
                i += 1

                target_atoms = list(residue.get_atoms())
                template_atoms = [template_residue[a.get_name()] for a in target_atoms]

                superimposer.set_atoms(target_atoms, template_atoms)
                superimposer.apply(template_residue)

                reconstructed_structure[model.id][chain.id].add(template_residue)

    return reconstructed_structure


def parse_residue_templates(templates_path: pathlib.Path) -> dict[str, Residue.Residue]:
    templates = {}
    for template in templates_path.iterdir():
        template_name = template.stem
        templates[template.stem] = next(
            PDBParser().get_structure(template_name, template).get_residues()
        )
    return templates


def main():
    args = parse_args()
    if not args.coarse_grain_input and not args.pdb_file_input:
        raise ValueError("Please provide pdb file or coarse grain input file")

    if args.pdb_file_input and args.coarse_grain_input:
        raise ValueError("Please provide only one input file")

    if not args.pdb_file_input:
        file_name = args.pdb_file_input.name.split(".")[0].upper()
        structure = read_pdb(args.coarse_grain_input, file_name)

        coarse_grain_file_name = f"{file_name}_coarse_grain.pdb"
        coarse_grain_file_path = args.coarse_grain_input.parent / coarse_grain_file_name
        save_structure(structure, coarse_grain_file_path, as_coarse_grain=True)

    templates = parse_residue_templates(TEMPLATES_PATH)
    coarse_grain_structure = PDBParser().get_structure(
        coarse_grain_file_name, coarse_grain_file_path
    )
    reconstructed_structure = reconstruct_pdb(coarse_grain_structure, templates)
    reco_file_name = f"{file_name}_reconstructed.pdb"
    save_structure(
        reconstructed_structure, args.coarse_grain_input.parent / reco_file_name
    )


if __name__ == "__main__":
    main()
