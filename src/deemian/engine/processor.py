from rdkit import Chem


def molecule_loader(molecules: dict):
    molecule_loader = dict()

    for molecule in molecules.keys():
        file_name = molecules[molecule]
        rdkit_mol = Chem.MolFromPDBFile(file_name)

        molecule_loader.update({molecule: rdkit_mol})

    return molecule_loader
