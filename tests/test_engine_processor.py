from rdkit import Chem

from deemian.engine import model


def test_molecule_loader(simple_yaml_parsed):
    workflow = model.Workflow.parse_obj(simple_yaml_parsed)
    molecule_id = list(workflow.load_molecule.keys())[0]

    workflow.molecule_loader()
    loaded_molecule = workflow.load_molecule[molecule_id]
    atom_num = loaded_molecule.GetNumAtoms()

    assert isinstance(loaded_molecule, Chem.rdchem.Mol)
    assert atom_num == 220


def test_selection_processor():
    # TODO: Create fixture for RDKit.mol of PDB example
    # TODO: Selection Model --> fill with selection_string, preposition, identifier
    # TODO: Add Selection method on Model that generate the corresponding selection (RDKit.mol)
    pass


def test_measurement_processor():
    pass
