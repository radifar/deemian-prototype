from rdkit import Chem

from deemian.engine import model, processor


def test_load_molecule_processor(simple_yaml_parsed):
    specification = model.Workflow.parse_obj(simple_yaml_parsed)
    molecules = specification.load_molecule
    molecule_id = [key for key in molecules.keys()][0]

    loaded_molecules = processor.molecule_loader(molecules)
    loaded_molecule = loaded_molecules[molecule_id]
    get_atom_num = loaded_molecule.GetNumAtoms()

    assert isinstance(loaded_molecule, Chem.rdchem.Mol)
    assert get_atom_num == 220


def test_selection_processor():
    # TODO: Create fixture for RDKit.mol of PDB example
    # TODO: Selection Model --> fill with selection_string, preposition, identifier
    # TODO: Add Selection method on Model that generate the corresponding selection (RDKit.mol)
    pass


def test_measurement_processor():
    pass
