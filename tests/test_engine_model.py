from deemian.engine import model


def test_interaction_model(simple_yaml_parsed):
    interaction = simple_yaml_parsed["interaction"]

    result = model.Interaction.parse_obj(interaction)

    assert result.measure == "nonpolar on protein"


def test_presentation_model(simple_yaml_parsed):
    presentation = simple_yaml_parsed["presentation"]

    result = model.Presentation.parse_obj(presentation)

    assert result.granularity == "atomic"
    assert result.field == ["atom_name", "residue_name", "chain", "distance"]
    assert result.output_file == "nonpolar_zinc_finger.txt"


def test_presentation_model_field_list(simple_yaml_parsed, monkeypatch):
    presentation = simple_yaml_parsed["presentation"]
    monkeypatch.setitem(presentation, "field", ["atom_name", "residue_name", "chain", "distance"])

    result = model.Presentation.parse_obj(presentation)

    assert result.field == ["atom_name", "residue_name", "chain", "distance"]


def test_specification_model(simple_yaml_parsed):
    interaction = simple_yaml_parsed["interaction"]
    interaction_model = model.Interaction.parse_obj(interaction)

    presentation = simple_yaml_parsed["presentation"]
    presentation_model = model.Presentation.parse_obj(presentation)

    expected_molecule = {"zinc_finger": "tests/data/1znm.pdb"}
    expected_selection = {"protein": "protein in zinc_finger"}

    result = model.Specification.parse_obj(simple_yaml_parsed)

    assert result.load_molecule == expected_molecule
    assert result.selection == expected_selection
    assert result.interaction == interaction_model
    assert result.presentation == presentation_model
