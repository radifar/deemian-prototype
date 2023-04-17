from deemian.engine import cfg


def test_parser_protein():
    """
    >>> import rich
    >>> from deemian.engine.cfg import parser
    >>> selection = parser("protein in zinc_finger")
    >>> rich.print(selection)
    start
    └── selection
        ├── selection_string
        │   └── protein
        ├── in
        └── zinc_finger
    >>> selection
    Tree(Token('RULE', 'start'),
      [Tree(Token('RULE', 'selection'),
        [Tree(Token('RULE', 'selection_string'),
          [Token('MACRO', 'protein')]),
        Token('PREPOSITION', 'in'),
        Token('IDENTIFIER', 'zinc_finger')])])
    """
    selection = "protein in zinc_finger"

    result = cfg.parser(selection)
    selection = result.children[0]
    selection_string, preposition, identifier = selection.children
    macro = selection_string.children[0]

    assert selection.data == "selection"
    assert selection_string.data == "selection_string"
    assert macro == "protein"
    assert preposition == "in"
    assert identifier == "zinc_finger"


def test_parser_selection_protein_or_metal():
    """
    >>> import rich
    >>> from deemian.engine.cfg import parser
    >>> selection = parser("protein OR metal in zinc_finger")
    >>> rich.print(selection)
    start
    └── selection
        ├── selection_string
        │   └── protein
        │   └── OR
        │   └── metal
        ├── in
        └── zinc_finger
    >>> selection
    Tree(Token('RULE', 'start'),
      [Tree(Token('RULE', 'selection'),
        [Tree(Token('RULE', 'selection_string'),
          [Token('MACRO', 'protein'),
          Token('LOGIC', 'OR'),
          Token('MACRO', 'metal')]),
        Token('PREPOSITION', 'in'),
        Token('IDENTIFIER', 'zinc_finger')])])
    """
    selection = "protein OR metal in zinc_finger"

    result = cfg.parser(selection)
    selection = result.children[0]
    selection_string, preposition, identifier = selection.children
    macro1, logic, macro2 = selection_string.children

    assert selection.data == "selection"
    assert selection_string.data == "selection_string"
    assert macro1 == "protein"
    assert logic == "OR"
    assert macro2 == "metal"
    assert preposition == "in"
    assert identifier == "zinc_finger"


def test_parser_selection_ligand():
    """
    >>> import rich
    >>> from deemian.engine.cfg import parser
    >>> selection = parser("ligand in zinc_finger")
    >>> rich.print(selection)
    start
    └── selection
        ├── selection_string
        │   └── ligand
        ├── in
        └── zinc_finger
    >>> selection
    Tree(Token('RULE', 'start'),
      [Tree(Token('RULE', 'selection'),
        [Tree(Token('RULE', 'selection_string'),
          [Token('MACRO', 'ligand')]),
        Token('PREPOSITION', 'in'),
        Token('IDENTIFIER', 'zinc_finger')])])
    """
    selection = "ligand in zinc_finger"

    result = cfg.parser(selection)
    selection = result.children[0]
    selection_string, preposition, identifier = selection.children
    macro = selection_string.children[0]

    assert selection.data == "selection"
    assert selection_string.data == "selection_string"
    assert macro == "ligand"
    assert preposition == "in"
    assert identifier == "zinc_finger"


def test_parser_measurement():
    """
    >>> import rich
    >>> from deemian.engine.cfg import parser
    >>> measurement = parser("nonpolar on protein")
    >>> rich.print(measurement)
    start
    └── measurement
        ├── interactions
        │   └── nonpolar
        └── target
            ├── on
            ├── protein
            ├── None
            └── None
    >>> measurement
    Tree(Token('RULE', 'start'),
      [Tree(Token('RULE', 'measurement'),
        [Tree(Token('RULE', 'interactions'),
          [Token('INTERACTION', 'nonpolar')]),
        Tree(Token('RULE', 'target'),
          [Token('PREPOSITION', 'on'),
          Token('IDENTIFIER', 'protein'),
          None,
          None])])])
    """
    measurement = "nonpolar on protein"

    result = cfg.parser(measurement)
    measurement = result.children[0]
    interaction, target = measurement.children

    assert measurement.data == "measurement"
    assert interaction.data == "interactions"
    assert interaction.children[0] == "nonpolar"
    assert target.data == "target"
    assert target.children[0] == "on"
    assert target.children[1] == "protein"
    assert target.children[2] is None
    assert target.children[3] is None


def test_parser_measurement_multi_interaction():
    """
    >>> import rich
    >>> from deemian.engine.cfg import parser
    >>> measurement = parser("nonpolar, electrostatic, hydrogen_bond on protein")
    >>> rich.print(measurement)
    start
    └── measurement
        ├── interactions
        │   ├── nonpolar
        │   ├── electrostatic
        │   └── hydrogen_bond
        └── target
            ├── on
            ├── protein
            ├── None
            └── None
    >>> measurement
    Tree(Token('RULE', 'start'),
      [Tree(Token('RULE', 'measurement'),
        [Tree(Token('RULE', 'interactions'),
          [Token('INTERACTION', 'nonpolar'),
          Token('INTERACTION', 'electrostatic'),
          Token('INTERACTION', 'hydrogen_bond')]),
        Tree(Token('RULE', 'target'),
          [Token('PREPOSITION', 'on'),
          Token('IDENTIFIER', 'protein'),
          None,
          None])])])
    """
    measurement = "nonpolar, electrostatic, hydrogen_bond on protein"

    result = cfg.parser(measurement)
    measurement = result.children[0]
    interaction, target = measurement.children

    assert interaction.children[0] == "nonpolar"
    assert interaction.children[1] == "electrostatic"
    assert interaction.children[2] == "hydrogen_bond"


def test_parser_measurement_multi_target():
    """
    >>> import rich
    >>> from deemian.engine.cfg import parser
    >>> measurement = parser("nonpolar between protein and ligand")
    >>> rich.print(measurement)
    start
    └── measurement
        ├── interactions
        │   └── nonpolar
        └── target
            ├── between
            ├── protein
            ├── and
            └── ligand
    >>> measurement
    Tree(Token('RULE', 'start'),
      [Tree(Token('RULE', 'measurement'),
        [Tree(Token('RULE', 'interactions'),
          [Token('INTERACTION', 'nonpolar')]),
        Tree(Token('RULE', 'target'),
          [Token('PREPOSITION', 'between'),
          Token('IDENTIFIER', 'protein'),
          Token('PREPOSITION', 'and'),
          Token('IDENTIFIER', 'ligand')])])])
    """
    measurement = "nonpolar between protein and ligand"

    result = cfg.parser(measurement)
    measurement = result.children[0]
    interaction, target = measurement.children

    assert target.children[0] == "between"
    assert target.children[1] == "protein"
    assert target.children[2] == "and"
    assert target.children[3] == "ligand"
