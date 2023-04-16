import pathlib

import pytest
import yaml


@pytest.fixture
def simple_yaml_parsed():
    """Parse simple YAML example and return as nested dictionary object"""

    this_file = pathlib.Path(__file__)
    this_dir = this_file.parent
    yaml_file = this_dir.joinpath("data", "simple-nonpolar.yaml")

    with open(yaml_file, "r") as f:
        yaml_loaded = yaml.safe_load(f)

    return yaml_loaded
