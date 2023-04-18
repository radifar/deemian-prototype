# Deemian

A Domain Specific Language for Deep Molecular Interaction Analysis. Warning this is still a prototype/MVP with minimum feature.

## Installation

First create a new environment using Conda, make sure it is using Python 3.10 as currently Deemian only tested on Python 3.10.

```bash
conda create -n deemian-prototype python=3.10
```

Then download the current snapshot of this repository in [zip format](https://github.com/radifar/deemian/archive/refs/heads/main.zip).
Extract this zip file in your local computer and enter the `deemian-main` directory, activate `deemian-prototype` environment, and install using pip.

```bash
cd deemian-main
conda activate deemian-prototype
pip install .
```

## Usage

- TODO

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`deemian` was created by Muhammad Radifar. It is licensed under the terms of the Apache License 2.0 license.

## Credits

`deemian` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).

The test data include the PDB data retrieved from RCSB.org under the
[CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/) license.
The PDBs being used are as follows

- [1ZNM](https://www.rcsb.org/structure/1znm): Viles, John H., et al. "Design, synthesis and structure of a zinc finger with an artificial β-turn." Journal of molecular biology 279.4 (1998): 973-986. <https://doi.org/10.1006/jmbi.1998.1764>
