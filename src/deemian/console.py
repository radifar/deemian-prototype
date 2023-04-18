import click

from . import __version__


@click.command()
@click.version_option(version=__version__)
def main():
    """Deemian - A Domain Specific Language for Deep Molecular Interaction Analysis"""
    click.echo("This is just a prototype!")
