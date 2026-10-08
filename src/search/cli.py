"""Command line interface. Commands are added phase by phase."""

from __future__ import annotations

import typer

from search import __version__

app = typer.Typer(no_args_is_help=True, help="Seed-based gene co-expression clustering.")


def _print_version(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version", callback=_print_version, is_eager=True, help="Print the version."
    ),
) -> None:
    """Search: seed-based gene co-expression clustering with named presets."""


@app.command()
def version() -> None:
    """Print the installed version."""
    typer.echo(__version__)


if __name__ == "__main__":  # pragma: no cover
    app()
