"""Compatibility wrapper for the installed ``neidspecmatch-fit`` command."""

from neidspecmatch.cli import fit_main


if __name__ == "__main__":
    raise SystemExit(fit_main())
