"""Compatibility wrapper for the installed ``neidspecmatch-crossval`` command."""

from neidspecmatch.cli import crossval_main


if __name__ == "__main__":
    raise SystemExit(crossval_main())
