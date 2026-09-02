"""Optional NEID Archive helpers with environment-driven credentials."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import stat
import tempfile

from .paths import safe_filename_component


USER_ENV = "NEID_ARCHIVE_USER"
PASSWORD_ENV = "NEID_ARCHIVE_PASSWORD"


def _ensure_private_directory(directory):
    directory = Path(directory).expanduser().absolute()
    if directory.is_symlink():
        raise ValueError("Archive session directory cannot be a symlink.")
    if directory.exists():
        if not directory.is_dir():
            raise ValueError("Archive session path must be a directory.")
    else:
        directory.mkdir(parents=True, mode=0o700)
    directory.chmod(0o700)
    if stat.S_IMODE(directory.stat().st_mode) != 0o700:
        raise PermissionError("Archive session directory must have mode 0700.")
    return directory


def _prepare_private_file(path):
    path = Path(path)
    if path.is_symlink():
        raise ValueError(f"Private archive file cannot be a symlink: {path.name}")
    if path.exists() and not path.is_file():
        raise ValueError(f"Private archive path is not a regular file: {path.name}")
    flags = os.O_WRONLY | os.O_CREAT
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    os.close(descriptor)
    path.chmod(0o600)
    if not stat.S_ISREG(path.stat().st_mode):
        raise ValueError(f"Private archive path is not regular: {path.name}")
    return path


def get_archive_credentials(user=None, password=None):
    """Return explicit/environment credentials or fail without echoing them."""
    user = user or os.environ.get(USER_ENV)
    password = password or os.environ.get(PASSWORD_ENV)
    missing = [
        name for name, value in ((USER_ENV, user), (PASSWORD_ENV, password))
        if not value
    ]
    if missing:
        raise RuntimeError(
            "Missing NEID Archive credentials. Set {} in the environment.".format(
                " and ".join(missing)
            )
        )
    return user, password


def login(directory, *, user=None, password=None, debug=False):
    """Log in and return the cookie path, using environment credentials."""
    try:
        from pyneid.neid import Neid
    except ImportError as exc:
        raise ImportError(
            "NEID Archive access requires optional dependencies: "
            "pip install neidspecmatch[archive]"
        ) from exc
    user, password = get_archive_credentials(user=user, password=password)
    directory = _ensure_private_directory(directory)
    cookie_path = _prepare_private_file(directory / "neid_archive_cookie.txt")
    debug_path = (
        _prepare_private_file(directory / "neid_archive.debug") if debug else None
    )
    kwargs = {"userid": user, "password": password, "cookiepath": str(cookie_path)}
    if debug_path is not None:
        kwargs["debugfile"] = str(debug_path)
    Neid.login(**kwargs)
    _prepare_private_file(cookie_path)
    if debug_path is not None:
        _prepare_private_file(debug_path)
    return cookie_path


def get_table_singlestar(
        starname, ra, dec, fmt="csv", download=False, directory=None,
        cookie_path=None, verbose=0):
    """Query NEID L2 products around one sky position."""
    try:
        from pyneid.neid import Neid
    except ImportError as exc:
        raise ImportError(
            "NEID Archive access requires optional dependencies: "
            "pip install neidspecmatch[archive]"
        ) from exc
    if fmt not in {"csv", "ipac"}:
        raise ValueError("fmt must be 'csv' or 'ipac'.")
    directory = Path(directory or ".").expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    cookie_path = Path(cookie_path or directory / "neid_archive_cookie.txt")
    suffix = "csv" if fmt == "csv" else "tbl"
    component = safe_filename_component(starname)
    outpath = directory / f"{component}_L2_q.{suffix}"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{component}-", suffix=f".{suffix}.tmp", dir=directory
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    if verbose:
        print(f"Querying {starname}")
    try:
        Neid.query_position(
            "l2", f"circle {ra} {dec} 0.5", cookiepath=str(cookie_path),
            format=fmt, outpath=str(temporary_path),
        )
        os.chmod(temporary_path, 0o600)
        os.replace(temporary_path, outpath)
    except BaseException:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass
        raise
    if fmt == "ipac" and download:
        Neid.download(
            str(outpath), "l2", "ipac", str(directory),
            cookiepath=str(cookie_path),
        )
    return outpath


def _parser():
    parser = argparse.ArgumentParser(description="Log in to the NEID Archive safely.")
    parser.add_argument("directory", help="Directory for the private session cookie")
    parser.add_argument("--debug", action="store_true", help="Write a pyneid debug log")
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    cookie_path = login(args.directory, debug=args.debug)
    print(f"Authenticated; private cookie written to {cookie_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
