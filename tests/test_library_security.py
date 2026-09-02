import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock
import zipfile
import hashlib
import stat
import types
import numpy as np
from astropy.io import fits

from neidspecmatch import config
from neidspecmatch import utils
from neidspecmatch import neid_archive
import neidspecmatch.neidspecmatch as core


class LibraryAndSecurityTests(unittest.TestCase):
    def make_library(self, root, count=78, manifest=False):
        root = Path(root)
        (root / "FITS").mkdir(parents=True)
        (root / config.DEFAULT_LIBRARY_CATALOG).write_text("OBJECT_ID\n", encoding="utf-8")
        for index in range(count):
            (root / "FITS" / f"star-{index:02d}.fits").write_bytes(b"")
        if manifest:
            files = [root / config.DEFAULT_LIBRARY_CATALOG]
            files.extend(sorted((root / "FITS").glob("*.fits")))
            (root / "library_manifest.json").write_text(json.dumps({
                "schema_version": utils.LIBRARY_MANIFEST_SCHEMA,
                "library_id": config.DEFAULT_LIBRARY_ID,
                "catalog": config.DEFAULT_LIBRARY_CATALOG,
                "fits_count": count,
                "source_archive_integrity": {
                    "algorithm": "md5",
                    "value": config.LIBRARY_ZIP_MD5,
                    "scope": "integrity_only_not_authentication",
                },
                "files": [
                    {
                        "path": filename.relative_to(root).as_posix(),
                        "size_bytes": filename.stat().st_size,
                        "sha256": hashlib.sha256(filename.read_bytes()).hexdigest(),
                    }
                    for filename in files
                ],
            }), encoding="utf-8")
        return root

    def test_validate_library_layout_and_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            library = self.make_library(Path(temp) / "library", manifest=True)
            result = utils.validate_library(library, require_manifest=True)
            self.assertEqual(result["fits_count"], 78)
            (library / "FITS" / "star-00.fits").unlink()
            with self.assertRaises(ValueError):
                utils.validate_library(library)

    def test_safe_zip_extraction_rejects_path_traversal(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            archive = temp / "bad.zip"
            destination = temp / "destination"
            destination.mkdir()
            with zipfile.ZipFile(archive, "w") as stream:
                stream.writestr("../escaped.txt", "bad")
            with self.assertRaises(ValueError):
                utils._safe_extract_zip(archive, destination)
            self.assertFalse((temp / "escaped.txt").exists())

    def test_safe_zip_rejects_symlink_member_and_resource_limits(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            archive = temp / "bad.zip"
            destination = temp / "destination"
            info = zipfile.ZipInfo("link")
            info.create_system = 3
            info.external_attr = (stat.S_IFLNK | 0o777) << 16
            with zipfile.ZipFile(archive, "w") as stream:
                stream.writestr(info, "target")
            with self.assertRaises(ValueError):
                utils._safe_extract_zip(archive, destination)

            archive = temp / "many.zip"
            with zipfile.ZipFile(archive, "w") as stream:
                stream.writestr("one", "1")
                stream.writestr("two", "2")
            with self.assertRaises(ValueError):
                utils._safe_extract_zip(
                    archive, temp / "many", max_members=1
                )
            with self.assertRaises(ValueError):
                utils._safe_extract_zip(
                    archive, temp / "large", max_uncompressed_bytes=1
                )

    def test_download_timeout_length_and_size_guards(self):
        class Response:
            def __init__(self, body, length):
                self.body = body
                self.headers = {"Content-Length": str(length)}

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return None

            def read(self, size):
                chunk, self.body = self.body[:size], self.body[size:]
                return chunk

        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "download"
            with mock.patch.object(
                    utils, "urlopen", return_value=Response(b"abc", 4)) as opened:
                with self.assertRaises(ValueError):
                    utils._stream_download("https://example.invalid", output,
                                           timeout=7, max_bytes=10)
            opened.assert_called_once_with("https://example.invalid", timeout=7)
            too_large = Path(temp) / "too-large"
            with mock.patch.object(
                    utils, "urlopen", return_value=Response(b"abc", 30)):
                with self.assertRaises(ValueError):
                    utils._stream_download(
                        "https://example.invalid", too_large, max_bytes=10
                    )
            wrong_pinned_size = Path(temp) / "wrong-pinned-size"
            with mock.patch.object(
                    utils, "urlopen", return_value=Response(b"abc", 3)):
                with self.assertRaises(ValueError):
                    utils._stream_download(
                        "https://example.invalid", wrong_pinned_size,
                        expected_bytes=4,
                    )

    def test_atomic_replace_restores_prior_library_on_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            source = temp / "source"
            destination = temp / "library"
            source.mkdir()
            destination.mkdir()
            (source / "new").write_text("new", encoding="utf-8")
            (destination / "old").write_text("old", encoding="utf-8")
            real_replace = os.replace
            calls = 0

            def fail_second(src, dst):
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise OSError("simulated install failure")
                return real_replace(src, dst)

            with mock.patch.object(utils.os, "replace", side_effect=fail_second):
                with self.assertRaises(OSError):
                    utils._atomic_replace_directory(source, destination)
            self.assertEqual((destination / "old").read_text(), "old")
            self.assertFalse((destination / "new").exists())

    def test_manifest_sizes_and_deep_hashes_detect_corruption(self):
        with tempfile.TemporaryDirectory() as temp:
            library = self.make_library(Path(temp) / "library", manifest=True)
            catalog = library / config.DEFAULT_LIBRARY_CATALOG
            original = catalog.read_bytes()
            catalog.write_bytes(bytes([original[0] ^ 1]) + original[1:])
            utils.validate_library(library, require_manifest=True, deep=False)
            with self.assertRaises(ValueError):
                utils.validate_library(library, require_manifest=True, deep=True)

    def test_release_owned_allowlist_is_complete_and_enforced(self):
        packaged, digest = utils.load_default_library_allowlist()
        self.assertEqual(packaged["library_id"], config.DEFAULT_LIBRARY_ID)
        self.assertEqual(packaged["catalog"], config.DEFAULT_LIBRARY_CATALOG)
        self.assertEqual(packaged["fits_count"], 78)
        self.assertRegex(digest, r"^[0-9a-f]{64}$")
        paths = [record["path"] for record in packaged["files"]]
        self.assertEqual(len(paths), 79)
        self.assertEqual(len(paths), len(set(paths)))
        self.assertIsInstance(
            utils._validated_header_identity_exceptions(
                packaged.get("header_identity_exceptions", [])
            ),
            dict,
        )

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "library"
            (root / "FITS").mkdir(parents=True)
            catalog = root / "catalog.csv"
            spectrum = root / "FITS" / "star.fits"
            catalog.write_bytes(b"catalog")
            spectrum.write_bytes(b"spectrum")
            allowlist = {
                "allowlist_id": "test-release",
                "catalog": "catalog.csv",
                "fits_count": 1,
                "files": [
                    {
                        "path": filename.relative_to(root).as_posix(),
                        "size_bytes": filename.stat().st_size,
                        "sha256": hashlib.sha256(filename.read_bytes()).hexdigest(),
                    }
                    for filename in (catalog, spectrum)
                ],
            }
            receipt = utils._verify_library_against_allowlist(root, allowlist)
            self.assertTrue(receipt["allowlist_sha256_verified"])
            spectrum.write_bytes(b"spectruN")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                utils._verify_library_against_allowlist(root, allowlist)

    def test_download_manifest_preserves_release_identity_exceptions(self):
        exceptions = [{
            "path": "FITS/star.fits",
            "catalog_object_id": "star",
            "fits_object": "archive typo",
            "gaia_source_id": "123456789012",
            "reason": "synthetic test exception",
        }]
        allowlist = {
            "allowlist_id": "test-release",
            "header_identity_exceptions": exceptions,
        }

        with tempfile.TemporaryDirectory() as temp:
            destination = Path(temp) / "library"

            def fake_download(url, filename, **kwargs):
                with zipfile.ZipFile(filename, "w") as stream:
                    stream.writestr(config.DEFAULT_LIBRARY_CATALOG, "OBJECT_ID\n")
                    stream.writestr("FITS/star.fits", b"")
                return config.LIBRARY_ZIP_MD5

            with (
                mock.patch.object(
                    utils, "load_default_library_allowlist",
                    return_value=(allowlist, "a" * 64),
                ),
                mock.patch.object(
                    utils, "_stream_download", side_effect=fake_download
                ),
                mock.patch.object(
                    utils, "_verify_library_against_allowlist"
                ),
                mock.patch.object(utils, "validate_library", return_value={}),
                mock.patch.object(utils, "build_library_manifest") as build,
            ):
                utils.get_library(
                    library_path=destination, overwrite=True, verbose=0
                )

        self.assertEqual(
            build.call_args.kwargs["header_identity_exceptions"], exceptions
        )

    def test_release_allowlist_rejects_extras_with_bounded_diagnostics(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "library"
            (root / "FITS").mkdir(parents=True)
            catalog = root / "catalog.csv"
            spectrum = root / "FITS" / "star.fits"
            catalog.write_bytes(b"catalog")
            spectrum.write_bytes(b"spectrum")
            allowlist = {
                "allowlist_id": "test-release",
                "catalog": "catalog.csv",
                "fits_count": 1,
                "files": [
                    {
                        "path": filename.relative_to(root).as_posix(),
                        "size_bytes": filename.stat().st_size,
                        "sha256": hashlib.sha256(filename.read_bytes()).hexdigest(),
                    }
                    for filename in (catalog, spectrum)
                ],
            }
            for index in range(30):
                (root / f"extra-{index:02d}").write_bytes(b"extra")
            with self.assertRaises(ValueError) as caught:
                utils._verify_library_against_allowlist(root, allowlist)
            message = str(caught.exception)
            self.assertIn("extra_count=30", message)
            self.assertLess(len(message), 1000)

    def test_science_receipt_deep_verifies_and_binds_loaded_references(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            library = self.make_library(temp / "library", manifest=True)
            target_file = temp / "target.fits"
            target_file.write_bytes(b"target")

            class Spectrum:
                def __init__(self, filename):
                    self.filename = str(filename)
                    self.f = np.ones(4)
                    self.f_sci = np.ones(4)
                    self.f_sky = np.ones(4)
                    self.e_sci = np.ones(4)
                    self.e_sky = np.ones(4)
                    self.error_scale = np.ones(4)
                    self.blaze_source = "l2"
                    self._neidspecmatch_loaded_state_sha256 = (
                        core._spectrum_state_fingerprint(self)
                    )

            target = Spectrum(target_file)
            references = [
                Spectrum(path) for path in sorted((library / "FITS").glob("*.fits"))
            ]
            provenance = core._library_manifest_provenance(
                library / config.DEFAULT_LIBRARY_CATALOG, deep=True
            )
            provenance.update(core._verify_reference_pool_against_manifest(
                target, references, library
            ))
            self.assertTrue(provenance["deep_verified"])
            self.assertTrue(provenance["reference_pool_verified"])

            references[0].f[0] = 2.0
            with self.assertRaises(ValueError):
                core._verify_reference_pool_against_manifest(
                    target, references, library
                )
            references = [
                Spectrum(path) for path in sorted((library / "FITS").glob("*.fits"))
            ]

            for field in (
                "f_sci", "f_sky", "e_sci", "e_sky", "error_scale"
            ):
                field_references = [
                    Spectrum(path)
                    for path in sorted((library / "FITS").glob("*.fits"))
                ]
                getattr(field_references[0], field)[0] = 2.0
                with self.assertRaises(ValueError, msg=field):
                    core._verify_reference_pool_against_manifest(
                        target, field_references, library
                    )

            unrelated_dir = temp / "unrelated"
            unrelated_dir.mkdir()
            unrelated = unrelated_dir / "star-00.fits"
            unrelated.write_bytes(b"different")
            references[0] = Spectrum(unrelated)
            with self.assertRaises(ValueError):
                core._verify_reference_pool_against_manifest(
                    target, references, library
                )

    def test_custom_library_manifest_and_symlink_rejection(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "custom"
            (root / "FITS").mkdir(parents=True)
            (root / "catalog-v2.csv").write_text(
                "OBJECT_ID,basenames\na,a.fits\nb,b.fits\n", encoding="utf-8"
            )
            (root / "FITS" / "a.fits").write_bytes(b"a")
            (root / "FITS" / "b.fits").write_bytes(b"b")
            manifest = utils.build_library_manifest(
                root, library_id="neid-drp-1.5-library",
                catalog_name="catalog-v2.csv", expected_fits=2,
            )
            self.assertTrue(manifest.is_file())
            result = utils.validate_library(
                root, require_manifest=True, deep=True,
                expected_library_id="neid-drp-1.5-library",
                catalog_name="catalog-v2.csv", expected_fits=2,
            )
            self.assertEqual(result["library_id"], "neid-drp-1.5-library")
            with self.assertRaises(ValueError):
                utils.validate_library(
                    root, require_manifest=True,
                    expected_library_id="wrong-library",
                    catalog_name="catalog-v2.csv",
                )
            manifest.unlink()
            manifest.symlink_to(root / "catalog-v2.csv")
            with self.assertRaises(ValueError):
                utils.validate_library(root, require_manifest=True)

    def test_checksum_failure_leaves_no_partial_install(self):
        with tempfile.TemporaryDirectory() as temp:
            destination = Path(temp) / "library"

            def fake_download(url, filename, **kwargs):
                Path(filename).write_bytes(b"not a zip")
                return "0" * 32

            with mock.patch.object(utils, "_stream_download", side_effect=fake_download):
                with self.assertRaises(ValueError):
                    utils.get_library(library_path=destination, verbose=0)
            self.assertFalse(destination.exists())

    def test_credentials_are_environment_driven_and_not_echoed(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(RuntimeError) as caught:
                neid_archive.get_archive_credentials()
        self.assertNotIn("password=", str(caught.exception).lower())
        with mock.patch.dict(os.environ, {
                neid_archive.USER_ENV: "archive-user",
                neid_archive.PASSWORD_ENV: "private-value",
        }, clear=True):
            user, password = neid_archive.get_archive_credentials()
        self.assertEqual(user, "archive-user")
        self.assertEqual(password, "private-value")

    def test_archive_login_enforces_private_regular_paths(self):
        calls = []

        class FakeNeid:
            @staticmethod
            def login(**kwargs):
                calls.append(kwargs)
                Path(kwargs["cookiepath"]).write_text("cookie", encoding="utf-8")
                if "debugfile" in kwargs:
                    Path(kwargs["debugfile"]).write_text("debug", encoding="utf-8")

        module = types.ModuleType("pyneid.neid")
        module.Neid = FakeNeid
        package = types.ModuleType("pyneid")
        package.neid = module
        with tempfile.TemporaryDirectory() as temp, mock.patch.dict(
            sys.modules, {"pyneid": package, "pyneid.neid": module}
        ):
            session = Path(temp) / "session"
            synthetic_password = "synthetic-" + "credential"
            cookie = neid_archive.login(
                session, user="archive-user", password=synthetic_password,
                debug=True,
            )
            self.assertEqual(stat.S_IMODE(session.stat().st_mode), 0o700)
            self.assertEqual(stat.S_IMODE(cookie.stat().st_mode), 0o600)
            self.assertEqual(
                stat.S_IMODE((session / "neid_archive.debug").stat().st_mode),
                0o600,
            )
            self.assertEqual(len(calls), 1)

            symlink_session = Path(temp) / "session-link"
            symlink_session.symlink_to(session, target_is_directory=True)
            with self.assertRaises(ValueError):
                neid_archive.login(
                    symlink_session, user="archive-user",
                    password=synthetic_password,
                )
            cookie.unlink()
            cookie.symlink_to(Path(temp) / "elsewhere")
            with self.assertRaises(ValueError):
                neid_archive.login(
                    session, user="archive-user", password=synthetic_password
                )

    def test_package_import_is_quiet(self):
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
        environment["MPLCONFIGDIR"] = tempfile.mkdtemp(prefix="neidsm-mpl-")
        completed = subprocess.run(
            [sys.executable, "-c", (
                "import builtins, matplotlib; "
                "before=matplotlib.rcParams['axes.linewidth']; "
                "original=builtins.__import__; "
                "blocked={'pyneid','astroquery','barycorrpy','emcee','radvel'}; "
                "builtins.__import__=lambda name,*a,**k: "
                "(_ for _ in ()).throw(ImportError(name)) if name.split('.')[0] in blocked "
                "else original(name,*a,**k); "
                "import neidspecmatch; "
                "assert matplotlib.rcParams['axes.linewidth']==before"
            )],
            check=True, capture_output=True, text=True, env=environment,
        )
        self.assertEqual(completed.stdout, "")

    def test_dedicated_fits_schema_and_object_mapping_validation(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "library"
            (root / "FITS").mkdir(parents=True)
            (root / config.DEFAULT_LIBRARY_CATALOG).write_text(
                "OBJECT_ID,basenames\nstar,star.fits\n", encoding="utf-8"
            )
            shape = (1, 4)
            header = fits.Header({
                "OBJECT": "star", "INSTRUME": "NEID", "OBSTYPE": "SCI",
                "DATALVL": 2, "OBS-MODE": "HR", "E_VER": "1.3.0",
                "DQLEVEL1": 0, "DQLEVEL2": 0,
            })
            hdus = [fits.PrimaryHDU(header=header)]
            for index in range(1, 18):
                hdus.append(fits.ImageHDU(data=np.ones(shape) if index in {
                    1, 2, 4, 5, 7, 8, 15, 16, 17
                } else np.zeros(shape)))
            fits.HDUList(hdus).writeto(root / "FITS" / "star.fits")
            result = utils.validate_library(
                root, expected_fits=1, validate_fits_schema=True
            )
            self.assertEqual(result["fits_schema"]["fits_count"], 1)
            self.assertEqual(
                result["fits_schema"]["mapped_objects"][0]["OBJECT_ID"],
                "star",
            )

    def test_header_identity_exception_requires_matching_gaia_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "library"
            (root / "FITS").mkdir(parents=True)
            catalog = root / "catalog.csv"
            catalog.write_text(
                "OBJECT,OBJECT_ID,basenames,simbadnames_x\n"
                "Gl_3378,GJ 3378,star.fits,"
                "Gaia DR3 999747822484300800|GJ 3378|\n",
                encoding="utf-8",
            )
            shape = (1, 4)
            header = fits.Header({
                "OBJECT": "GJ3778",
                "QOBJECT": "Gaia DR2 999747822484300800",
                "SCI-OBJ": "Gaia DR2 999747822484300800",
                "INSTRUME": "NEID",
                "OBSTYPE": "SCI",
                "DATALVL": 2,
                "OBS-MODE": "HR",
                "E_VER": "v1.5.3",
                "DQLEVEL1": 0,
                "DQLEVEL2": 0,
            })
            hdus = [fits.PrimaryHDU(header=header)]
            for index in range(1, 18):
                hdus.append(fits.ImageHDU(data=np.ones(shape) if index in {
                    1, 2, 4, 5, 7, 8, 15, 16, 17
                } else np.zeros(shape)))
            fits.HDUList(hdus).writeto(root / "FITS" / "star.fits")

            with self.assertRaisesRegex(ValueError, "no declared exception"):
                utils.validate_library_fits_schema(
                    root, catalog_name="catalog.csv"
                )

            exception = {
                "path": "FITS/star.fits",
                "catalog_object_id": "GJ 3378",
                "fits_object": "GJ3778",
                "gaia_source_id": "999747822484300800",
                "reason": "Archive OBJECT keyword omits one digit.",
            }
            result = utils.validate_library_fits_schema(
                root,
                catalog_name="catalog.csv",
                header_identity_exceptions=[exception],
            )
            self.assertEqual(
                result["mapped_objects"][0]["identity_status"],
                "manifest_exception_gaia_verified",
            )

            bad_exception = dict(exception, gaia_source_id="123456789012345678")
            with self.assertRaisesRegex(ValueError, "Gaia provenance"):
                utils.validate_library_fits_schema(
                    root,
                    catalog_name="catalog.csv",
                    header_identity_exceptions=[bad_exception],
                )

    def test_removed_fsr_generator_has_no_broken_runtime_path(self):
        self.assertFalse(hasattr(config, "PATH_FSR"))


if __name__ == "__main__":
    unittest.main()
