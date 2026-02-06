import csv
import os
import pickle
import platform
import subprocess
from glob import glob

from astropy.io import fits
import numpy as np

import neidspecmatch


def _format_absrv_suffix(value):
    """Return a filesystem-friendly suffix describing an absrv value."""
    if value is None:
        return "auto"
    sign = "p" if value >= 0 else "m"
    magnitude = f"{abs(value):06.1f}"
    return f"{sign}{magnitude}"


def _resolve_absrv_values(absrv, absrv_grid):
    """
    Determine which absrv values to use for SpecMatch runs.

    If absrv_grid is None, return a single-element list containing absrv.
    If absrv_grid is a (start, stop, step) tuple, build an inclusive sequence.
    Otherwise, treat absrv_grid as an iterable of explicit values.
    """
    if absrv_grid is None:
        return [absrv]

    values = []
    if isinstance(absrv_grid, tuple) and len(absrv_grid) == 3:
        start, stop, step = absrv_grid
        start = float(start)
        stop = float(stop)
        step = float(step)
        if step == 0:
            raise ValueError("absrv_grid step cannot be zero.")
        if step > 0 and start > stop:
            raise ValueError("For a positive step, start must be <= stop.")
        if step < 0 and start < stop:
            raise ValueError("For a negative step, start must be >= stop.")

        current = start
        if step > 0:
            while current <= stop or np.isclose(current, stop):
                values.append(current)
                current += step
        else:
            while current >= stop or np.isclose(current, stop):
                values.append(current)
                current += step

        if not np.isclose(values[-1], stop):
            values.append(stop)
    else:
        try:
            cleaned = []
            for val in absrv_grid:
                cleaned.append(float(val) if val is not None else None)
            values = cleaned
        except TypeError as exc:
            raise ValueError("absrv_grid must be iterable or a (start, stop, step) tuple.") from exc

    if not values:
        raise ValueError("absrv_grid produced no values; adjust the range or iterable.")

    return values


def process_neid_fits(
    input_dir,
    output_dir,
    filename,
    orders=('55', '101', '102', '103'),
    maxvsini=150,
    calibrate_feh=True,
    scaleres=1,
    deblazed=False,
    mode='HR',
    save_plot_data=True,
    absrv=None,
    absrv_grid=(-150, 150, 30),
    vsini=28,
    refine_absrv=True,
    max_refinement_iterations=3,
    notify_on_complete=False,
    notification_sound="Submarine"
):
    """
    Process all NEID FITS files in input_dir and run SpecMatch, saving output to output_dir.
    Each output will be in a subfolder named after the FITS filename (without extension).
    By default, when absrv is unknown (absrv=None), a grid of absolute RVs spanning
    -150 to +150 km/s in 10 km/s steps is sampled. Grid runs are consolidated into a
    single {basename} folder (with per-absrv order folders tagged by the absrv suffix)
    and only order 102 is evaluated for speed. Provide a numeric absrv to run a single
    fit with all requested orders, or adjust absrv_grid (tuple or iterable) to control
    the grid.

    Set notify_on_complete=True to send a desktop notification when processing finishes
    (macOS via osascript, Linux via notify-send). Set notification_sound to None to
    silence the chime, or provide a different macOS sound name / Linux sound-name hint.
    """
    path_df_lib = neidspecmatch.config.PATH_LIBRARY_DB
    path_df_lib_fits = neidspecmatch.config.PATH_LIBRARY_FITS

    # Ensure spectral library is downloaded
    neidspecmatch.get_library()

    # Find all FITS files
    fits_files = glob(os.path.join(input_dir, filename))

    if not fits_files:
        print(f"No FITS files found in {input_dir}")
        return
    absrv_values = _resolve_absrv_values(absrv, absrv_grid)
    grid_mode = absrv_grid is not None
    orders_to_use = tuple(str(o) for o in orders)

    if grid_mode:
        orders_to_use = ('102',)
        print("absrv grid mode detected: restricting runs to order 102 and consolidating outputs.")

    if grid_mode:
        grid_desc = ", ".join("auto" if v is None else f"{v:.1f}" for v in absrv_values)
        print(f"Running absrv grid with values: {grid_desc}")

    for filepath in fits_files:
        try:
            with fits.open(filepath) as hdul:
                targetname = hdul[0].header['OBJECT']
        except Exception as e:
            print(f"Skipping {filepath} due to error: {e}")
            continue

        base_name = os.path.splitext(os.path.basename(filepath))[0]
        base_outputdir = os.path.join(output_dir, base_name)
        if grid_mode:
            os.makedirs(base_outputdir, exist_ok=True)

        for rv_value in absrv_values:
            if grid_mode:
                outputdir = base_outputdir
                run_label = f"absrv_{_format_absrv_suffix(rv_value)}"
            else:
                run_label = None
                outputdir = base_outputdir
            os.makedirs(outputdir, exist_ok=True)
            absrv_label = "auto" if rv_value is None else f"{rv_value:.1f}"
            print(f"Processing {filepath} with absrv={absrv_label} into {outputdir}")

            neidspecmatch.run_specmatch_for_orders(
                targetfile=filepath,
                targetname=targetname,
                outputdirectory=outputdir,
                path_df_lib=path_df_lib,
                path_df_lib_fits=path_df_lib_fits,
                orders=orders_to_use,
                maxvsini=maxvsini,
                calibrate_feh=calibrate_feh,
                scaleres=scaleres,
                deblazed=deblazed,
                mode=mode,
                save_plot_data=save_plot_data,
                absrv=rv_value,
                vsini=vsini,
                refine_absrv=refine_absrv,
                max_refinement_iterations=max_refinement_iterations,
                run_label=run_label
            )

    if notify_on_complete:
        message = f"Finished processing {len(fits_files)} file(s) into {output_dir}."
        _send_notification(message, title="NEID SpecMatch", sound=notification_sound)


def _send_notification(message, title="NEID SpecMatch", sound=None):
    """Fire a best-effort desktop notification."""
    system = platform.system()
    try:
        if system == "Darwin":
            safe_message = message.replace('"', '\\"')
            safe_title = title.replace('"', '\\"')
            script = f'display notification "{safe_message}" with title "{safe_title}"'
            if sound:
                safe_sound = sound.replace('"', '\\"')
                script += f' sound name "{safe_sound}"'
            subprocess.run(["osascript", "-e", script], check=False)
        elif system == "Linux":
            args = ["notify-send"]
            if sound:
                args.extend(["--hint", f"string:sound-name:{sound}"])
            args.extend([title, message])
            subprocess.run(args, check=False)
            if sound and os.path.exists(sound):
                subprocess.run(["paplay", sound], check=False)
    except Exception as exc:
        print(f"Notification failed: {exc}")


def gather_neidl2_pickle_results(
    base_dir,
    output_csv="results_summary.csv",
    prefix="neidL2",
    verbose=True,
    refined_only=True
):
    """
    Traverse the given directory, extract .pkl results, and save them to a CSV.

    Each row = one neidL2 folder.
    Each column = teff/logg/feh/vsini per order.
    If refined_only=True, only include order folders ending in "_refined".
    """
    rows = []
    all_orders = set()

    # First pass: collect all unique order numbers
    for main_folder in sorted(os.listdir(base_dir)):
        main_path = os.path.join(base_dir, main_folder)
        if not os.path.isdir(main_path) or not main_folder.startswith(prefix):
            continue
        for subfolder in sorted(os.listdir(main_path)):
            subfolder_path = os.path.join(main_path, subfolder)
            if os.path.isdir(subfolder_path) and "_" in subfolder:
                parts = subfolder.split("_")
                if refined_only:
                    if len(parts) < 3 or parts[-1] != "refined":
                        continue
                    order = parts[-2]
                else:
                    if parts[-1] == "refined":
                        if len(parts) < 3:
                            continue
                        order = parts[-2]
                    else:
                        order = parts[-1]
                all_orders.add(order)

    all_orders = sorted(all_orders, key=lambda x: int(x) if x.isdigit() else float('inf'))

    # Build headers
    headers = ["Folder", "TIC"]
    for order in all_orders:
        headers += [f"teff_{order}", f"logg_{order}", f"feh_{order}", f"vsini_{order}"]

    # Second pass: extract data
    for main_folder in sorted(os.listdir(base_dir)):
        main_path = os.path.join(base_dir, main_folder)
        if not os.path.isdir(main_path) or not main_folder.startswith(prefix):
            continue

        row_data = {"Folder": main_folder, "TIC": ""}
        for order in all_orders:
            row_data[f"teff_{order}"] = ""
            row_data[f"logg_{order}"] = ""
            row_data[f"feh_{order}"] = ""
            row_data[f"vsini_{order}"] = ""

        any_subfolder_found = False

        for subfolder in sorted(os.listdir(main_path)):
            subfolder_path = os.path.join(main_path, subfolder)
            if os.path.isdir(subfolder_path) and "_" in subfolder:
                any_subfolder_found = True
                parts = subfolder.split("_")
                if refined_only:
                    if len(parts) < 3 or parts[-1] != "refined":
                        continue
                    order = parts[-2]
                    tic_name = "_".join(parts[:-2])
                else:
                    if len(parts) < 2:
                        continue
                    if parts[-1] == "refined":
                        if len(parts) < 3:
                            continue
                        order = parts[-2]
                        tic_name = "_".join(parts[:-2])
                    else:
                        order = parts[-1]
                        tic_name = "_".join(parts[:-1])

                if not row_data["TIC"]:
                    row_data["TIC"] = tic_name

                pkl_file = os.path.join(subfolder_path, f"{tic_name}_results.pkl")
                if os.path.exists(pkl_file):
                    try:
                        with open(pkl_file, "rb") as f:
                            data = pickle.load(f)
                        row_data[f"teff_{order}"] = data.get("teff", "")
                        row_data[f"logg_{order}"] = data.get("logg", "")
                        row_data[f"feh_{order}"] = data.get("feh", "")
                        row_data[f"vsini_{order}"] = data.get("vsini", "")
                        if verbose:
                            print(f"Loaded: {pkl_file}")
                    except Exception as e:
                        if verbose:
                            print(f"Error reading {pkl_file}: {e}")

        # Always keep the row even if empty
        rows.append([row_data.get(col, "") for col in headers])

    # Write CSV
    output_path = os.path.join(base_dir, output_csv)
    with open(output_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(rows)

    if verbose:
        print(f"\nDone! Results saved to {output_path}")

if __name__ == '__main__':
    # Set input directory containing all FITS files
    input_dir = '/Users/tehan/Documents/SURFSUP/NEID/NEID_Spectra'  # or wherever your .fits files are
    output_dir = '/Users/tehan/Documents/SURFSUP/NEID/NEIDSM_BF'
    filename = 'neidL2_20230905T090219.fits'
    vsini = 44.85
    # Example: sample absrv grid for manual inspection (set refine_absrv=False for speed)
    # process_neid_fits(
    #     input_dir,
    #     output_dir,
    #     filename=filename,
    #     orders=('55', '101', '102', '103'),
    #     maxvsini=150,
    #     calibrate_feh=True,
    #     scaleres=1,
    #     deblazed=False,
    #     mode='HR',
    #     save_plot_data=True,
    #     absrv=None,
    #     absrv_grid=(-150, 150, 30),
    #     vsini=vsini,
    #     refine_absrv=False,
    #     max_refinement_iterations=3,
    #     notify_on_complete=True,
    #     notification_sound="Submarine"
    # )
    # Example: sample absrv grid for manual inspection (set refine_absrv=False for speed)
    # process_neid_fits(
    #     input_dir,
    #     output_dir,
    #     filename=filename,
    #     orders=('55', '101', '102', '103'), # H alpha '80'
    #     maxvsini=150,
    #     calibrate_feh=True,
    #     scaleres=1,
    #     deblazed=False,
    #     mode='HR',
    #     save_plot_data=True,
    #     absrv=90,
    #     absrv_grid=None,
    #     vsini=vsini,
    #     refine_absrv=True,
    #     max_refinement_iterations=3,
    #     notify_on_complete=True,
    #     notification_sound="Submarine"
    # )

    gather_neidl2_pickle_results(output_dir)
