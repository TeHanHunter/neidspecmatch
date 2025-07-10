import neidspecmatch
import numpy as np
import os
from astropy.io import fits
from glob import glob
import pickle
import csv


def process_neid_fits(
    input_dir,
    output_dir,
    orders=('55', '101', '102', '103'),
    maxvsini=150,
    calibrate_feh=True,
    scaleres=1,
    deblazed=False,
    mode='HR',
    save_plot_data=True
):
    """
    Process all NEID FITS files in input_dir and run SpecMatch, saving output to output_dir.
    Each output will be in a subfolder named after the FITS filename (without extension).
    """
    path_df_lib = neidspecmatch.config.PATH_LIBRARY_DB
    path_df_lib_fits = neidspecmatch.config.PATH_LIBRARY_FITS

    # Ensure spectral library is downloaded
    neidspecmatch.get_library()

    # Find all FITS files
    fits_files = glob(os.path.join(input_dir, '*.fits'))

    if not fits_files:
        print(f"No FITS files found in {input_dir}")
        return

    for filepath in fits_files:
        try:
            with fits.open(filepath) as hdul:
                targetname = hdul[0].header['OBJECT']
        except Exception as e:
            print(f"Skipping {filepath} due to error: {e}")
            continue

        outputdir = os.path.join(output_dir, os.path.splitext(os.path.basename(filepath))[0])
        os.makedirs(outputdir, exist_ok=True)

        print(f"Processing {filepath} into {outputdir}")

        neidspecmatch.run_specmatch_for_orders(
            targetfile=filepath,
            targetname=targetname,
            outputdirectory=outputdir,
            path_df_lib=path_df_lib,
            path_df_lib_fits=path_df_lib_fits,
            orders=orders,
            maxvsini=maxvsini,
            calibrate_feh=calibrate_feh,
            scaleres=scaleres,
            deblazed=deblazed,
            mode=mode,
            save_plot_data=save_plot_data
        )


def gather_neidl2_pickle_results(
    base_dir,
    output_csv="results_summary.csv",
    prefix="neidL2",
    verbose=True
):
    """
    Traverse the given directory, extract .pkl results, and save them to a CSV.

    Each row = one neidL2 folder.
    Each column = teff/logg/feh/vsini per order.
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
                order = subfolder.split("_")[-1]
                all_orders.add(order)

    all_orders = sorted(all_orders, key=lambda x: int(x) if x.isdigit() else x)

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
                if len(parts) >= 2:
                    tic_name = "_".join(parts[:-1])
                    if not row_data["TIC"]:
                        row_data["TIC"] = tic_name
                    order = parts[-1]
                    pkl_prefix = "_".join(parts[:-1])
                    pkl_file = os.path.join(subfolder_path, f"{pkl_prefix}_results.pkl")
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
    output_dir = '/Users/tehan/Documents/SURFSUP/NEID/NEIDSM'
    process_neid_fits(
            input_dir,
            output_dir,
            orders=('55', '101', '102', '103'),
            maxvsini=150,
            calibrate_feh=True,
            scaleres=1,
            deblazed=False,
            mode='HR',
            save_plot_data=True
    )
    gather_neidl2_pickle_results(output_dir)