import neidspecmatch
import numpy as np
import os
from astropy.io import fits
from glob import glob

if __name__ == '__main__':
    # Set input directory containing all FITS files
    input_dir = '/Users/tehan/Documents/SURFSUP/NEID/NEID_Spectra'  # or wherever your .fits files are
    output_dir = '/Users/tehan/Documents/SURFSUP/NEID/NEIDSM'
    orders = ['55', '101', '102', '103']
    path_df_lib = neidspecmatch.config.PATH_LIBRARY_DB
    path_df_lib_fits = neidspecmatch.config.PATH_LIBRARY_FITS
    maxvsini = 150
    calibrate_feh = True
    scaleres = 1
    deblazed = False
    mode = 'HR'
    save_plot_data = True

    # Make sure the spectral library is available
    neidspecmatch.get_library()

    # Process all .fits files in the input directory
    for filepath in glob(os.path.join(input_dir, '*.fits')):
        try:
            with fits.open(filepath) as hdul:
                targetname = hdul[0].header['OBJECT']
        except Exception as e:
            print(f"Skipping {filepath} due to error: {e}")
            continue

        outputdir = os.path.join(output_dir, targetname)
        os.makedirs(outputdir, exist_ok=True)

        print(f"Processing {filepath} for target {targetname}")
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