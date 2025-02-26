import neidspec
import neidspecmatch
import neidspecmatch.config
import argparse
import pandas as pd

if __name__ == '__main__':
    # Make sure library is available, if not, download it
    neidspecmatch.get_library()
    # for i in range(55,104):
    for i in range(102,103):
        order = str(i)
        df_lib = pd.read_csv(neidspecmatch.config.PATH_LIBRARY_DB)
        HLS = neidspec.NEIDSpecList(filelist=neidspecmatch.config.LIBRARY_FITSFILES)
        outputdir = neidspecmatch.config.PATH_LIBRARY_CROSSVAL
        plot_results = True
        calibrate_feh = False
        scaleres = 2.

        # Run cross validation for orders
        neidspecmatch.run_crossvalidation_for_orders(order=order,
                                                     df_lib=df_lib,
                                                     HLS=HLS,
                                                     outputdir=outputdir,
                                                     plot_results=plot_results,
                                                     calibrate_feh=calibrate_feh,
                                                     scaleres=scaleres)

    # parser = argparse.ArgumentParser(description="NEID SpecMatch: Cross Validation")
    # parser.add_argument("order", type=int, default=101, help="Order to use for cross validation, e.g., order 101")
    # parser.add_argument("--df_lib", type=str, default=neidspecmatch.config.PATH_LIBRARY_DB,
    #                     help="Path to stellar library csv, e.g., 20201008_specmatch_nir_library.csv")
    # parser.add_argument("--HLS", type=str, default=neidspecmatch.config.LIBRARY_FITSFILES,
    #                     help="List of stellar library fits files, e.g., neidspecmatch.config.LIBRARY_FITSFILES")
    # parser.add_argument("--savefolder", type=str, default=neidspecmatch.config.PATH_LIBRARY_CROSSVAL,
    #                     help="Specify foldername to save (e.g. o17_crossval)")
    # parser.add_argument("--plot_results", default=True, help="Save cross validation summary plots", action="store_true")
    # parser.add_argument("--calibrate_feh", default=True, help="Calibrate the Fe/H", action="store_true")
    # parser.add_argument("--scaleres", type=float, default=1., help="Residual Scaling Factor")
    # argv = ["", neidspecmatch.config.PATH_LIBRARY_DB, neidspecmatch.config.LIBRARY_FITSFILES,
    #         neidspecmatch.config.PATH_LIBRARY_CROSSVAL, True, True, 1.]
    #
    # args = parser.parse_args(argv[1:])
