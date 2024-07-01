import neidspec
import neidspecmatch
import argparse

if __name__=='__main__':
    # parser = argparse.ArgumentParser(description="HPF SpecMatch: Match")
    # parser.add_argument("filename",type=str,help="filename to reduce")
    # parser.add_argument("object",type=str,help="HPF object to reduce (SIMBAD or TIC Queryable)")
    # parser.add_argument("--savefolder",type=str,default="specmatch_results",help="Specify foldername to save (e.g. results_123)")
    # parser.add_argument("--orders",type=int,default=[4,5,6,14,15,16,17],help="Orders to use for HPF SpecMatch, e.g., --orders 4 5 6",nargs='+')
    # parser.add_argument("--vsinimax",type=int,default=100.,help="Maximum vsini to fit for in km/s") #TODO: vsini limit
    # parser.add_argument("--calibrate_feh",default=False,help="Calibrate the Fe/H",action="store_true")
    # parser.add_argument("--scaleres",type=float,default=1.,help="Residual Scaling Factor")
    #
    # args = parser.parse_args()
    #
    # Make sure library is availabe, if not, download it
    neidspecmatch.get_library()

    filename = '/home/tehan/Downloads/neidL2_20230706T110904.fits'
    targetname = 'TOI 2145'
    outputdir = '/home/tehan/Downloads/specmatch_results'
    orders = ['55', '101', '102']
    path_df_lib = neidspecmatch.config.PATH_LIBRARY_DB
    maxvsini = 100
    calibrate_feh = False
    scaleres = 1

    # Run specmatch for orders
    neidspecmatch.run_specmatch_for_orders(targetfile=filename,
                                           targetname=targetname,
                                           outputdirectory=outputdir,
                                           path_df_lib=path_df_lib,
                                           orders=orders,
                                           maxvsini=maxvsini,
                                           calibrate_feh=calibrate_feh,
                                           scaleres=scaleres)
