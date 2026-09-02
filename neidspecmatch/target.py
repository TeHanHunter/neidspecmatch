import configparser
import os
from pathlib import Path
import tempfile

from platformdirs import user_data_path

from . import bary
from .paths import safe_filename_component

DIRNAME = os.path.dirname(__file__)
PATH_TARGETS = os.fspath(user_data_path("neidspecmatch") / "targets")

class Target(object):
    """
    Simple target class. Capable of querying SIMBAD. Can calculate barycentric corrections.
    
    EXAMPLE:
        H = HPFSpectrum(fitsfiles[1])
        H.plot_order(14,deblazed=True)
        T = Target('G 9-40')
        T.calc_barycentric_velocity(H.jd_midpoint, 'Kitt Peak National Observatory')
        T = Target('G 9-40')
    """
    
    def __init__(self, name, config_folder=PATH_TARGETS, verbose=False,
                 obsname=None, allow_network=False):
        self.verbose = verbose
        requested_name = str(name)
        self.config_folder = os.fspath(Path(config_folder).expanduser())
        component = safe_filename_component(requested_name)
        self.config_filename = os.fspath(Path(self.config_folder) / f"{component}.config")
        if name=='Teegarden':
            name = "Teegarden's star"
        if name=='HR8926-4':
            name = 'GL87'
        if name=='GJ_1151':
            name = 'GJ 1151'
        if name=='GJ_324_A':
            name = 'GJ_324A'
        if name=='HD_68988':
            name = 'HD 68988'
        if name=='NLTT_51984':
            name = 'GJ_9751'
        self.name = name
        try:
            self.data = self.from_file(verbose=verbose)
        except FileNotFoundError:
            if not allow_network:
                raise
            if verbose:
                print('Target cache does not exist; querying a catalog.')
            if 'TIC' in name:
                if verbose:
                    print('Querying TIC for data')
                self.data = self.query_tic(name)
            else:
                if verbose:
                    print('Querying SIMBAD for data')
                try:
                    import barycorrpy
                except ImportError as exc:
                    raise ImportError(
                        "Catalog resolution requires optional dependencies: "
                        "pip install neidspecmatch[archive]"
                    ) from exc
                self.data, self.warning = barycorrpy.utils.get_stellar_data(name)
            self.to_file(self.data)
        self.ra = self.data['ra']
        self.dec = self.data['dec']
        self.pmra = self.data['pmra']
        self.pmdec = self.data['pmdec']
        self.px = self.data['px']
        self.epoch = self.data['epoch']
        if self.data['rv'] is None:
            self.rv = 0.
        else:
            self.rv = self.data['rv']/1000.# if self.data['rv'] < 1e20 else 0.
        self.obsname = obsname


    def query_tic(self,ticname):
        """
        Query the TESS Input Catalog for data
        """
        try:
            from astroquery.mast import Catalogs
        except ImportError as exc:
            raise ImportError(
                "TIC queries require optional dependencies: "
                "pip install neidspecmatch[archive]"
            ) from exc
        name = ticname.replace('-',' ').replace('_',' ')
        df = Catalogs.query_object(name, radius=0.0003, catalog="TIC").to_pandas()[0:1]
        data = {}
        data['ra'] = df.ra.values[0]
        data['dec'] = df.dec.values[0]
        data['pmra'] = df.pmRA.values[0]
        data['pmdec'] = df.pmDEC.values[0]
        data['px'] = df.plx.values[0]
        data['epoch'] = 2451545.0
        data['rv'] = 0.
        return data

    def from_file(self,verbose=False):
        if verbose:
            print('Reading from file {}'.format(self.config_filename))
        filename = Path(self.config_filename)
        if not filename.is_file():
            raise FileNotFoundError(filename)
        config = configparser.ConfigParser()
        with filename.open(encoding='utf-8') as stream:
            config.read_file(stream)
        data = dict(config.items('targetinfo'))
        required = {'ra', 'dec', 'pmra', 'pmdec', 'px', 'epoch', 'rv'}
        missing = sorted(required.difference(data))
        if missing:
            raise ValueError(
                "Malformed target cache; missing keys: {}".format(', '.join(missing))
            )
        for key, value in data.items():
            data[key] = None if value.strip().lower() == 'none' else float(value)
        return data

    def to_file(self,data):
        folder = Path(self.config_folder)
        folder.mkdir(parents=True, exist_ok=True, mode=0o700)
        try:
            folder.chmod(0o700)
        except OSError:
            pass
        if self.verbose:
            print('Saving to file {}'.format(self.config_filename))
        config = configparser.ConfigParser()
        config.add_section('targetinfo')
        for key in data.keys():
            config.set('targetinfo',key,str(data[key]))
            if self.verbose:
                print(key, data[key])
        descriptor, temporary_name = tempfile.mkstemp(
            prefix='.target-', suffix='.tmp', dir=folder
        )
        try:
            os.chmod(temporary_name, 0o600)
            with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
                config.write(stream)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary_name, self.config_filename)
            os.chmod(self.config_filename, 0o600)
        except BaseException:
            try:
                os.close(descriptor)
            except OSError:
                pass
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
            raise
        if self.verbose:
            print('Done')
        
    def calc_barycentric_velocity(self, jdtime, obs=None):
        """
        OUTPUT:
            BJD_TDB
            berv in km/s
        
        EXAMPLE:
            bjd, berv = bary.bjdbrv(H.jd_midpoint,T.ra,T.dec,obsname='Kitt Peak National Observatory',
                           pmra=T.pmra,pmdec=T.pmdec,rv=T.rv,parallax=T.px,epoch=T.epoch)
        """
        #bjd, berv = bary.bjdbrv(jdtime,self.ra,self.dec,obsname=self.obsname,
        #                           pmra=self.pmra,pmdec=self.pmdec,rv=self.rv,parallax=self.px,epoch=self.epoch)
        obs = obs or self.obsname
        if obs is None:
            raise ValueError(
                "An observatory name is required. For NEID L2 products, "
                "prefer the DRP SSBJD/SSBRV header values."
            )
        bjd, berv = bary.bjdbrv(jdtime, self.ra, self.dec, obsname=obs,
                                pmra=self.pmra, pmdec=self.pmdec, rv=self.rv, parallax=self.px, epoch=self.epoch)
        return bjd, berv/1000.
    
    def __repr__(self):
        return "{}, ra={:0.4f}, dec={:0.4f}, pmra={}, pmdec={}, rv={:0.4f}, px={:0.4f}, epoch={}".format(self.name,
                                            self.ra,self.dec,self.pmra,self.pmdec,self.rv,self.px,self.epoch)
