import os
import astropy.time
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pickle
import re
import pandas as pd
import hashlib
import json
from importlib import resources
from pathlib import Path, PurePosixPath
import shutil
import stat
import tempfile
from urllib.request import urlopen
import warnings
import zipfile
import neidspecmatch.config as config
from .paths import sha256_file
norm_mean = lambda x: x/np.nanmean(x)

def pickle_dump(filename, obj, *, verbose=False):
    """
    write obj to to filename and save
    
    INPUT:
        filename - name to save pickle file (str)
        obj - object to write to file
    
    OUTPUT:
        obj saved as pickle file
        
    EXAMPLE:
        pickle_dump()
    
    """
    warnings.warn(
        "Pickle is a legacy Python-specific format. Prefer JSON for durable "
        "NEIDSpecMatch results.",
        DeprecationWarning,
        stacklevel=2,
    )
    with open(filename, "wb") as savefile:
        pickle.dump(obj, savefile, protocol=pickle.HIGHEST_PROTOCOL)
    if verbose:
        print("Saved to {}".format(filename))

def pickle_load(filename, python3=True, *, trusted=False):
    """
    load obj from filename
    
    INPUT:
        filename - name of saved pickle file (str)
        obj - object to write to file
    
    OUTPUT:
        load obj from pickle file
        
    EXAMPLE:
        pickle_load()
    """
    if not trusted:
        raise ValueError(
            "Refusing to load pickle without trusted=True. Pickle can execute "
            "arbitrary code; load only files from a trusted source."
        )
    warnings.warn(
        "Loading legacy pickle; caller asserted that the file is trusted.",
        RuntimeWarning,
        stacklevel=2,
    )
    with open(filename, "rb") as openfile:
        return pickle.load(openfile, encoding='latin1')


def jd2datetime(times):
    """
    convert jd to iso utc date/time
    
    INPUT:
        times in jd (array)
    
    OUTPUT:
        date/time in utc (array)
        
    EXAMPLE:
        
    """
    return np.array([astropy.time.Time(time,format="jd",scale="utc").datetime for time in times])

def iso2jd(times):
    """
    convert iso utc to jd date/time
    
    INPUT:
        date/time in iso uts (array)
    
    OUTPUT:
        date/time in jd (array)
        
    EXAMPLE:
        
    """
    return np.array([astropy.time.Time(time,format="iso",scale="utc").jd for time in times])

def make_dir(dirname,verbose=True):
    """
    create new directory
    
    INPUT:
        dirname - name of sirectory (str)
        verbose - include print statements (bool)
    
    OUTPUT:
        creates folder named dirname
        
    EXAMPLE:
        make_dir("folder")
        
    """
    try:
        os.makedirs(dirname)
        if verbose==True: print("Created folder:",dirname)
    except OSError:
        if verbose==True: print(dirname,"already exists.")

def vac2air(wavelength,P,T,input_in_angstroms=True):
    """
    Convert vacuum wavelengths to air wavelengths

    INPUT:
        wavelength - in A if input_in_angstroms is True, else nm
        P - in Torr
        T - in Celsius

    OUTPUT:
        Wavelength in air in A if input_in_angstroms is True, else nm
    """
    if input_in_angstroms:
        nn = n_air(P,T,wavelength/10.)
    else:
        nn = n_air(P,T,wavelength)
    return wavelength/(nn+1.)

def n_air(P,T,wavelength):
    """
    The edlen equation for index of refraction of air with pressure
    
    INPUT:
        P - pressure in Torr
        T - Temperature in Celsius
        wavelength - wavelength in nm
        
    OUTPUT:
        (n-1)_tp - see equation 1, and 2 in REF below.
        
    REF:
        http://iopscience.iop.org/article/10.1088/0026-1394/30/3/004/pdf

    EXAMPLE:
        nn = n_air(763.,20.,500.)-n_air(760.,20.,500.)
        (nn/(nn + 1.))*3.e8
    """
    wavenum = 1000./wavelength # in 1/micron
    # refractivity is (n-1)_s for a standard air mixture
    refractivity = ( 8342.13 + 2406030./(130.-wavenum**2.) + 15997./(38.9-wavenum**2.))*1e-8
    return ((P*refractivity)/720.775) * ( (1.+P*(0.817-0.0133*T)*1e-6) / (1. + 0.0036610*T) )

def get_cmap_colors(cmap='jet',p=None,N=10):
    """
    get colormap instance
    
    INPUT:
        cmap - colormap instance
        p -
        N -
    
    OUTPUT:
        
        
    EXAMPLE:
        get_cmap_colors()
    """
    cm = plt.get_cmap(cmap)
    if p is None:
        return [cm(i) for i in np.linspace(0,1,N)]
    else:
        normalize = matplotlib.colors.Normalize(vmin=min(p), vmax=max(p))
        colors = [cm(normalize(value)) for value in p]
        return colors

def ax_apply_settings(ax,ticksize=None):
    """
    Apply axis settings that I keep applying
    """
    ax.minorticks_on()
    if ticksize is None:
        ticksize=12
    ax.tick_params(pad=3,labelsize=ticksize)
    ax.grid(lw=0.3,alpha=0.3)

def ax_add_colorbar(ax,p,cmap='jet',tick_width=1,tick_length=3,direction='out',pad=0.02,minorticks=False,*kw_args):
    """
    Add a colorbar to a plot (e.g., of lines)

    INPUT:
        ax - axis to put the colorbar
        p  - parameter that will be used for the scaling (min and max)
        cmap - 'jet', 'viridis'

    OUTPUT:
        cax - the colorbar object

    NOTES:
        also see here:
            import matplotlib.pyplot as plt
            sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=0, vmax=1))
            # fake up the array of the scalar mappable. Urgh...
            sm._A = []
            plt.colorbar(sm)

    EXAMPLE:
        cmap = 'jet'
        colors = get_cmap_colors(N=len(bjds),cmap=cmap)
        fig, ax = plt.subplots(dpi=200)
        for i in range(len(bjds)):
            ax.plot(vin[i],ccf_sub[i],color=colors[i],lw=1)
        ax_apply_settings(ax,ticksize=14)
        ax.set_xlabel('Velocity [km/s]',fontsize=20)
        ax.set_ylabel('Relative Flux ',fontsize=20)
        cx = ax_add_colorbar(ax,obs_phases,cmap=cmap,pad=0.02)
        ax_apply_settings(cx)
        cx.axes.set_ylabel('Phase',fontsize=20,labelpad=0)


        cbar.set_clim(1.4,4.1)
        cbar.set_ticks([1.5,2.0,2.5,3.0,3.5,4.0])
    """
    cax, _ = matplotlib.colorbar.make_axes(ax,pad=pad,*kw_args)
    normalize = matplotlib.colors.Normalize(vmin=np.nanmin(p), vmax=np.nanmax(p))
    cm = plt.get_cmap(cmap)
    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cm, norm=normalize)
    if minorticks:
        cax.minorticks_on()
    cax.axes.tick_params(width=tick_width,length=tick_length,direction=direction)
    return cax, cbar

def get_indices_of_items(arr,items):
    return np.where(pd.DataFrame(arr).isin(items))[0]


def remove_items_from_list(l,bad_items):
    ibad = np.where(pd.DataFrame(l).isin(bad_items))[0]
    return np.delete(l,ibad)

def savefigure(fig,savename,s1='{}',p1='',s2='{}',p2='',dpi=200):
    """
    Handy function to save figures and append suffixes to filenames
    
    EXAMPLE:
        savefigure(fig,'MASTER_FLATS/COMPARE_PLOTS/testing.png',s1='_o{}',p1=5,s2='_spi{}',p2=14)
    """
    fp = FilePath(savename)
    make_dir(fp.directory)
    fp.add_suffix(s1.format(p1))
    fp.add_suffix(s2.format(p2))
    fig.tight_layout()
    fig.savefig(fp._fullpath,dpi=dpi)
    print('Saved figure to: {}'.format(fp._fullpath))

def grep_date(string,intype="isot",outtype='iso'):
    """
    A function to extract date from string.

    INPUT:
        string: string
        intype: "isot" - 20181012T001823
                "iso"  - 20181012
        outtype: "iso" - iso
                 "datetime" - datetime

    OUTPUT:
        string with the date
    """
    if intype == "isot":
        date = re.findall(r'\d{8}T\d{6}',string)[0]
    elif intype == "iso":
        date = re.findall(r'\d{8}',string)[0]
    else:
        print("intype has to be 'isot' or 'iso'")
    if outtype == 'iso':
        return date
    elif outtype == 'datetime':
        return pd.to_datetime(date).to_pydatetime()
    else:
        print('outtype has to be "iso" or "datetime"')

def grep_dates(strings,intype="isot",outtype='iso'):
    """
    A function to extract date from strings

    INPUT:
        string: string
        intype: "isot" - 20181012T001823
                "iso"  - 20181012
        outtype: "iso" - iso
                 "datetime" - datetime

    OUTPUT:
        string with the date

    EXAMPLE:
        df = grep_dates(files,intype="isot",outtype='series')
        df['2018-06-26':].values
    """
    if outtype=='series':
        dates = [grep_date(i,intype=intype,outtype='datetime') for i in strings]
        return pd.Series(index=dates,data=strings)
    else:
        return [grep_date(i,intype,outtype) for i in strings]

def replace_dir(files,old,new):
    for i,f in enumerate(files):
        files[i] = files[i].replace(old,new)
    return files

def get_header_df(fitsfiles,keywords=["OBJECT","DATE-OBS"],verbose=True):
    """
    A function to read headers and returns a pandas dataframe

    INPUT:
    fitsfiles - a list of fitsfiles
    keywords - a list of header keywords to read

    OUTPUT:
    df - pandas dataframe with column names as the passed keywords
    """
    headers = []
    for i,name in enumerate(fitsfiles):
        if verbose: print(i,name)
        head = astropy.io.fits.getheader(name)
        values = [name]
        for key in keywords:
            values.append(head[key])
        headers.append(values)
    df_header = pd.DataFrame(headers,columns=["filename"]+keywords)
    return df_header

def plot_rv_model(bjd,RV,e_RV,P,T0,K,e,w,nbjd=None,nRV=None,ne_RV=None,title="",fig=None,ax=None,bx=None):
    """
    Plot RVs and the rv model on top. Also plot residuals
    """
    bjd_model = np.linspace(bjd[0]-2,bjd[-1]+2,10000)
    rv_model = get_rv_curve(bjd_model,P,T0,e=e,omega=w,K=K,plot=False)
    rv_obs = get_rv_curve(bjd,P,T0,e,w,K,plot=False)
    
    if nbjd is not None:
        rv_obs_bin = get_rv_curve(nbjd,P,T0,e,w,K,plot=False)

    # Plotting
    if fig is None and ax is None and bx is None:
        fig, (ax, bx) = plt.subplots(dpi=200,nrows=2,gridspec_kw={"height_ratios":[5,2]},figsize=(6,4),sharex=True)
    
    ax.errorbar(bjd,RV,e_RV,elinewidth=1,mew=0.5,capsize=5,marker="o",lw=0,markersize=6,alpha=0.5,
            label="Unbin, median errorbar={:0.1f}m/s".format(np.median(e_RV)))
    ax.plot(bjd_model,rv_model,label="Expected Orbit",lw=1,color="crimson")
    res = RV - rv_obs
    bx.errorbar(bjd,res,e_RV,elinewidth=1,mew=0.5,capsize=5,marker="o",lw=0,markersize=6,alpha=0.5,label=r'Residual: $\sigma$={:0.2f}m/s'.format(np.std(res)))
    
    if nbjd is not None:
        ax.errorbar(nbjd,nRV,ne_RV,elinewidth=1,mew=0.5,capsize=5,marker="h",lw=0,color="crimson",
            markersize=8,label="Bin, median errorbar={:0.1f}m/s".format(np.median(ne_RV)))
        bx.errorbar(nbjd,nRV-rv_obs_bin,ne_RV,elinewidth=1,mew=0.5,capsize=5,marker="h",lw=0,markersize=8,alpha=0.5,color="crimson")
    ax.legend(fontsize=8,loc="upper right")
    bx.legend(fontsize=8,loc="upper right")

    for xx in [ax,bx]:
        ax_apply_settings(xx,ticksize=12)
        xx.set_ylabel("RV [m/s]",fontsize=16)

    bx.set_xlabel("Date [UT]",labelpad=0,fontsize=16)
    ax.set_title(title)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.05)
    
def plot_rv_model_phased(bjd,RV,e_RV,P,T0,K,e,w,nbjd=None,nRV=None,ne_RV=None,title="",fig=None,ax=None,bx=None):
    """
    Plot RVs and the rv model on top. Also plot residuals
    """
    bjd_model = np.linspace(bjd[0]-2,bjd[-1]+2,10000)
    rv_model = get_rv_curve(bjd_model,P,T0,e=e,omega=w,K=K,plot=False)
    rv_obs = get_rv_curve(bjd,P,T0,e,w,K,plot=False)
    
    # Phases:
    df_phase = get_phases_sorted(bjd,P,T0).sort_values("time")
    df_phase_model = get_phases_sorted(bjd_model,P,T0,rvs=rv_model).sort_values("time")

    if nbjd is not None:
        rv_obs_bin = get_rv_curve(nbjd,P,T0,e,w,K,plot=False)
        df_phase_bin = get_phases_sorted(nbjd,P,T0).sort_values("time")

    # Plotting
    if fig is None and ax is None and bx is None:
        fig, (ax, bx) = plt.subplots(dpi=200,nrows=2,gridspec_kw={"height_ratios":[5,2]},figsize=(6,4),sharex=True)
    
    ax.errorbar(df_phase.phases,RV,e_RV,elinewidth=1,mew=0.5,capsize=5,marker="o",lw=0,markersize=6,alpha=0.5,
            label="Unbin, median errorbar={:0.1f}m/s".format(np.median(e_RV)))
    ax.plot(df_phase_model.phases,rv_model,label="Expected Orbit",lw=0,marker="o",markersize=2,color="crimson")
    res = RV - rv_obs
    bx.errorbar(df_phase.phases,res,e_RV,elinewidth=1,mew=0.5,capsize=5,marker="o",lw=0,markersize=6,alpha=0.5,label=r'Residual: $\sigma$={:0.2f}m/s'.format(np.std(res)))
    
    if nbjd is not None:
        ax.errorbar(df_phase_bin.phases,nRV,ne_RV,elinewidth=1,mew=0.5,capsize=5,marker="h",lw=0,color="crimson",
            markersize=8,label="Bin, median errorbar={:0.1f}m/s".format(np.median(ne_RV)),alpha=0.9)
        bx.errorbar(df_phase_bin.phases,nRV-rv_obs_bin,ne_RV,elinewidth=1,mew=0.5,capsize=5,marker="h",lw=0,markersize=8,alpha=0.9,color="crimson")
    ax.legend(fontsize=8,loc="upper right")
    bx.legend(fontsize=8,loc="upper right")

    for xx in [ax,bx]:
        ax_apply_settings(xx,ticksize=12)
        xx.set_ylabel("RV [m/s]",fontsize=16)

    bx.set_xlabel("Orbital phase",labelpad=0,fontsize=16)
    ax.set_title(title)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.05)

def get_phases_sorted(t, P, t0,rvs=None,rvs_err=None,sort=True,centered_on_0=True,tdur=None):
    """
    Get a sorted pandas dataframe of phases, times (and Rvs if supplied)
    
    INPUT:
    t  - times in jd
    P  - period in days
    t0 - time of periastron usually
    
    OUTPUT:
    df - pandas dataframe with columns:
     -- phases (sorted)
     -- time - time
     -- rvs  - if provided
    
    NOTES:
    Useful for RVs.    
    """
    phases = np.mod(t - t0,P)
    phases /= P
    df = pd.DataFrame(zip(phases,t),columns=['phases','time'])
    if rvs is not None:
        df['rvs'] = rvs
    if rvs_err is not None:
        df['rvs_err'] = rvs_err
    if centered_on_0:
        _p = df.phases.values
        m = df.phases.values > 0.5
        _p[m] = _p[m] - 1.
        df['phases'] = _p
    if tdur is not None:
        df["intransit"] = np.abs(df.phases) < tdur/(2.*P)
        print("Found {} in transit".format(len(df[df['intransit']])))
    if sort:
        df = df.sort_values('phases').reset_index(drop=True)
    return df

def get_rv_curve(times_jd,P,tc,e,omega,K,plot=True,ax=None,verbose=True,plot_tnow=True):
    """
    A function to plot an RV curve as a function of time (not phased)
    
    INPUT:
        times_jd: times in jd
        P: orbital period in days
        tc: transit center in jd
        e: eccentricity
        omega: periastron in degrees
        K: RV semi-amplitude in m/s
    
    OUTPUT:
        rv: array of RVs at times times_jd
    """
    try:
        import radvel
    except ImportError as exc:
        raise ImportError(
            "get_rv_curve requires optional dependencies: "
            "pip install neidspecmatch[utilities]"
        ) from exc
    t_peri = radvel.orbit.timetrans_to_timeperi(tc=tc,
                                                per=P,
                                                ecc=e,
                                                omega=np.deg2rad(omega))
    rvs = radvel.kepler.rv_drive(times_jd,[P,
                                            t_peri,
                                            e,
                                            np.deg2rad(omega),
                                            K])
    if verbose:
        print("Assuming:")
        print("P {}d".format(P))
        print("tc {}".format(tc))
        print("e {}".format(e))
        print("omega {}deg".format(omega))
        print("K {}m/s".format(K))
    if plot:
        if ax is None:
            fig, ax = plt.subplots(figsize=(12,8))
        times = jd2datetime(times_jd)
        ax.plot(times,rvs)
        ax.set_xlabel("Time")
        ax.set_ylabel("RV [m/s]")
        ax.grid(lw=0.5,alpha=0.3)
        ax.minorticks_on()
        if plot_tnow:
            t_now = Time(datetime.datetime.utcnow())
            ax.axvline(t_now.datetime,color="red",label="Time now")
        xlim = ax.get_xlim()
        ax.hlines(0,xlim[0],xlim[1],alpha=0.5,color="k",lw=1)
        for label in (ax.get_xticklabels()):
            label.set_fontsize(10)
    return rvs


def filter_simbadnames(l):
    """
    Filter SIMBAD names
    """
    _n = find_str_in_list(l,'DR2')
    if _n=='':
        _n = find_str_in_list(l,'DR1')
    if _n=='':
        _n = find_str_in_list(l,'HD')
    if _n=='':
        _n = find_str_in_list(l,'GJ')
    return _n

def get_simbad_object_names(name,as_str=False):
    """
    Get all Simbad names for an object
    """
    try:
        from astroquery.simbad import Simbad
    except ImportError as exc:
        raise ImportError(
            "SIMBAD queries require optional dependencies: "
            "pip install neidspecmatch[archive]"
        ) from exc
    result_table = Simbad.query_objectids(name).to_pandas()
    names = result_table.ID.values.astype(str)
    if as_str:
        names = [i.replace(' ','_') for i in names]
        return '|'.join(names)
    else:
        return names

def find_str_in_list(l,string):
    """
    Return first element that matches a string in a list
    """
    for element in l:
        if string in element:
            return element
    return ''

LIBRARY_MANIFEST_SCHEMA = 2
DEFAULT_MAX_ZIP_MEMBERS = 512
DEFAULT_MAX_EXPANDED_BYTES = 64 * 1024 ** 3
DEFAULT_MAX_DOWNLOAD_BYTES = 16 * 1024 ** 3
DEFAULT_LIBRARY_ALLOWLIST_RESOURCE = "default_library_allowlist.json"


def _library_data_files(path, fits_files, catalog_name):
    return [path / catalog_name, *fits_files]


def _manifest_file_records(path, fits_files, catalog_name):
    records = []
    for filename in _library_data_files(path, fits_files, catalog_name):
        records.append({
            "path": filename.relative_to(path).as_posix(),
            "size_bytes": filename.stat().st_size,
            "sha256": sha256_file(filename),
        })
    return records


def load_default_library_allowlist():
    """Load and validate the release-owned default-library SHA-256 allowlist.

    Unlike the Zenodo-published MD5, these per-file digests are shipped with
    the NEIDSpecMatch release and authenticate the exact catalog/FITS payload
    that this code was reviewed against.
    """
    resource = resources.files("neidspecmatch.data").joinpath(
        DEFAULT_LIBRARY_ALLOWLIST_RESOURCE
    )
    payload = resource.read_bytes()
    allowlist = json.loads(payload.decode("utf-8"))
    if not isinstance(allowlist, dict):
        raise ValueError("Default library allowlist must be a JSON object.")
    if int(allowlist.get("schema_version", -1)) != LIBRARY_MANIFEST_SCHEMA:
        raise ValueError("Default library allowlist schema is unsupported.")
    expected_fields = {
        "library_id": config.DEFAULT_LIBRARY_ID,
        "catalog": config.DEFAULT_LIBRARY_CATALOG,
        "fits_count": 78,
        "source_url": config.URL_LIBRARY,
        "source_archive_size_bytes": config.LIBRARY_ZIP_SIZE_BYTES,
    }
    for field, expected in expected_fields.items():
        if allowlist.get(field) != expected:
            raise ValueError(
                f"Default library allowlist has unexpected {field!r}."
            )
    integrity = allowlist.get("source_archive_integrity")
    if not isinstance(integrity, dict) or (
            integrity.get("algorithm") != "md5"
            or integrity.get("value") != config.LIBRARY_ZIP_MD5
            or integrity.get("scope") != "integrity_only_not_authentication"):
        raise ValueError("Default library archive-integrity metadata is invalid.")
    if not str(allowlist.get("allowlist_id", "")).strip():
        raise ValueError("Default library allowlist has no identifier.")
    return allowlist, hashlib.sha256(payload).hexdigest()


def _verify_library_against_allowlist(
        library_path, allowlist, *, allow_installed_manifest=False):
    """Verify an extracted library against a release-owned file allowlist."""
    path = Path(library_path)
    _require_regular_path(path, "library directory", directory=True)
    if not isinstance(allowlist, dict):
        raise ValueError("Library allowlist must be a mapping.")
    records = allowlist.get("files")
    if not isinstance(records, list) or not records:
        raise ValueError("Library allowlist is missing per-file records.")
    expected = {}
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("Library allowlist contains an invalid record.")
        relative_text = str(record.get("path", ""))
        relative = PurePosixPath(relative_text)
        if (
                not relative_text or relative.is_absolute()
                or ".." in relative.parts or "." in relative.parts):
            raise ValueError(f"Unsafe library allowlist path: {relative_text!r}")
        if relative_text in expected:
            raise ValueError("Library allowlist contains duplicate file paths.")
        expected[relative_text] = record
    catalog = str(allowlist.get("catalog", ""))
    fits_count = int(allowlist.get("fits_count", -1))
    expected_paths = {catalog}
    expected_paths.update(
        name for name in expected if PurePosixPath(name).parts[:1] == ("FITS",)
    )
    if set(expected) != expected_paths or len(expected_paths) != fits_count + 1:
        raise ValueError(
            "Library allowlist must contain exactly one catalog and the declared FITS files."
        )
    actual = set()
    for item in path.rglob("*"):
        if item.is_symlink():
            raise ValueError(f"Library payload cannot contain a symbolic link: {item}")
        if item.is_file():
            _require_regular_path(item, "library payload file")
            relative = item.relative_to(path).as_posix()
            if allow_installed_manifest and relative == "library_manifest.json":
                continue
            actual.add(relative)
    if actual != set(expected):
        missing = sorted(set(expected).difference(actual))
        extra = sorted(actual.difference(expected))
        preview_limit = 20
        raise ValueError(
            "Library payload does not match the release allowlist "
            f"(missing_count={len(missing)}, extra_count={len(extra)}, "
            f"missing_preview={missing[:preview_limit]}, "
            f"extra_preview={extra[:preview_limit]})."
        )
    for relative, record in expected.items():
        filename = path / relative
        if filename.stat().st_size != int(record.get("size_bytes", -1)):
            raise ValueError(f"Release-allowlist size mismatch: {relative}")
        expected_digest = str(record.get("sha256", ""))
        if not re.fullmatch(r"[0-9a-f]{64}", expected_digest):
            raise ValueError(f"Invalid release-allowlist SHA-256: {relative}")
        if sha256_file(filename) != expected_digest:
            raise ValueError(f"Release-allowlist SHA-256 mismatch: {relative}")
    return {
        "allowlist_id": str(allowlist.get("allowlist_id")),
        "allowlist_sha256_verified": True,
        "file_count": len(expected),
    }


def _require_regular_path(path, description, *, directory=False):
    path = Path(path)
    if path.is_symlink():
        raise ValueError(f"{description} cannot be a symbolic link: {path}")
    predicate = path.is_dir if directory else path.is_file
    if not predicate():
        raise FileNotFoundError(f"Missing {description}: {path}")
    mode = path.stat().st_mode
    expected = stat.S_ISDIR(mode) if directory else stat.S_ISREG(mode)
    if not expected:
        raise ValueError(f"{description} must be a regular {'directory' if directory else 'file'}: {path}")


def build_library_manifest(
        library_path, *, library_id, catalog_name=config.DEFAULT_LIBRARY_CATALOG,
        expected_fits=None, source_url=None, source_archive_integrity=None,
        release_allowlist=None, header_identity_exceptions=None):
    """Build a deterministic per-file manifest for any explicit library.

    Source/archive metadata are optional provenance.  They are never used as a
    substitute for the catalog and per-file SHA-256 records.
    """
    path = Path(config.resolve_library_path(library_path))
    _require_regular_path(path, "library directory", directory=True)
    if not str(library_id).strip():
        raise ValueError("library_id must be non-empty.")
    if Path(catalog_name).name != str(catalog_name):
        raise ValueError("catalog_name must be a basename, not a path.")
    catalog = path / str(catalog_name)
    fits_dir = path / "FITS"
    _require_regular_path(catalog, "library catalog")
    _require_regular_path(fits_dir, "library spectra directory", directory=True)
    fits_files = sorted(fits_dir.glob("*.fits"))
    for filename in fits_files:
        _require_regular_path(filename, "library FITS file")
    if expected_fits is not None and len(fits_files) != int(expected_fits):
        raise ValueError(
            f"Expected {expected_fits} library FITS files; found {len(fits_files)}."
        )
    manifest = {
        "schema_version": LIBRARY_MANIFEST_SCHEMA,
        "library_id": str(library_id),
        "catalog": str(catalog_name),
        "fits_count": len(fits_files),
        "files": _manifest_file_records(path, fits_files, str(catalog_name)),
    }
    if source_url is not None:
        manifest["source_url"] = str(source_url)
    if source_archive_integrity is not None:
        manifest["source_archive_integrity"] = dict(source_archive_integrity)
    if release_allowlist is not None:
        receipt = dict(release_allowlist)
        if (
                not str(receipt.get("id", "")).strip()
                or not re.fullmatch(r"[0-9a-f]{64}", str(receipt.get("sha256", "")))):
            raise ValueError("Release allowlist receipt is invalid.")
        manifest["release_allowlist"] = receipt
    if header_identity_exceptions is not None:
        manifest["header_identity_exceptions"] = list(
            header_identity_exceptions
        )
    manifest_path = path / "library_manifest.json"
    if manifest_path.is_symlink():
        raise ValueError("Library manifest cannot be a symbolic link.")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=".library_manifest-", suffix=".tmp", dir=path
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(manifest, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, manifest_path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise
    return manifest_path


def validate_library(
        library_path=None, expected_fits=None, require_manifest=False, *,
        deep=False, validate_fits_schema=False,
        expected_library_id=None, catalog_name=None):
    """Validate the catalog/FITS layout and its reproducibility manifest.

    Manifest file sizes are checked whenever a manifest exists.  ``deep=True``
    additionally recomputes every catalog/FITS SHA-256.  FITS structure and
    catalog-to-basename mapping are an explicit, slower validation step.
    """
    path = Path(config.resolve_library_path(library_path))
    if path.is_symlink():
        raise ValueError(f"Library path cannot be a symbolic link: {path}")
    _require_regular_path(path, "library directory", directory=True)
    manifest_path = path / "library_manifest.json"
    if manifest_path.is_symlink():
        raise ValueError("Library manifest cannot be a symbolic link.")
    manifest = None
    if manifest_path.is_file():
        _require_regular_path(manifest_path, "library manifest")
        with manifest_path.open(encoding="utf-8") as stream:
            manifest = json.load(stream)
    manifest_catalog = manifest.get("catalog") if isinstance(manifest, dict) else None
    selected_catalog = (
        str(catalog_name) if catalog_name is not None
        else str(manifest_catalog or config.DEFAULT_LIBRARY_CATALOG)
    )
    if Path(selected_catalog).name != selected_catalog:
        raise ValueError("Library catalog must be a basename.")
    catalog = path / selected_catalog
    fits_dir = path / "FITS"
    _require_regular_path(catalog, "library catalog")
    _require_regular_path(fits_dir, "library spectra directory", directory=True)
    fits_files = sorted(fits_dir.glob("*.fits"))
    for filename in fits_files:
        _require_regular_path(filename, "library FITS file")
    if expected_fits is not None and len(fits_files) != int(expected_fits):
        raise ValueError(
            f"Expected {expected_fits} library FITS files in {fits_dir}; "
            f"found {len(fits_files)}."
        )
    if require_manifest and not manifest_path.is_file():
        raise FileNotFoundError(f"Missing library manifest: {manifest_path}")
    if manifest is not None:
        if int(manifest.get("schema_version", -1)) != LIBRARY_MANIFEST_SCHEMA:
            raise ValueError(
                "Library manifest schema is unsupported; reinstall the library "
                "to generate per-file integrity records."
            )
        if not str(manifest.get("library_id", "")).strip():
            raise ValueError("Library manifest has no non-empty library_id.")
        if (
            expected_library_id is not None
            and str(manifest.get("library_id")) != str(expected_library_id)
        ):
            raise ValueError("Library manifest has an unexpected library_id.")
        if str(manifest.get("catalog")) != selected_catalog:
            raise ValueError("Library manifest catalog does not match the selected catalog.")
        if int(manifest.get("fits_count", -1)) != len(fits_files):
            raise ValueError("Library manifest FITS count does not match the files.")
        archive_integrity = manifest.get("source_archive_integrity")
        if archive_integrity is not None:
            if not isinstance(archive_integrity, dict):
                raise ValueError("Library source archive integrity metadata is invalid.")
            if not all(
                str(archive_integrity.get(key, "")).strip()
                for key in ("algorithm", "value", "scope")
            ):
                raise ValueError("Library source archive integrity metadata is incomplete.")
        release_allowlist = manifest.get("release_allowlist")
        if release_allowlist is not None:
            if not isinstance(release_allowlist, dict) or (
                    not str(release_allowlist.get("id", "")).strip()
                    or not re.fullmatch(
                        r"[0-9a-f]{64}", str(release_allowlist.get("sha256", ""))
                    )):
                raise ValueError("Library release allowlist receipt is invalid.")
        records = manifest.get("files")
        if not isinstance(records, list):
            raise ValueError("Library manifest is missing per-file records.")
        expected_paths = {
            filename.relative_to(path).as_posix(): filename
            for filename in _library_data_files(path, fits_files, selected_catalog)
        }
        record_paths = [record.get("path") for record in records]
        if len(record_paths) != len(set(record_paths)):
            raise ValueError("Library manifest contains duplicate file paths.")
        if set(record_paths) != set(expected_paths):
            raise ValueError("Library manifest file list does not match the library.")
        for record in records:
            filename = expected_paths[record["path"]]
            if filename.stat().st_size != int(record.get("size_bytes", -1)):
                raise ValueError(f"Library file size mismatch: {record['path']}")
            digest = record.get("sha256", "")
            if not re.fullmatch(r"[0-9a-f]{64}", str(digest)):
                raise ValueError(f"Invalid SHA-256 record: {record['path']}")
            if deep and sha256_file(filename) != digest:
                raise ValueError(f"Library SHA-256 mismatch: {record['path']}")
    schema_result = None
    if validate_fits_schema:
        schema_result = validate_library_fits_schema(
            path, catalog_name=selected_catalog,
            header_identity_exceptions=(
                manifest.get("header_identity_exceptions", [])
                if manifest is not None else []
            ),
        )
    return {
        "library_path": str(path),
        "catalog": str(catalog),
        "fits_path": str(fits_dir),
        "fits_count": len(fits_files),
        "library_id": manifest.get("library_id") if manifest else None,
        "catalog_name": selected_catalog,
        "manifest": manifest,
        "fits_schema": schema_result,
    }


def _catalog_basename(row):
    for column in ('basenames', 'basename', 'filename'):
        value = row.get(column)
        if isinstance(value, str) and value.strip():
            basename = os.path.basename(value.strip())
            if basename != value.strip().replace('\\', '/') .split('/')[-1]:
                raise ValueError(f"Unsafe catalog FITS path: {value!r}")
            return basename
    raise ValueError("Every library row needs a FITS basename column.")


def _normalized_object_name(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _catalog_object_aliases(row):
    aliases = set()
    for column in (
            "OBJECT", "OBJECT_ID", "ID_NAME", "simbadnames_x",
            "simbadnames_y"):
        value = row.get(column)
        if not isinstance(value, str):
            continue
        aliases.update(
            _normalized_object_name(item)
            for item in value.split("|") if item.strip()
        )
    return aliases


def _gaia_source_ids(value):
    return set(re.findall(
        r"(?i)gaia\s*(?:dr[23]\s*)?(\d{12,})", str(value)
    ))


def _validated_header_identity_exceptions(exceptions):
    if exceptions is None:
        return {}
    if not isinstance(exceptions, list):
        raise ValueError("header_identity_exceptions must be a list.")
    required = {
        "path", "catalog_object_id", "fits_object", "gaia_source_id",
        "reason",
    }
    parsed = {}
    for exception in exceptions:
        if not isinstance(exception, dict) or set(exception) != required:
            raise ValueError(
                "Each header identity exception must contain exactly "
                f"{sorted(required)}."
            )
        path_text = str(exception["path"])
        path = PurePosixPath(path_text)
        if (
                path.is_absolute() or ".." in path.parts
                or path.parts[:1] != ("FITS",) or len(path.parts) != 2
                or path.suffix.lower() != ".fits"):
            raise ValueError(
                f"Unsafe header identity exception path: {path_text!r}"
            )
        if path_text in parsed:
            raise ValueError("Duplicate header identity exception path.")
        gaia_source_id = str(exception["gaia_source_id"])
        if not re.fullmatch(r"\d{12,}", gaia_source_id):
            raise ValueError("Invalid Gaia source ID in identity exception.")
        for field in ("catalog_object_id", "fits_object", "reason"):
            if not str(exception[field]).strip():
                raise ValueError(
                    f"Empty {field} in header identity exception."
                )
        parsed[path_text] = dict(exception)
    return parsed


def validate_library_fits_schema(
        library_path=None, *, catalog_name=config.DEFAULT_LIBRARY_CATALOG,
        header_identity_exceptions=None):
    """Deeply validate NEID FITS structure and catalog/object mapping."""
    from astropy.io import fits

    path = Path(config.resolve_library_path(library_path))
    if Path(catalog_name).name != str(catalog_name):
        raise ValueError("catalog_name must be a basename.")
    catalog_path = path / str(catalog_name)
    table = pd.read_csv(catalog_path)
    if 'OBJECT_ID' not in table:
        raise ValueError("Library catalog is missing OBJECT_ID.")
    if table['OBJECT_ID'].isna().any() or table['OBJECT_ID'].astype(str).duplicated().any():
        raise ValueError("Library OBJECT_ID values must be nonempty and unique.")
    mapping = []
    identity_exceptions = _validated_header_identity_exceptions(
        header_identity_exceptions
    )
    used_identity_exceptions = set()
    seen_basenames = set()
    for _, row in table.iterrows():
        basename = _catalog_basename(row)
        if basename in seen_basenames:
            raise ValueError(f"Duplicate catalog FITS basename: {basename}")
        seen_basenames.add(basename)
        filename = path / 'FITS' / basename
        _require_regular_path(filename, "catalog-mapped library FITS file")
        with fits.open(filename, memmap=False) as hdus:
            if len(hdus) <= 17:
                raise ValueError(f"NEID L2 FITS has too few HDUs: {basename}")
            header = hdus[0].header
            required_headers = {
                'INSTRUME': 'NEID', 'OBSTYPE': 'SCI', 'OBS-MODE': 'HR',
            }
            for keyword, expected in required_headers.items():
                if str(header.get(keyword, '')).strip().upper() != expected:
                    raise ValueError(
                        f"NEID library FITS has invalid {keyword}: {basename}"
                    )
            try:
                data_level = int(header.get('DATALVL'))
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"NEID library FITS has invalid DATALVL: {basename}") from exc
            if data_level != 2:
                raise ValueError(f"NEID library FITS is not Level 2: {basename}")
            if not str(header.get('E_VER', '')).strip():
                raise ValueError(f"NEID library FITS has no DRP E_VER: {basename}")
            for keyword in ('DQLEVEL1', 'DQLEVEL2'):
                try:
                    dq_value = int(header.get(keyword))
                except (TypeError, ValueError, OverflowError) as exc:
                    raise ValueError(
                        f"NEID library FITS has invalid {keyword}: {basename}"
                    ) from exc
                if dq_value < 0 or (dq_value & 0b11) in {2, 3}:
                    raise ValueError(
                        f"NEID library FITS fails {keyword}: {basename}"
                    )
            fits_object = str(header.get('OBJECT', '')).strip()
            if not fits_object:
                raise ValueError(f"NEID L2 FITS has no OBJECT: {basename}")
            aliases = _catalog_object_aliases(row)
            identity_status = "catalog_alias"
            if _normalized_object_name(fits_object) not in aliases:
                relative_path = f"FITS/{basename}"
                exception = identity_exceptions.get(relative_path)
                if exception is None:
                    raise ValueError(
                        "NEID L2 OBJECT does not match any catalog alias and "
                        f"has no declared exception: {basename}"
                    )
                if (
                        str(exception["catalog_object_id"])
                        != str(row["OBJECT_ID"])
                        or str(exception["fits_object"]) != fits_object):
                    raise ValueError(
                        f"Header identity exception does not match {basename}."
                    )
                catalog_gaia_ids = set()
                for value in row.values:
                    catalog_gaia_ids.update(_gaia_source_ids(value))
                header_gaia_ids = (
                    _gaia_source_ids(header.get("QOBJECT", ""))
                    | _gaia_source_ids(header.get("SCI-OBJ", ""))
                )
                gaia_source_id = str(exception["gaia_source_id"])
                if (
                        gaia_source_id not in catalog_gaia_ids
                        or gaia_source_id not in header_gaia_ids):
                    raise ValueError(
                        "Header identity exception lacks matching catalog and "
                        f"FITS Gaia provenance: {basename}"
                    )
                used_identity_exceptions.add(relative_path)
                identity_status = "manifest_exception_gaia_verified"
            arrays = []
            for index in (1, 2, 4, 5, 7, 8, 15, 16, 17):
                data = hdus[index].data
                if data is None or np.ndim(data) != 2:
                    raise ValueError(
                        f"NEID L2 FITS HDU {index} is not a 2-D array: {basename}"
                    )
                arrays.append(data)
            if len({array.shape for array in arrays}) != 1:
                raise ValueError(f"NEID L2 FITS array shapes differ: {basename}")
        mapping.append({
            'OBJECT_ID': str(row['OBJECT_ID']),
            'basename': basename,
            'fits_object': fits_object,
            'identity_status': identity_status,
        })
    unused_exceptions = set(identity_exceptions).difference(
        used_identity_exceptions
    )
    if unused_exceptions:
        raise ValueError(
            "Unused header identity exceptions: "
            f"{sorted(unused_exceptions)}"
        )
    actual = {item.name for item in (path / 'FITS').glob('*.fits')}
    if actual != seen_basenames:
        raise ValueError("Catalog-to-FITS mapping is not one-to-one and complete.")
    return {'mapped_objects': mapping, 'fits_count': len(mapping)}


def _safe_extract_zip(
        archive, destination, *, max_members=DEFAULT_MAX_ZIP_MEMBERS,
        max_uncompressed_bytes=DEFAULT_MAX_EXPANDED_BYTES):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    if any(destination.iterdir()):
        raise ValueError("ZIP extraction destination must be empty.")
    with zipfile.ZipFile(archive, 'r') as zip_ref:
        members = zip_ref.infolist()
        if len(members) > int(max_members):
            raise ValueError(
                f"Library archive has {len(members)} members; limit is {max_members}."
            )
        total_size = sum(member.file_size for member in members)
        if total_size > int(max_uncompressed_bytes):
            raise ValueError(
                "Library archive expands beyond the configured size limit."
            )
        names = [member.filename for member in members]
        if len(names) != len(set(names)):
            raise ValueError("Library archive contains duplicate member names.")
        for member in members:
            pure_name = PurePosixPath(member.filename)
            if pure_name.is_absolute() or '..' in pure_name.parts:
                raise ValueError(f"Unsafe path in library archive: {member.filename!r}")
            unix_type = (member.external_attr >> 16) & 0o170000
            if unix_type == stat.S_IFLNK:
                raise ValueError(f"Symlink in library archive: {member.filename!r}")
            if unix_type not in (0, stat.S_IFREG, stat.S_IFDIR):
                raise ValueError(f"Special file in library archive: {member.filename!r}")
            if member.flag_bits & 0x1:
                raise ValueError(f"Encrypted ZIP member is unsupported: {member.filename!r}")
            candidate = (destination / Path(*pure_name.parts)).resolve()
            if destination != candidate and destination not in candidate.parents:
                raise ValueError(f"Unsafe path in library archive: {member.filename!r}")
            if member.is_dir():
                candidate.mkdir(parents=True, exist_ok=True)
                continue
            candidate.parent.mkdir(parents=True, exist_ok=True)
            with zip_ref.open(member, 'r') as source, candidate.open('xb') as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)


def _stream_download(
        url, destination, chunk_size=1024 * 1024, *, timeout=60,
        max_bytes=DEFAULT_MAX_DOWNLOAD_BYTES, expected_bytes=None):
    # MD5 is only an integrity check against the value published by Zenodo;
    # it is not cryptographic source authentication.
    digest = hashlib.md5()
    total = 0
    with urlopen(url, timeout=timeout) as response, Path(destination).open('xb') as stream:
        content_length = response.headers.get('Content-Length')
        expected_length = int(content_length) if content_length is not None else None
        if expected_bytes is not None:
            if expected_length is None:
                raise ValueError(
                    "Authenticated library download requires Content-Length."
                )
            if expected_length != int(expected_bytes):
                raise ValueError(
                    "Library download size does not match the pinned Zenodo record."
                )
        if expected_length is not None and expected_length > int(max_bytes):
            raise ValueError("Library download exceeds the configured size limit.")
        while True:
            chunk = response.read(chunk_size)
            if not chunk:
                break
            total += len(chunk)
            if total > int(max_bytes):
                raise ValueError("Library download exceeded the configured size limit.")
            stream.write(chunk)
            digest.update(chunk)
    if expected_length is not None and total != expected_length:
        raise ValueError(
            f"Incomplete library download: expected {expected_length} bytes, got {total}."
        )
    if expected_bytes is not None and total != int(expected_bytes):
        raise ValueError("Library download byte count does not match the pinned record.")
    return digest.hexdigest()


def _atomic_replace_directory(source, destination):
    """Atomically install ``source``, restoring a prior directory on error."""
    source = Path(source)
    destination = Path(destination)
    _require_regular_path(source, "replacement library directory", directory=True)
    if destination.is_symlink():
        raise ValueError("Library destination cannot be a symbolic link.")
    if destination.exists() and not destination.is_dir():
        raise ValueError("Library destination must be a directory when it exists.")
    if not destination.exists():
        os.replace(source, destination)
        return
    placeholder = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}.backup-", dir=destination.parent
    ))
    placeholder.rmdir()
    backup = placeholder
    os.replace(destination, backup)
    try:
        os.replace(source, destination)
    except BaseException:
        try:
            os.replace(backup, destination)
        except BaseException as restore_error:
            raise RuntimeError(
                f"Library replacement failed and recovery is at {backup}."
            ) from restore_error
        raise
    shutil.rmtree(backup)


def get_library(overwrite=False, library_path=None, verbose=1):
    """Download and verify the public library into user-writable storage.

    Nothing is downloaded at import time.  ``library_path`` overrides
    ``NEIDSPECMATCH_LIBRARY`` and the platform-specific user-data default.
    """
    destination = Path(config.resolve_library_path(library_path))
    release_allowlist, allowlist_sha256 = load_default_library_allowlist()
    if destination.is_symlink():
        raise ValueError("Library destination cannot be a symbolic link.")
    if destination.exists() and not destination.is_dir():
        raise ValueError("Library destination must be a directory when it exists.")
    if destination.exists() and not overwrite:
        result = validate_library(
            destination, expected_fits=78, require_manifest=True, deep=True,
            expected_library_id=config.DEFAULT_LIBRARY_ID,
            catalog_name=config.DEFAULT_LIBRARY_CATALOG,
        )
        result["release_allowlist"] = _verify_library_against_allowlist(
            destination, release_allowlist, allow_installed_manifest=True
        )
        if verbose:
            print('Library already present at {}'.format(destination))
        return result

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
            prefix="neidspecmatch-library-", dir=destination.parent) as temp_name:
        temporary = Path(temp_name)
        archive = temporary / (config.DEFAULT_LIBRARY_ID + ".zip")
        extracted = temporary / "extracted"
        extracted.mkdir()
        if verbose:
            print('Downloading library from {}'.format(config.URL_LIBRARY))
        checksum = _stream_download(
            config.URL_LIBRARY, archive,
            expected_bytes=config.LIBRARY_ZIP_SIZE_BYTES,
        )
        if checksum.lower() != config.LIBRARY_ZIP_MD5.lower():
            raise ValueError(
                "Library download checksum mismatch: expected {}, got {}.".format(
                    config.LIBRARY_ZIP_MD5, checksum)
            )
        _safe_extract_zip(archive, extracted)

        candidates = [extracted]
        candidates.extend(path for path in extracted.iterdir() if path.is_dir())
        source = next(
            (path for path in candidates
             if (path / config.DEFAULT_LIBRARY_CATALOG).is_file()),
            None,
        )
        if source is None:
            raise ValueError("Downloaded archive does not contain the expected catalog.")
        if source != extracted and set(extracted.iterdir()) != {source}:
            raise ValueError(
                "Downloaded library archive must contain one exact library root."
            )
        _verify_library_against_allowlist(source, release_allowlist)
        build_library_manifest(
            source,
            library_id=config.DEFAULT_LIBRARY_ID,
            catalog_name=config.DEFAULT_LIBRARY_CATALOG,
            expected_fits=78,
            source_url=config.URL_LIBRARY,
            source_archive_integrity={
                "algorithm": "md5",
                "value": checksum,
                "scope": "integrity_only_not_authentication",
            },
            release_allowlist={
                "id": release_allowlist["allowlist_id"],
                "sha256": allowlist_sha256,
            },
        )
        validate_library(
            source, expected_fits=78, require_manifest=True, deep=True,
            expected_library_id=config.DEFAULT_LIBRARY_ID,
            catalog_name=config.DEFAULT_LIBRARY_CATALOG,
            validate_fits_schema=True,
        )
        if destination.exists() and not overwrite:
            raise FileExistsError(destination)
        _atomic_replace_directory(source, destination)

    result = validate_library(
        destination, expected_fits=78, require_manifest=True, deep=True,
        expected_library_id=config.DEFAULT_LIBRARY_ID,
        catalog_name=config.DEFAULT_LIBRARY_CATALOG,
    )
    result["release_allowlist"] = _verify_library_against_allowlist(
        destination, release_allowlist, allow_installed_manifest=True
    )
    if verbose:
        print('Library installed at {}'.format(destination))
    return result
