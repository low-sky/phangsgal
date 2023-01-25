from astropy.coordinates import SkyCoord, Angle, FK5
import astropy.units as u
from astropy.units import Quantity
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import warnings
import numpy as np
from astropy.utils.data import get_pkg_data_filename
import os
import glob
import astropy.utils.introspection as intro


__all__ = ['PhangsGalaxy']

def _parse_galtable(galobj, name):
    # TODO: Generalize
    if os.environ.get('PHANGSDATA') is not None:
        table_dir = os.environ.get('PHANGSDATA') + '/Archive/Products/sample_tables/'
        fl = glob.glob(table_dir + '*.fits')
        # TODO: This isn't right
        table_name = fl[-1]
    else:
        parent = os.path.dirname(intro.resolve_name('phangs').__file__)
        table_name = get_pkg_data_filename('data/phangs_sample_table_v1p6.fits',
                                        package='phangs')
    galtable = Table.read(table_name)
    hits = [x for x in galtable if name.lower() in x['alias']]
    # import pdb; pdb.set_trace()
    if len(hits) > 0:
        if len(hits) > 1:
            exact = np.zeros(len(hits), dtype=np.bool)
            for i, h in enumerate(hits):
                exact[i] = np.any(
                    [n.strip() == name for n in h['alias'].split(';')])
            if np.sum(exact) == 1:
                thisobj = hits[(np.where(exact))[0][0]]
            else:
                raise Exception("More than one exact match in galbase")
        else:
            thisobj = hits[0]
        galobj.name = thisobj['name'].strip()
        galobj.center_position = SkyCoord(
            thisobj['orient_ra'], thisobj['orient_dec'], frame='fk5',
            unit='degree')

        propdict = {'vsys':'orient_vlsr',
                    'distance':'dist',
                    'inclination':'orient_incl',
                    'position_angle':'orient_posang',
                    'mw_av':'mwext_sf11',
                    'm_star':'props_mstar',
                    'm_hi':'props_mhi',
                    'sfr':'props_sfr',
                    'Re':'size_reff',
                    'R25':'size_r25'}
        for tag in propdict:
            setattr(galobj, tag, (thisobj[propdict[tag]]
                                  * galtable[propdict[tag]].unit)) 
        
        galobj.provenance = 'PhangsTable'
        galobj.morph = thisobj['morph_string']
        galobj.has_alma = bool(thisobj['survey_alma_status'])
        galobj.has_astrosat = bool(thisobj['survey_astrosat_status'])
        galobj.has_dense = bool(thisobj['survey_dense_status'])
        galobj.has_muse = bool(thisobj['survey_muse_status'])
        galobj.has_galex = bool(thisobj['survey_galex_status'])
        galobj.has_halpha = bool(thisobj['survey_halpha_status'])
        galobj.has_herschel = bool(thisobj['survey_herschel_status'])
        galobj.has_hi = bool(thisobj['survey_hi_status'])
        galobj.has_hst = bool(thisobj['survey_hst_status'])
        galobj.has_irac = bool(thisobj['survey_irac_status'])
        galobj.morph_T = thisobj['morph_t']
        galobj.morph_bar = bool(thisobj['morph_bar'])

        galobj.table = thisobj
        return True

class PhangsGalaxy(object):
    '''
    Parameters
    ----------
    name : str
        Name of the galaxy.the
    params : dict, optional
        Optionally provide custom parameter values as a dictionary.
    '''

    def __init__(self, name, params=None):

        self.name = name
# An astropy coordinates structure
        self.center_position = None
# This is the preferred name in a database.
        self.canonical_name = None
# With units
        self.distance = None
        self.inclination = None
        self.position_angle = None
        self.redshift = None
        self.vsys = None
        self.PA_is_kinematic = False
        self.provenance = None
        
        if params is not None:
            if not isinstance(params, dict):
                raise TypeError("params must be a dictionary.")

            required_params = ["center_position", "distance", "inclination",
                               "position_angle", "vsys"]
            optional_params = ["canonical_name", "redshift"]

            keys = params.keys()
            for par in required_params:
                if par not in keys:
                    raise ValueError("params is missing the required key"
                                     " {}".format(par))
                setattr(self, par, params[par])

            for par in optional_params:
                if par in keys:
                    setattr(self, par, params[par])
        else:
            if not _parse_galtable(self, name):
                
                try:
                    from astroquery.ned import Ned
                    t = Ned.query_object(name)
                    if len(t) == 1:
                        self.canonical_name = t['Object Name'][0]
                        self.velocity = t['Velocity'][0] * u.km / u.s
                        self.center_position = \
                            SkyCoord(t['RA(deg)'][0], t['DEC(deg)'][0],
                                     frame='fk5',
                                     unit='degree')
                        self.redshift = t['Redshift'][0]
                        self.provenance = 'NED'
                except:
                    warnings.warn("Unsuccessful query to NED")
                    pass

    def __repr__(self):
        return "Galaxy {0} at RA={1}, DEC={2}".format(self.name,
                                                      self.center_position.ra,
                                                      self.center_position.dec)

    def skycoord_grid(self, header=None, wcs=None):
        '''
        Return a grid of RA and Dec values.
        '''
        if header is not None:
            w = WCS(header)
        elif wcs is not None:
            w = wcs
        else:
            raise ValueError("header or wcs must be given.")
        w = WCS(header)
        naxis = w.pixel_shape
        ymat, xmat = np.indices((naxis[1], naxis[0]))
        ramat, decmat = w.celestial.wcs_pix2world(xmat, ymat, 0)
        return SkyCoord(ramat, decmat, unit=(u.deg, u.deg))

    def radius(self, skycoord=None, ra=None, dec=None,
               header=None, returnXY=False):
        if skycoord:
            PAs = self.center_position.position_angle(skycoord)
            Offsets = skycoord
        elif isinstance(header, fits.Header):
            Offsets = self.skycoord_grid(header=header)
            PAs = self.center_position.position_angle(Offsets)
        elif np.any(ra) and np.any(dec):
            Offsets = SkyCoord(ra, dec, unit=(u.deg, u.deg))
            PAs = self.center_position.position_angle(Offsets)
        else:
            warnings.warn('You must specify either RA/DEC, a header or a '
                          'skycoord')
        GalPA = PAs - self.position_angle
        GCDist = Offsets.separation(self.center_position)
        # Transform into galaxy plane
        Rplane = self.distance * np.tan(GCDist)
        Xplane = Rplane * np.cos(GalPA)
        Yplane = Rplane * np.sin(GalPA)
        Xgal = Xplane
        Ygal = Yplane / np.cos(self.inclination)
        Rgal = np.sqrt(Xgal**2 + Ygal**2)
        if returnXY:
            return (Xgal.to(u.pc), Ygal.to(u.pc))
        else:
            return Rgal.to(u.pc)

    def position_angles(self, skycoord=None, ra=None, dec=None,
                        header=None):
        X, Y = self.radius(skycoord=skycoord, ra=ra, dec=dec,
                           header=header, returnXY=True)

        return Angle(np.arctan2(Y, X))

    def to_center_position_pixel(self, wcs=None, header=None):

        if header is not None:
            wcs = WCS(header)

        if wcs is None:
            raise ValueError("Either wcs or header must be given.")

        return self.center_position.to_pixel(wcs)

    def rotation_curve(self, radius, return_bounds=False,
                       extend=False):
        if type(radius) is not Quantity:
            warnings.warn('Radius must be specified with astropy.unit')
            return(None)
        table_name = get_pkg_data_filename('data/RCtable_Nov2019.fits',
                                           package='phangs')
        rctable = Table.read(table_name)
        idx = (rctable['Galaxy'] == self.name.lower())

        if ~np.any(idx):
            warnings.warn('No rotation curve data found for ' + self.name)
            return(None)


        dist_table_name = get_pkg_data_filename('data/RC_assumed_distance.txt',
                                                package='phangs')
        dist_table = Table.read(dist_table_name, format='ascii')
        dist_idx = (dist_table['ID'] == self.name.lower())
        if ~np.any(dist_idx):
            warnings.warn('No distance foundf for ' + self.name)
            warnings.warn('Not rescaling distances in rotation curve')
            distance_factor = 1.0
        else:
            distance_factor = np.float(self.distance.to(u.Mpc).value
                                       / dist_table[dist_idx]['Distance'].data[0])
        rsample = rctable[idx]['Radius'] * distance_factor
        vsample = rctable[idx]['Vrot']
        vrot_lower = rctable[idx]['Vrot_lower']
        vrot_upper = rctable[idx]['Vrot_upper']

        if radius.unit.is_equivalent(u.m):
            r = radius.to(rsample.unit).value
        elif radius.unit.is_equivalent(u.rad):
            r = (radius.to(u.rad) * self.distance.to(rsample.unit)).value
        else:
            warnings.warn("Radius units must be equivalent to angle or distance")
        if extend:
            extendval = vsample.data[-1]
        else:
            extendval = np.nan
        vinterp = np.interp(r, np.r_[0, rsample.data], np.r_[0, vsample.data],
                            left=np.nan, right=extendval)
        if return_bounds:
            if extend:
                extendval_lower = vrot_lower.data[-1]
                extendval_upper = vrot_upper.data[-1]
            else:
                extendval_lower = np.nan
                extendval_upper = np.nan
            vr_lower = np.interp(r, np.r_[0, rsample.data],
                                 np.r_[0, vrot_lower.data],
                                 left=np.nan, right=extendval_lower)
            vr_upper = np.interp(r, np.r_[0, rsample.data],
                                 np.r_[0, vrot_upper.data],
                                 left=np.nan, right=extendval_upper)
            return(vinterp * vsample.unit,
                   vr_lower * vsample.unit,
                   vr_upper * vsample.unit)
        else:
            return(vinterp * vsample.unit)


    def los_velocity(self, skycoord=None, ra=None, dec=None,
                    header=None, extend=False):
            radius = self.radius(skycoord=skycoord, ra=ra, dec=dec, header=header)
            x, y = self.radius(skycoord=skycoord, ra=ra, dec=dec, header=header,
                            returnXY=True)
            phi = np.arctan2(y, x)
            vrot = self.rotation_curve(radius, extend=extend)
            vlos = self.vsys + vrot * np.cos(phi) * np.sin(self.inclination)
            return(vlos)
