from spectral_cube import SpectralCube, Projection
import os
import astropy.units as u
from astropy.io import fits
from astropy.wcs import wcs
from radio_beam import Beam
from .phangs_galaxies import PhangsGalaxy

__all__ = ['convolve_to_resolution']

def convolve_to_resolution(filename, name=None,
                           resolution=100 * u.pc, **kwargs):
    """
    Convolves an FITS image to a given linear resolution from a given distance sample table
    """    
    
    hdr = fits.getheader(filename)  
    w = wcs.WCS(hdr)
    if w.naxis == 3:
        data = SpectralCube.read(filename)
    if w.naxis == 2:
        data = Projection.from_hdu(fits.open(filename)[0])
    gal = PhangsGalaxy(data.header['OBJECT'])
    target_resn = (resolution / gal.distance * u.rad).to(u.arcsec)
    target_beam = Beam(major=target_resn,
                       minor=target_resn,
                       pa=0 * u.deg)
    data_out = data.convolve_to(target_beam)
    if name is None:
        resn_str = '_' + resolution.to_string().replace(' ','')
        name = filename.replace('.fits', resn_str+'.fits')
    data_out.write(name, overwrite=True)
    
