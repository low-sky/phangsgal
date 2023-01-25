from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from scipy.optimize import newton
from astropy import units as u
from .phangs_galaxies import PhangsGalaxy

__all__ = ['predict_metallicity',
           'predict_alphaCO10',
           'alphaCO_from_coords']

def predict_metallicity(
        Mstar, calibrator='O3N2(PP04)', MZR='Sanchez+19',
        Rgal=None, Re=None, gradient=None, **kwargs):
    """
    Predict log(O/H) based on a given mass-metallicity relation.

    This function predicts gas phase oxygen abundances log(O/H)
    from galaxy stellar mass values, based on a given
    mass-metallicity relation for a particular calibrator.

    This function can also predict a galactocentric-radius-dependent
    metallicity, provided that `Rgal` and `Re` are both given,
    and a particular metallicity gradient is specified.

    Parameters
    ----------
    Mstar : number, `~numpy.ndarray`, `~astropy.units.Quantity` object
        Stellar mass value(s), un units of solar mass
    calibrator : {'O3N2(PP04)', 'R23(T04)', 'R23(R11)', 't2(P12)'}
        Metallicity calibrator to adopt. Default is 'O3N2(PP04)'
    MZR : {'Sanchez+19', 'Sanchez+17', 'Sanchez+13'}
        Mass-metallicity relation to adopt. Default is 'Sanchez+19'
    Rgal : number, ndarray, Quantity object
        Galactocentric radius/radii, in units of kilo-parsec
    Re : number, ndarray, Quantity object
        Effective radius/radii, in units of kilo-parsec
    gradient : {'Sanchez+14', number}
        Metallicity gradient to adopt, in units of dex/Re
        'Sanchez+14' means adopting a fixed gradient of -0.1 dex/Re

    Returns
    -------
    logOH : number or `~numpy.ndarray`
        Predicted gas phase oxygen abundance(s)
    """

    calibrators = ['O3N2(PP04)', 'R23(T04)', 'R23(R11)', 't2(P12)']
    if calibrator not in calibrators:
        raise ValueError(
            "Available choices for `calibrator` are: "
            "{}".format(calibrators))

    if hasattr(Mstar, 'unit'):
        x = np.log10(Mstar.to('Msun').value) - 8
    else:
        x = np.log10(Mstar) - 8

    # Mass-metallicity relation
    if MZR == 'Sanchez+19':  # SAMI
        if calibrator == 'O3N2(PP04)':
            a = 8.73  # +- 0.03
            b = 0.010  # +- 0.002
        elif calibrator == 'R23(T04)':
            a = 8.84  # +- 0.02
            b = 0.007  # +- 0.001
        elif calibrator == 'R23(R11)':
            a = 8.48  # +- 0.02
            b = 0.004  # +- 0.001
        elif calibrator == 't2(P12)':
            a = 8.84  # +- 0.02
            b = 0.008  # +- 0.001
        c = 3.5
        logOH_Re = a + b*(x-c)*np.exp(-(x-c))

    elif MZR == 'Sanchez+17':  # CALIFA
        if calibrator == 'O3N2(PP04)':
            a = 8.76  # +- 0.06
            b = 0.005  # +- 0.037
        elif calibrator == 'R23(T04)':
            a = 8.92  # +- 0.04
            b = 0.008  # +- 0.029
        elif calibrator == 'R23(R11)':
            a = 8.54  # +- 0.03
            b = 0.003  # +- 0.020
        elif calibrator == 't2(P12)':
            a = 8.85  # +- 0.01
            b = 0.007  # +- 0.061
        c = 3.5
        logOH_Re = a + b*(x-c)*np.exp(-(x-c))

    elif MZR == 'Sanchez+13':
        # this calibrator likely underestimates log(O/H)
        # for low mass (<10^9) galaxies
        if calibrator == 'O3N2(PP04)':
            a = 8.74  # +- 0.01
            b = 0.018  # +- 0.007
        else:
            raise ValueError("")
        c = 3.5  # +- 0.3
        logOH_Re = a + b*(x-c)*np.exp(-(x-c))

    else:
        raise ValueError("Invalid input value for `MZR`!")

    if (Rgal is None) or (Re is None) or (gradient is None):
        return logOH_Re

    if hasattr(Mstar, 'unit'):
        Rgal_ = Rgal.to('kpc').value
    else:
        Rgal_ = Rgal
    if hasattr(Mstar, 'unit'):
        Re_ = Re.to('kpc').value
    else:
        Re_ = Re

    # metallicity gradient
    if gradient == 'Sanchez+14':
        if not calibrator == 'O3N2(PP04)':
            print("Beware: "
                  "Sanchez+14 adopted the O3N2(PP04) calibrator")
        alpha_logOH = -0.10  # dex/Re
        logOH = (logOH_Re + alpha_logOH * (Rgal_/Re_ - 1))
    else:
        alpha_logOH = gradient  # dex/Re
        logOH = (logOH_Re + alpha_logOH * (Rgal_/Re_ - 1))

    return logOH


def predict_alphaCO10(
        prescription='constant', iterative=False,
        PHANGS_Zprime=None,
        N12_Zprime=None, N12_WCO10GMC=None,
        B13_Zprime=None, B13_SigmaGMC=None, B13_Sigmatot=None,
        B13_WCO10GMC=None, B13_WCO10kpc=None, B13_Sigmakpc=None,
        A16_logOH_PC09=None,
        A17_logOH_PP04=None, A17_Mstar=None,
        A17_SFR=None, A17_redshift=None,
        suppress_error=False,
        **kwargs):
    """
    Predict CO-to-H2 conversion factor following a given prescription.

    Available prescriptions are either existing suggestions in the
    literature, or motivated by discussions in the X_CO hacking group.
    Detailed explanations on each prescription can be found below.

    + prescription='constant'
        To return the Milky Way conversion factor value
        alpha_CO(1-0) = 4.35 Msun/pc^2/(K*km/s).
        Reference: Bolatto et al. (2013), ARA&A, 51, 207

    + prescription='PHANGS'
        alpha_CO(1-0) = 4.35 * Z'^{-1.6} Msun/pc^2/(K*km/s).
        Requred input parameters:
        1. `PHANGS_Zprime` [unitless]
            metallicity in terms of Z'

    + prescription='Narayanan+12'
        Required input parameters:
        1. `N12_Zprime` [unitless]
            metallicity in terms of Z'
        2. `N12_WCO10GMC` [K*km/s]
            CO(1-0) integrated intensity in GMCs
        Note: an additional multiplicative factor of 1.35 is included
        here to correct for Helium contribution (it is not included
        in the original formula in this paper).
        Reference: Narayanan et al. (2012), MNRAS, 421, 3127 

    + prescription='Bolatto+13', iterative=False
        Required input parameters:
        1. `B13_Zprime` [unitless]
            metallicity in terms of Z'
        2. `B13_SigmaGMC` [Msun/pc^2]
            typical gas surface density of GMCs
        3. `B13_Sigmatot` [Msun/pc^2]
            kpc-scale surface density of total gas + stars
        Reference: Bolatto et al. (2013), ARA&A, 51, 207

    + prescription='Bolatto+13', iterative=True
        Similar to the previous prescription, but to iteratively
        solve for alpha_CO based on observables.
        Required input parameters:
        1. `B13_Zprime` [unitless]
            metallicity in terms of Z'
        2. `B13_WCO10GMC` [K*km/s]
            typical CO(1-0) intensity of GMCs
        3. `B13_WCO10kpc` [K*km/s]
            kpc-scale average CO(1-0) intensity
        4. `B13_Sigmakpc` [Msun/pc^2]
            kpc-scale surface density of HI gas + stars
        Reference: Bolatto et al. (2013), ARA&A, 51, 207

    + prescription='Amorin+16'
        Required input parameters:
        1. `A16_logOH_PC09` [unitless]
            metallicity in terms of PC09 N2 calibrator
            (see Perez-Montero & Contini 2009)
        Reference: Amorin et al. (2016), A&A, 588, A23

    + prescription='Accurso+17'
        Required input parameters:
        1. `A17_logOH_PP04` [unitless]
            metallicity in terms of PP04 O3N2 calibrator
            (see Pettini & Pagel 2004)
        2. `A17_Mstar` [Msun] (optional)
            galaxy global stellar mass
        3. `A17_SFR` [Msun/yr] (optional)
            galaxy global star formation rate
        4. `A17_redshift` [unitless] (optional)
            galaxy redshift
        Note: the secondary dependence on log\Delta(MS) is turned off
        if any of `A17_Mstar`, `A17_SFR`, or `A17_redshift` is None.
        Reference: Accurso et al. (2017), MNRAS, 470, 4750

    Parameters
    ----------
    prescription : string, optional
        The alphaCO prescription to adopt (default='constant')
        See available choices and explanations above in the docstring.
    iterative : bool, optional
        Whether to solve for alphaCO iteratively.
        See descriptions above in the docstring.
    suppress_error : bool, optional
        Whether to suppress error if iteration fail to converge.
        Default is to not suppress.
    [other available keyword arguments]
        See descriptions above in the docstring. Note that these
        inputs can be either scalar or `~numpy.ndarray` object.
        If they do not carry units, a default unit is assumed for
        each of them (see above); otherwise the units they carry
        need to be equivalent to the default unit.
    **kwargs
        Keywords to be passed to `~scipy.optimize.newton`.
        Useful only when using the iterative prescriptions.

    Returns
    -------
    alphaCO : `~astropy.units.Quantity` object
        Predicted CO-to-H2 conversion factor, carrying along a unit of
        Msun/pc^2/(K*km/s).
    """

    alphaCO_MW = 4.35  # Milky Way value

    prescriptions = ['constant', 'PHANGS', 'Narayanan+12',
                     'Bolatto+13', 'Amorin+16', 'Accurso+17']

    if prescription not in prescriptions:
        raise ValueError(
            "Available choices for `prescription` are: "
            "{}".format(prescriptions))

    elif prescription == 'constant':
        alphaCO = np.array(alphaCO_MW)

    elif prescription == 'PHANGS':
        alphaCO = (
            np.array(alphaCO_MW) * np.atleast_1d(PHANGS_Zprime)**-1.6)

    elif prescription == 'Narayanan+12':
        if (N12_Zprime is None) or (N12_WCO10GMC is None):
            raise ValueError(
                "Narayanan+12 prescription needs "
                "both `N12_Zprime` and `N12_WCO10GMC`")
        if hasattr(N12_WCO10GMC, 'unit'):
            WCO10GMC = N12_WCO10GMC.to('K km / s').value
        else:
            WCO10GMC = N12_WCO10GMC
        alphaCO = 10.7 * np.atleast_1d(WCO10GMC)**-0.32
        alphaCO[alphaCO > 6.3] = 6.3
        alphaCO = alphaCO / N12_Zprime**0.65
        alphaCO = alphaCO * 1.35  # include Helium contribution

    elif prescription == 'Bolatto+13':
        if not iterative:
            if ((B13_Zprime is None) or (B13_SigmaGMC is None) or
                    (B13_Sigmatot is None)):
                raise ValueError(
                    "Bolatto+13 prescription (non-iterative) needs "
                    "`B13_Zprime`, `B13_SigmaGMC`, and `B13_Sigmatot`"
                )
            # CO faint correction
            if hasattr(B13_SigmaGMC, 'unit'):
                SigGMC100 = B13_SigmaGMC.to('100 Msun / pc^2').value
            else:
                SigGMC100 = B13_SigmaGMC / 100
            alphaCO = 2.9 * np.exp(0.4 / B13_Zprime / SigGMC100)
            # starburst correction
            if hasattr(B13_Sigmatot, 'unit'):
                Sigtot100 = B13_Sigmatot.to('100 Msun / pc^2').value
            else:
                Sigtot100 = B13_Sigmatot / 100
            f_SB = np.atleast_1d(Sigtot100)**-0.5
            f_SB[f_SB > 1] = 1
            alphaCO = alphaCO * f_SB
        else:
            # iteratively solve for alphaCO based on input observables
            if ((B13_Zprime is None) or (B13_WCO10GMC is None) or
                    (B13_WCO10kpc is None) or (B13_Sigmakpc is None)):
                raise ValueError(
                    "Bolatto+13 prescription (iterative) needs "
                    "`B13_Zprime`, `B13_WCO10GMC`, `B13_WCO10kpc`, "
                    "and `B13_Sigmakpc`")
            if hasattr(B13_WCO10GMC, 'unit'):
                WG = B13_WCO10GMC.to('K km / s').value
            else:
                WG = B13_WCO10GMC
            if hasattr(B13_WCO10kpc, 'unit'):
                WK = B13_WCO10kpc.to('K km / s').value
            else:
                WK = B13_WCO10kpc
            if hasattr(B13_Sigmakpc, 'unit'):
                SK = B13_Sigmakpc.to('Msun / pc2').value
            else:
                SK = B13_Sigmakpc
            Zp, WG, WK, SK = np.broadcast_arrays(
                B13_Zprime, WG, WK, SK)
            alphaCO = []
            x0 = alphaCO_MW
            for Zp_, WG_, WK_, SK_ in zip(
                    Zp.ravel(), WG.ravel(), WK.ravel(), SK.ravel()):
                if not ((Zp_ > 0) and (WG_ > 0) and
                        (WK_ > 0) and (SK_ > 0)):
                    alphaCO += [np.nan]
                    continue

                def func(x):
                    return predict_alphaCO10(
                        prescription='Bolatto+13',
                        iterative=False,
                        B13_Zprime=Zp_,
                        B13_SigmaGMC=WG_*x,
                        B13_Sigmatot=(WK_*x+SK_)).value - x
                if suppress_error:
                    try:
                        alphaCO += [newton(func, x0, **kwargs)]
                    except RuntimeError as e:
                        print(e)
                        alphaCO += [np.nan]
                else:
                    alphaCO += [newton(func, x0, **kwargs)]
            alphaCO = np.array(alphaCO).reshape(Zp.shape)

    elif prescription == 'Amorin+16':
        if A16_logOH_PC09 is None:
            raise ValueError(
                "Amorin+16 prescription needs `A16_logOH_PC09`")
        alphaCO = 10 ** (
            0.68 - 1.45*np.atleast_1d(A16_logOH_PC09 - 8.7))

    elif prescription == 'Accurso+17':
        if A17_logOH_PP04 is None:
            raise ValueError(
                "Accurso+17 prescription needs at least "
                "`A17_logOH_PP04`")
        alphaCO = 10 ** (
            14.752 - 1.623*np.atleast_1d(A17_logOH_PP04))
        if ((A17_Mstar is not None) and (A17_SFR is not None) and
                (A17_redshift is not None)):
            if hasattr(A17_Mstar, 'unit'):
                Mstar = A17_Mstar.to('Msun').value
            else:
                Mstar = A17_Mstar
            if hasattr(A17_SFR, 'unit'):
                SFR = A17_SFR.to('Msun/yr').value
            else:
                SFR = A17_SFR
            log_sSFR = np.log10(SFR * 1e9 / Mstar)
            log_sSFR_MS = (
                -1.12 + 1.14*A17_redshift - 0.19*A17_redshift**2 -
                (0.3 + 0.13*A17_redshift) * (np.log10(Mstar) - 10.5))
            log_Delta_MS = log_sSFR - log_sSFR_MS
            alphaCO = alphaCO * 10 ** (0.062 * log_Delta_MS)

    if alphaCO.size == 1:
        alphaCO = alphaCO.item()
    return alphaCO * u.Unit('Msun s / (pc2 K km)')


def alphaCO_from_coords(galaxy, skycoord=None, ra=None, dec=None,
                        header=None, logOH=None, prescription='Accurso+17',
                        line_ratio=0.65, A17_logOH_PP04=None, 
                        A17_redshift=None, A17_SFR=None, A17_Mstar=None,
                        gradient='Sanchez+14', PHANGS_Zprime=None,
                        B13_Zprime=None,
                         **kwargs):

    """
    Predicts alphaCO(1-0) from coordinates given a PHANGS galaxy object

    Parameters
    ----------
    galaxy : phangs.PhangsGalaxy
        PhangsGalaxy object initialized for the galaxy to return coordinates

    Keywords
    --------
    skycoord : astropy.SkyCoord
        SkyCoord object(s) from astropy from which the values of alpha_CO are generated.
    ra : np.array
        Right ascension of samples
    dec : np.array
        Declination of samples
    header : astropy.fits.Header
        FITS header specifying the map for the alpha_CO values to be returned
    line_ratio : float
        Line ratio of the CO line relative to CO(1-0) used for rescaling.  Assumes 0.65.
    Extra keywords are passed to predict_alphaCO10() or predict_metallicity()
    """
    
    radius = galaxy.radius(skycoord=skycoord, ra=ra, dec=dec, header=header)
    sfr = galaxy.sfr

    if A17_SFR is None:
        A17_SFR = sfr

    if logOH is None:
        Reff = (galaxy.Re / u.rad * galaxy.distance).to(u.kpc)
        logOH = predict_metallicity(galaxy.m_star, Rgal=radius, Re=Reff,
                                    gradient=gradient, **kwargs)
    if A17_logOH_PP04 is None:
        A17_logOH_PP04 = logOH
    if PHANGS_Zprime is None:
        PHANGS_Zprime = 1e1**(logOH - 8.69)
    if A17_Mstar is None:
        A17_Mstar = galaxy.m_star
    if A17_redshift is None:
        A17_redshift = 0
    alphaCO = predict_alphaCO10(prescription=prescription, 
                                A17_logOH_PP04=A17_logOH_PP04,
                                PHANGS_Zprime=PHANGS_Zprime,
                                A17_Mstar=A17_Mstar,
                                A17_SFR=sfr,
                                A17_redshift=A17_redshift, **kwargs) / line_ratio
    return(alphaCO)
    
# import phangs
# g = phangs.PhangsGalaxy('ngc3137')
# from astropy.io import fits
# hdr = fits.getheader('ngc3137_12m+7m+tp_co21_pbcorr_round_k.fits')
# aco = phangs.XCO.alphaCO_from_coords(g, header=hdr)
