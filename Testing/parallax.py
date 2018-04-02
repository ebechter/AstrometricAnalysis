---

import numpy as np

from astropy.time import Time

import ephem as ep

import de421

import jplephem

def get_parallax_offset(ra,dec,parallax,epoch,ra_units='detector'):

    """

    Takes the star's ra, dec, parallax, and a given epoch.

    Returns the parallax shift in star's position, in arcseconds (sexagesimal seconds for raoff; decimal seconds for decoff)

    INPUTS:

    ra          -- RA position of star, in degrees

    dec         -- DEC position of star, in degrees

    parallax    -- parallax of star, in MILLIarcseconds (1000.0*1/distance_pc)

    epoch       -- epoch (decimal years) to compute parallax offset (scalar or monotonically increasing vector)

    ra_units    -- USE EITHER:

                   'time': raoff is in seconds of time in RA (need to multiply by 15 to get same units as decoff)

                   'detector': [default] 1 unit displacement same in both horizontal/vertical direction

                               (factor of cos(dec) applied to raoff; multiply by 15 to get same units)

    OUTPUTS:

    raoff       -- parallax shift in star's position, in milliarcseconds

    decoff      -- parallax shift in star's position, in milliarcseconds

    """

    # Convert from deg to rad

    ra_rad=np.radians(ra)

    dec_rad=np.radians(dec)



    # Use jplephem with de421 to determine position of Earth geocenter wrt SS Barycenter

    eph=jplephem.Ephemeris(de421)

    t=Time(epoch,format='jyear',scale='ut1')

    JD=t.jd

    barycenter = eph.position('earthmoon', JD)

    moonvector = eph.position('moon', JD)

    earthPos = (barycenter - moonvector * eph.earth_share) / 1.49597871e+8

    usex=np.float(earthPos[0])

    usey=np.float(earthPos[1])

    usez=np.float(earthPos[2])



    # Compute offsets (in milliarcseconds)

    delta_alpha = parallax * 1./np.cos(dec_rad) * (usex*np.sin(ra_rad) - usey*np.cos(ra_rad))

    delta_delta = parallax * (usex*np.cos(ra_rad)*np.sin(dec_rad) + usey*np.sin(ra_rad)* np.sin(dec_rad) - usez*np.cos(dec_rad))

    # Convert RA seconds as necessary

    if ra_units=='time':

        delta_alpha = delta_alpha / 15.0

    elif ra_units=='detector':

        delta_alpha = delta_alpha * np.cos(dec_rad)

    else:

        raise Exception('ra_units set to unhandled option: {}'.format(ra_units))



    # Set outputs

    raoff=delta_alpha

    decoff=delta_delta



    return (raoff,decoff)