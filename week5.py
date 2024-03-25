# Import the packages:
from astropy.io import fits
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import CartesianRepresentation, CartesianDifferential

# Open the file
hdul = fits.open('week5.fits.gz')

# Identify the ‘good’ rows that have finite values in all columns:
frows = np.logical_and.reduce((
    np.isfinite(hdul[1].data['parallax']),
    np.isfinite(hdul[1].data['pmra']),
    np.isfinite(hdul[1].data['pmdec']),
    np.isfinite(hdul[1].data['radial_velocity']),
    np.isfinite(hdul[1].data['phot_g_mean_mag']),
    np.isfinite(hdul[1].data['phot_bp_mean_mag']),
    np.isfinite(hdul[1].data['phot_rp_mean_mag'])
))

# Import only the ‘good’ data into the SkyCoord class, and set units.
eq_coords = SkyCoord(frame='icrs',
    ra=hdul[1].data[frows]['ra'] * u.deg,
    dec=hdul[1].data[frows]['dec'] * u.deg,
    distance=(1000.0 / hdul[1].data[frows]['parallax']) * u.parsec,
    pm_ra_cosdec=hdul[1].data[frows]['pmra'] * u.mas / u.yr,
    pm_dec=hdul[1].data[frows]['pmdec'] * u.mas / u.yr,
    radial_velocity=hdul[1].data[frows]['radial_velocity'] * u.km / u.s
)

# Calculate the absolute magnitude and color index:
absmag = hdul[1].data[frows]['phot_g_mean_mag'] + 5.0 - (5.0 * np.log10(eq_coords.distance.to_value(u.pc)))
colour = hdul[1].data[frows]['phot_bp_mean_mag'] - hdul[1].data[frows]['phot_rp_mean_mag']

# Apply a coordinate transformation to galactic position coordinates, and Cartesian velocity coordinates:
gal_coords = eq_coords.galactic
gal_coords.differential_type = CartesianDifferential