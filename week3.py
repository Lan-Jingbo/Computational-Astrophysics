from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math

hdul = fits.open('week4.fits')
temp = hdul[1].data
plt.scatter('TEFF', 'TEFF')
# ra_srad = hdul[1].data['ra'] * (math.pi / 180.0) * np.cos(hdul[1].data['dec'] * (math.pi / 180.0))
# dec_srad = hdul[1].data['dec'] * (math.pi / 180.0)
# plt.hist2d(ra_srad, dec_srad, bins=(180, 90))

# distance = 1000.0 / hdul[1].data['parallax']
# absmag = hdul[1].data['phot_g_mean_mag'] + 5.0 - (5.0 * np.log10(distance))

# ax = plt.gca()
# plt.scatter(distance, absmag)

# def plot_lf(absmag, distance, dmin, dmax, mbin):
#     farr = np.logical_and(distance >= dmin, distance <= dmax)
#     data = absmag[farr]
#     mmin = math.floor(min(data) / mbin) * mbin
#     mmax = math.ceil(max(data) / mbin) * mbin
#     magbins = np.arange(mmin, mmax, mbin)
#     fm = np.zeros(int(magbins.shape[0]))
#     for ms in magbins:
#         marr = np.logical_and(data >= ms, data < ms + mbin)
#         idx = int(math.floor((ms + mbin * 0.5 - mmin) / mbin))
#         fm[idx] = np.count_nonzero(marr)
#     plt.semilogy(magbins + mbin * 0.5, fm)