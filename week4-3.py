# Import the packages
from astropy.io import fits
from astropy.table import QTable, Table, Column
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Open the fits table
hdul = fits.open('week4.fits')

# Copy this into an astropy table, and split out some of the abundances into individual columns.
# Convert parallax to distance modulus as well (the log10 will choke on some numbers so it is good to do it here)
# Then, we drop the columns that we don’t intend to use.
tdata = Table(hdul[1].data, copy=False)
tdata.add_column(tdata['X_H'][:, np.where(hdul[3].data['ELEM_SYMBOL'] == 'O')[1][0]], name='O_H')
tdata.add_column(tdata['X_H_ERR'][:, np.where(hdul[3].data['ELEM_SYMBOL'] == 'O')[1][0]], name='O_H_ERR')
tdata.add_column(tdata['X_H'][:, np.where(hdul[3].data['ELEM_SYMBOL'] == 'C')[1][0]], name='C_H')
tdata.add_column(tdata['X_H_ERR'][:, np.where(hdul[3].data['ELEM_SYMBOL'] == 'C')[1][0]], name='C_H_ERR')
tdata.add_column(tdata['X_H'][:, np.where(hdul[3].data['ELEM_SYMBOL'] == 'Fe')[1][0]], name='Fe_H')
tdata.add_column(tdata['X_H_ERR'][:, np.where(hdul[3].data['ELEM_SYMBOL'] == 'Fe')[1][0]], name='Fe_H_ERR')
tdata.add_column(tdata['X_H'][:, np.where(hdul[3].data['ELEM_SYMBOL'] == 'Ni')[1][0]], name='Ni_H')
tdata.add_column(tdata['X_H_ERR'][:, np.where(hdul[3].data['ELEM_SYMBOL'] == 'Ni')[1][0]], name='Ni_H_ERR')

tdata.add_column(5.0 * np.log10(1000.0 / tdata['GAIAEDR3_PARALLAX']) - 5.0, name='GAIAEDR3_DISTANCE_MOD')

tdata.keep_columns(['APOGEE_ID', 'RA', 'DEC', 'J', 'J_ERR', 'H', 'H_ERR', 'K', 'K_ERR', 'TEFF',
                    'TEFF_ERR', 'LOGG', 'LOGG_ERR', 'M_H', 'M_H_ERR', 'ALPHA_M', 'ALPHA_M_ERR',
                    'GAIAEDR3_SOURCE_ID', 'GAIAEDR3_PARALLAX', 'GAIAEDR3_PARALLAX_ERROR',
                    'GAIAEDR3_PMRA', 'GAIAEDR3_PMRA_ERROR', 'GAIAEDR3_PMDEC', 'GAIAEDR3_PMDEC_ERROR',
                    'GAIAEDR3_DR2_RADIAL_VELOCITY', 'GAIAEDR3_DR2_RADIAL_VELOCITY_ERROR',
                    'GAIAEDR3_PHOT_G_MEAN_MAG', 'GAIAEDR3_PHOT_BP_MEAN_MAG', 'GAIAEDR3_PHOT_RP_MEAN_MAG',
                    'O_H', 'O_H_ERR', 'C_H', 'C_H_ERR', 'Fe_H', 'Fe_H_ERR', 'Ni_H', 'Ni_H_ERR',
                    'GAIAEDR3_DISTANCE_MOD'])

# Convert this to a pandas dataframe and drop all rows that have any NaN value.
pdata = tdata.to_pandas()
pdata.dropna(inplace=True)

# Write the dataframe into a file
t2 = Table.from_pandas(pdata)
t2.write('APOGEE_AllStarLite_Filtered.fits', format='fits', overwrite=True)