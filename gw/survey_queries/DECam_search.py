# -*- coding: utf-8 -*-
"""DECam_search.py

"""

import urllib.request
from numpy.core.defchararray import startswith
from astropy import units as u
from pyvo.dal import sia
#from astroquery.hips2fits import hips2fits
#from astropy import wcs as astropy_wcs
import astropy.io.fits as pf

# List of available services at:
# https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/DataAccessInterfaces/SimpleImageAccessSIA/SimpleImageAccessSIA.html
#COLLECTION = 'coadd_all'
COLLECTION = 'ls_dr9'
COLLECTION = 'coadd/decaps_dr2'
COLLECTION = 'delve_dr2'
DECAM_SERVICE = f'https://datalab.noirlab.edu/sia/{COLLECTION}'


#url = 'https://alasky.cds.unistra.fr/hips-image-services/hips2fits#hips=CDS%2FP%2FDECaPS%2FDR2%2Fg&ra=119.51849999999999&dec=-27.298400000000004&fov=0.4&projection=TAN&format=fits'
#url = 'https://alasky.cds.unistra.fr/hips-image-services/hips2fits?hips=CDS%2FP%2FDECaPS%2FDR2%2Fg&width=1200&height=900&fov=0.4&projection=TAN&coordsys=icrs&rotation_angle=0.0&ra=119.51849999999999&dec=-27.298400000000004&format=fits'
#url = 'https://alasky.cds.unistra.fr/hips-image-services/hips2fits?hips=CDS%2FP%2FPanSTARRS%2FDR1%2Fg&width=1200&height=900&fov=0.4&projection=TAN&coordsys=icrs&rotation_angle=0.0&ra=48.837009&dec=41.611606&format=fits'
#urllib.request.urlretrieve(url, 'test.fits')

# hips2fits.timeout = 180

# w = astropy_wcs.WCS(header={
#     'NAXIS1': 4100,         # Width of the output fits/image
#     'NAXIS2': 4100,         # Height of the output fits/image
#     'WCSAXES': 2,           # Number of coordinate axes
#     'CRPIX1': 2050,       # Pixel coordinate of reference point
#     'CRPIX2': 2050,        # Pixel coordinate of reference point
#     'CDELT1': -0.000108056,        # [deg] Coordinate increment at reference point
#     'CDELT2': 0.000108056,         # [deg] Coordinate increment at reference point
#     'CUNIT1': 'deg',        # Units of coordinate increment and value
#     'CUNIT2': 'deg',        # Units of coordinate increment and value
#     'CTYPE1': 'RA---TAN',   # galactic longitude, Mollweide's projection
#     'CTYPE2': 'DEC--TAN',   # galactic latitude, Mollweide's projection
#     'CRVAL1': 48.837009,          # [deg] Coordinate value at reference point
#     'CRVAL2': 41.611606,          # [deg] Coordinate value at reference point
# })
# hips = 'CDS/P/PanSTARRS/DR1/g'
# result = hips2fits.query_with_wcs(
#    hips=hips,
#    wcs=w
# )

# result.writeto('test.fits', overwrite=True)

def search_for_DECam(self, survey):

    from . import generate_FOV_grid, LCO_INSTRUMENTS, SURVEYS

    connect = sia.SIAService(DECAM_SERVICE)
    table = connect.search(pos = (self.coord.ra.deg, self.coord.dec.deg), size = (LCO_INSTRUMENTS['sinistro'].fov.to(u.deg).value, LCO_INSTRUMENTS['sinistro'].fov.to(u.deg).value), verbosity=2).to_table()

    with open('table.csv','w') as out:
        for c in table.columns:
            out.write(c+',')
        out.write('\n')
        for i in range(len(table[table.columns[0]])):
            for c in table.columns:
                out.write(str(table[c][i])+',')
            out.write('\n')
    
    exit()


    for flt in self.filters:
        
        sel = (table['prodtype'] == 'image') & (startswith(table['obs_bandpass'].astype(str), flt))

        if len(table[sel]) != 0:
            out_name = f'{OUT_FOLDER}/{survey}/{self.obj}_{flt}.fits'
            urllib.request.urlretrieve(table[sel]['access_url'][0], out_name)
            self.templates_paths[f'{survey}_{flt}'] = out_name
            with pf.open(out_name, mode='update') as hdu:
                hdu[0].header['MJD-OBS'] = hdu[0].header['MJD_MEAN']
                hdu[0].header['DATE-OBS'] = hdu[0].header['DATEOBS']
                hdu[0].header['DAY-OBS'] = hdu[0].header['DATE-OBS'].split('T')[0]
                hdu[0].header['FILTER'] = flt+'p' #note, this works only with griz
                hdu[0].header['OBJECT'] = self.obj
                hdu[0].header['WCSERR'] = 0
                hdu[0].header['RA'] = self.coord.ra.deg
                hdu[0].header['DEC'] = self.coord.dec.deg

                hdu.flush()
        else:
            print(f'Coordinates {self.coord.ra.deg} {self.coord.dec.deg} not in DECam with filter {flt}.')