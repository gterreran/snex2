# -*- coding: utf-8 -*-
"""SDSS_search.py

"""

from astroquery.sdss import SDSS
from astropy import units as u
from astropy.time import Time
from astropy.wcs import WCS
import astropy.io.fits as pf
import numpy as np
import logging

logger = logging.getLogger(__name__)

def search_for_SDSS(self, survey):
    '''
    This function searches and downloads SDSS templates.
    Together with the data, it will fetch also the variance images

    '''
    
    # Instruments are sorted by FOV size.
    # :self.instruments[0]: is the instrument with the biggest FOV
    # and therefore requires the biggest template. We start creating
    # this, and instead of repeating the procedure, we just scale down
    # this first template to match the smaller requirements. 

    from . import generate_FOV_grid, LCO_INSTRUMENTS, SURVEYS
    
    inst = self.instruments[0]

    # Generate grid of points covering the FOV
    grid = generate_FOV_grid(_center = self.coord, _fov = LCO_INSTRUMENTS[inst].fov, _step=SURVEYS[survey].skycell_size)

    all_pointings = []

    img_data_list = {}
    weight_data_list = {}

    for f,flt in enumerate(self.filters):

        img_data_list[flt] = []
        weight_data_list[flt] = []

        # For the first filter, we need to collect the references for all the required SDSS fields
        # We don't need to repeat the query for the following filters
        # The SDSS fields references are stored in :all_pointings:
        if f == 0:
            while True:
                # SDSS queries radius must be less than 3.0 arcmin.
                # Initially :grid[0]: will be the center of the field.
                # After each iteration, the :grid: array is trimmed down and its first element will be
                # the point on the grid most distant from the center
                logger.info(f'Querying the SDSS database for coordinates {grid[0].ra.deg} {grid[0].dec.deg}.')
                xid = SDSS.query_region(grid[0], radius=3*u.arcmin, spectro=False)
                if xid is None:
                    logger.info(f'\nThe provided coordinates {grid[0].ra.deg},{grid[0].dec.deg} do not appear to be in the SDSS footprint.\n')
                    return
                
                logger.info('Coordinates in SDSS footprint.')
                new_pointings = [list(el) for el in np.unique(xid['run','camcol','field']) if list(el) not in all_pointings]

                all_pointings = all_pointings + new_pointings

                for run, camcol, field in new_pointings:
                    #fetching the fits as an astropy hdu
                    logger.info(f'Retrieving file with run={run}, camcol={camcol}, field={field} and filter {flt}.')
                    im = SDSS.get_images(run=run, camcol=camcol, field=field, band=flt, cache=True)[0]

                    sky_subtracted_original_image_hdu, weight_image_hdu = elaborate_SDSS(im, flt, run, camcol, field)
                    img_data_list[flt].append(sky_subtracted_original_image_hdu)
                    weight_data_list[flt].append(weight_image_hdu)
                    
            
                for tile in img_data_list[flt]:
                    grid = [p for p in grid if not p.contained_by(WCS(tile.header))]
                
                if len(grid) == 0:
                    break
            
        else:
            #Note that we changed the list to all_pointings now that we have them all
            for run, camcol, field in all_pointings:
                #fetching the fits as an astropy hdu
                im = SDSS.get_images(run=run, camcol=camcol, field=field, band=flt, cache=True)[0]

                sky_subtracted_original_image_hdu, weight_image_hdu = elaborate_SDSS(im, flt, run, camcol, field)
                img_data_list[flt].append(sky_subtracted_original_image_hdu)
                weight_data_list[flt].append(weight_image_hdu)
    
    return img_data_list, weight_data_list



def elaborate_SDSS(_im, _flt, _run, _camcol, _field):
    '''
    SDSS images require a little bit of manipulation.
    This function will produce a ready to be used image
    as well as its corresponding weight image.

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :param _im: SDSS image HDUList
    :type _im: `~astropy.io.fits.HDUList`

    :param _flt: SDSS filter
    :type _flt: string

    :param _run: Run number
    :type _run: int

    :param _camcol: Column in the imaging camera
    :type _camcol: int

    :param _field: Field number
    :type _field: int

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :return: HDU of the image, HDU of the weight image
    :rtype: list

    '''

    from scipy import interpolate

    header = _im[0].header

    #Although we do not write out each single tile, this edits in the header will
    #make it easier to keep track of each tile in the merged fits file
    header['RADESYSa'] = header['RADECSYS']
    del header['RADECSYS']
    header['MJD-OBS'] = Time(f"{header['DATE-OBS']} {header['TAIHMS']}", format='iso', scale='tai').mjd
    header['DATE-OBS'] = 'T'.join(Time(header['MJD-OBS'], format='mjd', scale='utc').iso.split())
    #For some reason the field info is not stored in the header. There is a "FRAME" key, but it's different
    header['FIELD'] = _field

    #retrieve gain and darkVariance from another online table
    q = SDSS.query_sql(f'SELECT DISTINCT gain_{_flt}, darkVariance_{_flt} FROM field WHERE run={_run} AND camcol={_camcol}')
    gain = q[f'gain_{_flt}'][0]
    darkVariance = q[f'darkVariance_{_flt}'][0]

    #im[2] contains the sky
    sky, x, y = _im[2].data[0]
    x_old = np.arange(np.shape(sky)[0])
    y_old = np.arange(np.shape(sky)[1])

    # `interp2d` is deprecated in SciPy 1.10 and will be removed in SciPy 1.14.0
    # se We are using `interpolate.RegularGridInterpolator` instead
    # NOTE: `interp2d` has x and y inverted but doesn't need a `meshgrid` to produce the output
    # This is the code if `interp2d` is still preferred despite of its deprecation
    # tck  = interpolate.interp2d(y_old, x_old, sky, kind=method)
    # sky_image = tck(x, y)

    tck = interpolate.RegularGridInterpolator((x_old,y_old), sky, method="linear", bounds_error=False, fill_value = None)

    # need to invert x and in y to produce the meshgrid with the right shape
    x_grid, y_grid = np.meshgrid(y, x, indexing='ij')

    sky_image = tck((x_grid, y_grid)) #in ADU

    #im[1] contains a row with factors nanomaggies per counts for each column
    calib_image = np.asarray([_im[1].data.tolist()]*np.shape(_im[0].data)[0])

    sky_subtracted_original_image = _im[0].data / calib_image #in ADU

    variance_image = ((_im[0].data / calib_image + sky_image) / gain  + darkVariance)# in ADU * calib_image **2

    sky_subtracted_original_image_hdu = pf.PrimaryHDU(sky_subtracted_original_image, header=header)
    #to reproject and combine, we just need the weight, not the variance
    weight_image_hdu = pf.PrimaryHDU(1/variance_image, header=header)

    return sky_subtracted_original_image_hdu, weight_image_hdu

def track_skycells_in_header(_header, _data):
    _header['NTILES'] = (len(_data), 'Number of SDSS tiles merged')
    for i,im in enumerate(_data, start=1):
        _header[f'RUN{i:>03}'] = (im.header['RUN'], f'Run number of tile {i}')
        _header[f'CAMCO{i:>03}'] = (im.header['CAMCOL'], f'Column in the imaging camera of tile {i}')
        _header[f'FIELD{i:>03}'] = (im.header['FIELD'], f'Field number of tile {i}')
        _header[f'MJD{i:>03}'] = (im.header['MJD-OBS'], f'MJD of tile {i}')
    
    return _header