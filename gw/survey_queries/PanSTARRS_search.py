# -*- coding: utf-8 -*-
"""PanSTARRS_search.py

"""

from astropy.table import Table
from astropy import units as u
from astropy.time import Time
from astropy.wcs import WCS
import astropy.io.fits as pf

# PS1 urls
PS1_TAB = 'https://ps1images.stsci.edu/cgi-bin/ps1filenames.py'

# PS1 provides a cutout API. We probably prefer to use the whole
# skycell, but we can easily change that with this bool variable
SKYCELL = True
if SKYCELL:
    PS1_CUT = 'https://ps1images.stsci.edu/'
else:
    PS1_CUT = 'https://ps1images.stsci.edu/cgi-bin/fitscut.cgi'


def search_for_PS1(query, survey):
    '''
    This function searches and downloads PS1 templates.
    Together with the data, it will fetch also the weight images (i.e. 1/variance)
    and the mask images.

    '''

    from . import generate_FOV_grid, LCO_INSTRUMENTS, SURVEYS

    # Instruments are sorted by FOV size.
    # :query.instruments[0]: is the instrument with the biggest FOV
    # and therefore requires the biggest template. We start creating
    # this, and instead of repeating the procedure, we just scale down
    # this first template to match the smaller requirements. 
    inst = query.instruments[0]

    # Generate grid of points covering the FOV
    grid = generate_FOV_grid(_center = query.coord, _fov = LCO_INSTRUMENTS[inst].fov, _step=SURVEYS[survey].skycell_size)

    # NOTE: wt are actually variance maps
    data_lists = {'stack':{}, 'mask':{}, 'wt':{}}

    while True:
        # Initially :grid[0]: will be the center of the field.
        # After each iteration, the :grid: array is trimmed down
        # and its first element will be the point on the grid
        # most distant from the center

        # Creating the url to query the PS1 database.
        # NOTE: PS1 needs the size in pixels, considering 0.25 arcsec/pixel.
        # This first query returns a Table with a list of files
        print(f'Querying the PS1 database for coordinates {grid[0].ra.deg} {grid[0].dec.deg}.')
        tableurl = f'{PS1_TAB}?ra={grid[0].ra.deg}&dec={grid[0].dec.deg}&size={LCO_INSTRUMENTS[inst].fov.to(u.arcsec).value/0.25}&format=fits&filters={query.filters}&type=stack,stack.wt,stack.mask'
        table = Table.read(tableurl, format='ascii')

        # If the table is empty, the images are not in PS1
        if len(table) == 0:
            print(f'\nThe provided coordinates {grid[0].ra.deg},{grid[0].dec.deg} do not appear to be in the PS1 footprint.\n')
            return
        
        print('Coordinates in PS1 footprint.')
        for tab in table:
            flt = tab['filter']
            filename = tab['filename']
            filetype = tab['type'].split('.')[-1]

            # Initializing dictionaries
            if flt not in data_lists[filetype]:
                data_lists[filetype][flt]=[]
                data_lists[filetype][flt]=[]
                data_lists[filetype][flt]=[]
            
            if SKYCELL:
                url = f'{PS1_CUT}{filename}'
            else:
                url = f'{PS1_CUT}?ra={grid[0].ra.deg}&dec={grid[0].dec.deg}&size=6000&format=fits&red={filename}'

            print(f'Retrieving file {filename}')
            with pf.open(url) as hdulist:
                # The skycell headers need some adjustments
                if SKYCELL:
                    #hdulist[1].header['TIMESYS'] = 'TAI'
                    hdulist[1].header['RADESYSa'] = 'FK5'
                    hdulist[1].header['DATE-OBS'] = 'T'.join(Time(hdulist[1].header['MJD-OBS'], format='mjd', scale='utc').iso.split())
                    hdulist[1].header['PC1_1'] = hdulist[1].header['PC001001']
                    hdulist[1].header['PC1_2'] = hdulist[1].header['PC001002']
                    hdulist[1].header['PC2_1'] = hdulist[1].header['PC002001']
                    hdulist[1].header['PC2_2'] = hdulist[1].header['PC002002']
                    del hdulist[1].header['PC001001']
                    del hdulist[1].header['PC001002']
                    del hdulist[1].header['PC002001']
                    del hdulist[1].header['PC002002']
                    
                    # The algorithm for decompressing fits file data
                    # doesn't work on PS1 skycells since astropy version 5.3
                    # The easiest way to fix this, is to remove the
                    # 'ZBLANK' key (assuming it's equal to 'BLANK')
                    if 'ZBLANK' in hdulist[1].header:
                        del hdulist[1].header['ZBLANK']

                    # The image flux scaling is also non-standard for the stack images
                    if filetype != 'mask':
                        hdulist[1].data = hdulist[1].header['BOFFSET'] + hdulist[1].header['BSOFTEN'] * (10**(0.4*hdulist[1].data) - 10**(-0.4*hdulist[1].data))
                        del hdulist[1].header['BLANK']

                    # Getting rid of the edge nans
                    fixed_hdu = trim_nan(hdulist[1])
                else:
                    fixed_hdu = trim_nan(hdulist[0])

                del fixed_hdu.header['RADESYS']
            
            # Since the weight images are actually variance maps, we need
            # to take the reciprocal of it
            if filetype == 'wt':
                fixed_hdu.data = 1/fixed_hdu.data

            data_lists[filetype][flt].append(fixed_hdu)

        for skycell in data_lists['stack'][query.filters[0]]:
            #trimming down the grid, leaving the points currently not in any skycell
            grid = [p for p in grid if not p.contained_by(WCS(skycell.header))]
            
        if len(grid) == 0:
            break

    return data_lists['stack'], data_lists['wt']


def trim_nan(_hdu):
    '''
    When downloading an image from PS1, the image is padded with Nan.
    This function, trim the image, removing the np.nan array, and preserving
    the WCS.

    NOTE: the flux of the full skycell stack images are scaled. Astropy
    rescale the data automatically as soon as the data attribute of the
    hdu is accessed. See https://docs.astropy.org/en/stable/io/fits/usage/image.html
    for more details.

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :param _hdu: HUDList containing the data to trim
    :type _hdu: `~astropy.io.fits.HDUList`
    
    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :return: HDUList with edges trimmed
    :rtype: `~astropy.io.fits.HDUList`

    '''

    from astropy.nddata import Cutout2D
    import numpy as np
    import math

    # Getting all the columns with only np.nan in them
    good_cols = np.all(np.isnan(_hdu.data),axis=1)

    # Retrieving the first and the last index of good_cols == False,
    # i.e. the first and last columns containing (at least partially)
    # real data
    for i in range(len(good_cols)):
        if not good_cols[i]:
            # Adding 10 pixel padding to account for edge artifact
            first_col = i+10
            break
    
    for i in range(len(good_cols)-1,-1,-1):
        if not good_cols[i]:
            # We want to include the index, so we add 1, butting the
            # edge right after the column. Subtracting 10 pixel padding
            # to account for edge artifact.
            last_col = i+1-10
            break


    # Doing the same with the rows
    good_rows = np.all(np.isnan(_hdu.data),axis=0)
    for i in range(len(good_rows)):
        if not good_rows[i]:
            first_row = i+10
            break
            
    for i in range(len(good_rows)-1,-1,-1):
        if not good_rows[i]:
            last_row = i+1-10
            break

    size_x = last_col-first_col
    size_y = last_row-first_row
    cen_x = math.floor((last_col+first_col)/2.)
    cen_y = math.floor((last_row+first_row)/2.)
    
    # No idea why x and y are inverted between the center and the size
    # ¯\_(ツ)_/¯
    cut = Cutout2D(_hdu.data, (cen_y,cen_x), (size_x,size_y), wcs=WCS(_hdu.header))

    out_hdu = pf.PrimaryHDU()
    out_hdu.data = cut.data
    out_hdu.header = _hdu.header.copy()
    out_hdu.header.update(cut.wcs.to_header())

    return out_hdu

def track_skycells_in_header(_header, _data):
    # Keeping track of the single skycells that are being combined.
    _header['NSKYCELL'] = (len(_data), 'Number of SDSS skycells merged')
    for i,im in enumerate(_data, start=1):
        _header[f'STK_TY{i:>02}'] = (im.header['STK_TYPE'], f'type of stack {i}')
        _header[f'STK_ID{i:>02}'] = (im.header['STK_ID'], f'type of stack {i}')
        _header[f'SKYCEL{i:>02}'] = (im.header['SKYCELL'], f'type of stack {i}')
        _header[f'TES_ID{i:>02}'] = (im.header['TESS_ID'], f'type of stack {i}')
    
    return _header