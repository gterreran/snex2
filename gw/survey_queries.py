# -*- coding: utf-8 -*-
"""survey_queries.py

"""

from astropy.table import Table
from astropy import units as u
import astropy.io.fits as pf
from astropy.time import Time
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
import astropy.units as u
import os, math
import numpy as np
from astropy.wcs import WCS

# Folder definitions
OUT_FOLDER = os.getenv('GWTEMPLATES_FOLDER')
#OUT_FOLDER = '/supernova/data/extdata/GW_templates/O4'

# PS1 urls
PS1_TAB = 'https://ps1images.stsci.edu/cgi-bin/ps1filenames.py'

# PS1 provides a cutout API. We probably prefer to use the whole
# skycell, but we can easily change that with this bool variable
SKYCELL = True
if SKYCELL:
    PS1_CUT = 'https://ps1images.stsci.edu/'
else:
    PS1_CUT = 'https://ps1images.stsci.edu/cgi-bin/fitscut.cgi'

# Specifications for different instruments
# The FOV here is thought as the diameter of a circle that would circumscribe the LCO image.
# This is effectively bigger than the actual FOV of each instrument, but it will guarantee
# that the full FOV is covered given any rotation 
INSTRUMENT_LIST = ['sinistro', 'muscat', 'qhy']
#FOV = {'sinistro': 26.5*np.sqrt(2)*u.arcmin, 'muscat': 9.1*np.sqrt(2)*u.arcmin, 'qhy': np.sqrt(1.9**2 + 1.2**2)*u.deg}
# The FOV of 'qhy' is rectangular, se we calculate the side of square having the same diagonal.
FOV = {'sinistro': 26.5*u.arcmin, 'muscat': 9.1*u.arcmin, 'qhy': np.sqrt((1.9**2 + 1.2**2)/2)*u.deg}
RESOLUTION = {'sinistro':0.389*u.arcsec, 'muscat': 0.27*u.arcsec, 'qhy':0.74*u.arcsec}

SURVEYS_SKYCELL_SIZES = {'PS1': 0.4*u.deg, 'SDSS': 10*u.arcmin}



class survey_request:
    '''
    This class will query the surveys looking for the optical templates.
    Each function is associated with a survey. First it will check if
    the coordinates are in the footprint of that survey. Each function
    will return True or False accordingly. In addition, if the coordinates
    are indeed in the footprint, then the function will store the hdu of
    the template in the 'hdu' attribute.

    '''

    def __init__(self, _obj, _coord, _filters, _inst):
        '''
        Initializing the class with the name of the object :_obj:,
        the coordinates :_coord: and the filters :_filters: to search.
        The template images will be resampled to match the resolution
        of the instrument (or list of instruments) provided by :_inst:.
        The template FOV will also be as big as to include the whole
        FOV of that instrument.
        The empty dictionary :templates_paths: is also created.

        *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        :_obj: Name of the object, mainly for file naming purposes
        :type _obj: string

        :_coord: coordinates of the object, corresponding to the center
                 of the search
        :type _coord: astropy.coordinates.sky_coordinate.SkyCoord

        :_filters: filters to search. If searching for multiple filters
                   use a single long string, e.g. 'gri'
        :type _filters: string

        :_inst: instrument name, or list of instrument names. In order 
                to be valid, the name will be checked against the global
                variable :INSTRUMENT_LIST:, defined above.
        :type _inst: string or list

        '''
        self.obj = _obj
        self.coord = _coord
        self.filters = _filters

        # Making the :instruments: list and checking the instruments exist.
        if isinstance(_inst, str):
            if _inst == 'all':
                self.instruments = INSTRUMENT_LIST
            else:
                self.instruments = [_inst]
        else:
            self.instruments = _inst
        for i in self.instruments:
            if i not in INSTRUMENT_LIST:
                raise ValueError(f'Instrument `{i}` does not appear to be currently supported. If this is a new instrument, you will need to included it in the INSTRUMENT_LIST, as well as update the FOV and RESOLUTION dictionaries.')
        
        # Sorting the instruments in decreasing FOV size.
        self.instruments.sort(key = lambda x: FOV[x], reverse=True)

        self.templates_paths = {}

    
    # Defining a decorator. Because decorators are cool
    def set_survey(_survey):
        '''
        This decorator set the name of survey for the following function, through the
        :_survey: variable, and initialize the :self.templates_paths: dictionary with
        a '--' string for that specific survey. The idea is that by initializing the
        dictionary, there will be a record that the search has been performed.
        Otherwise the dictionary will not have any entry for that survey.

        *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        :_survey: Name of the survey
        :type _center: string

        *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        :return: Decorated function
        :rtype: function

        '''
        def decorate(function):
            def initialize_survey(self):
                # Initializing all the filters
                for flt in self.filters:
                    for inst in self.instruments:
                        self.templates_paths[f'{inst}_{_survey}_{flt}'] = '--'
                print(f'Searching for {_survey} templates.\n')
                return function(self, _survey)
            return initialize_survey
        return decorate


#######   ---PS1---   #######
    @set_survey('PS1')
    def search_for_PS1(self, survey):
        '''
        This function searches and downloads PS1 templates.
        Together with the data, it will fetch also the weight images (i.e. 1/variance)
        and the mask images.

        '''

        # Instruments are sorted by FOV size.
        # :self.instruments[0]: is the instrument with the biggest FOV
        # and therefore requires the biggest template. We start creating
        # this, and instead of repeating the procedure, we just scale down
        # this first template to match the smaller requirements. 
        inst = self.instruments[0]

        # Generate grid of points covering the FOV
        grid = generate_FOV_grid(_center = self.coord, _fov = FOV[inst], _step=SURVEYS_SKYCELL_SIZES[survey])

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
            tableurl = f'{PS1_TAB}?ra={grid[0].ra.deg}&dec={grid[0].dec.deg}&size={FOV[inst].to(u.arcsec).value/0.25}&format=fits&filters={self.filters}&type=stack,stack.wt,stack.mask'
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

            for skycell in data_lists['stack'][self.filters[0]]:
                #trimming down the grid, leaving the points currently not in any skycell
                grid = [p for p in grid if not p.contained_by(WCS(skycell.header))]
                
            if len(grid) == 0:
                break
        
        # We now trim the list of skycells according to each instrument FOV
        # For the first instrument this step should be redundant since we
        # selected the skycell in order to tile the FOV already.
        for inst in self.instruments:
            new_grid = generate_FOV_grid(_center = self.coord, _fov = FOV[inst], _step=SURVEYS_SKYCELL_SIZES[survey])
            for flt in data_lists['stack']:
                # Since enumerate() returns a generator and generators can't be reversed,
                # we need to convert it to a list first
                for skycell_idx, skycell in reversed(list(enumerate(data_lists['stack'][flt]))):
                    contained = False
                    for p in new_grid:
                        if p.contained_by(WCS(skycell.header)):
                            contained = True
                            break
                    # If no point of the grid is contained by the skycell,
                    # drop it from all the data lists
                    if not contained:
                        for filetype in data_lists:
                            data_lists[filetype][flt].pop(skycell_idx)
                
                print(f'Reprojecting and combining skycells of filer {flt} for {inst}.')
                array, footprint, header = reproject_and_combine(data_lists['stack'][flt], data_lists['wt'][flt], inst)

                # Keeping track of the single skycells that are being combined.
                header['NSKYCELL'] = (len(data_lists['stack'][flt]), 'Number of SDSS skycells merged')
                for i,im in enumerate(data_lists['stack'][flt], start=1):
                    header[f'STK_TY{i:>02}'] = (im.header['STK_TYPE'], f'type of stack {i}')
                    header[f'STK_ID{i:>02}'] = (im.header['STK_ID'], f'type of stack {i}')
                    header[f'SKYCEL{i:>02}'] = (im.header['SKYCELL'], f'type of stack {i}')
                    header[f'TES_ID{i:>02}'] = (im.header['TESS_ID'], f'type of stack {i}')

                # Populating the the header with the keys needed by the pipeline
                header['RA'] = self.coord.ra.deg
                header['DEC'] = self.coord.dec.deg
                header['CAT-RA'] = self.coord.ra.deg
                header['CAT-DEC'] = self.coord.dec.deg
                header['MJD-OBS'] = (np.mean([im.header['MJD-OBS'] for im in data_lists['stack'][flt]]), 'Average MJD')
                header['DATE-OBS'] = 'T'.join(Time(header['MJD-OBS'], format='mjd', scale='utc').iso.split())
                header['DAY-OBS'] = header['DATE-OBS'].split('T')[0].replace('-', '')
                header['TELESCOP'] = 'PS1'
                header['INSTRUME'] = 'GPC1'
                header['OBJECT'] = self.obj
                header['FILTER'] = flt+'p' # NOTE: this works only with griz
                header['AIRMASS'] = 1
                header['WCSERR'] = 0
                header['PIXSCALE'] = (RESOLUTION[inst].value, 'arcsec/pixel')
                header['BUNIT'] = ('counts', 'ADUs')

                # Writing out the final combined template
                out_name = f'{OUT_FOLDER}/{survey}/{inst}_{self.obj}_{flt}'
                self.templates_paths[f'{inst}_{survey}_{flt}'] = out_name+'.fits'
                pf.writeto(out_name+'.fits', array, header, overwrite=True)

                # Hotpants needs the variance, so we need to write out 1/footprint
                header['BUNIT'] = ('counts2', 'ADUs^2')
                pf.writeto(out_name+'.var.fits', 1/footprint, header, overwrite=True)
        



#######   ---SDSS---   #######
    @set_survey('SDSS')
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
        inst = self.instruments[0]

        # Generate grid of points covering the FOV
        grid = generate_FOV_grid(_center = self.coord, _fov = FOV[inst], _step=SURVEYS_SKYCELL_SIZES[survey])

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
                    print(f'Querying the SDSS database for coordinates {grid[0].ra.deg} {grid[0].dec.deg}.')
                    xid = SDSS.query_region(grid[0], radius=3*u.arcmin, spectro=False)
                    if xid is None:
                        print(f'\nThe provided coordinates {grid[0].ra.deg},{grid[0].dec.deg} do not appear to be in the SDSS footprint.\n')
                        return
                    
                    print('Coordinates in SDSS footprint.')
                    new_pointings = [list(el) for el in np.unique(xid['run','camcol','field']) if list(el) not in all_pointings]

                    all_pointings = all_pointings + new_pointings

                    for run, camcol, field in new_pointings:
                        #fetching the fits as an astropy hdu
                        print(f'Retrieving file with run={run}, camcol={camcol}, field={field} and filter {flt}.')
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

            
        # We now trim the list of skycells according to each instrument FOV
        # For the first instrument this step should be redundant since we
        # selected the skycell in order to tile the FOV already.
        for inst in self.instruments:
            new_grid = generate_FOV_grid(_center = self.coord, _fov = FOV[inst], _step=SURVEYS_SKYCELL_SIZES[survey])
            for flt in img_data_list:
                # Since enumerate() returns a generator and generators can't be reversed,
                # we need to convert it to a list first
                for skycell_idx, skycell in reversed(list(enumerate(img_data_list[flt]))):
                    contained = False
                    for p in new_grid:
                        if p.contained_by(WCS(skycell.header)):
                            contained = True
                            break
                    # If no point of the grid is contained by the skycell,
                    # drop it from all the data lists
                    if not contained:
                        img_data_list[flt].pop(skycell_idx)
                        weight_data_list[flt].pop(skycell_idx)
                
                print(f'Reprojecting and combining skycells of filer {flt} for {inst}.')
                array, footprint, header = reproject_and_combine(img_data_list[flt], weight_data_list[flt], inst)



                header['NTILES'] = (len(img_data_list[flt]), 'Number of SDSS tiles merged')
                for i,im in enumerate(img_data_list[flt], start=1):
                    header[f'RUN{i:>03}'] = (im.header['RUN'], f'Run number of tile {i}')
                    header[f'CAMCO{i:>03}'] = (im.header['CAMCOL'], f'Column in the imaging camera of tile {i}')
                    header[f'FIELD{i:>03}'] = (im.header['FIELD'], f'Field number of tile {i}')
                    header[f'MJD{i:>03}'] = (im.header['MJD-OBS'], f'MJD of tile {i}')

                header['RA'] = self.coord.ra.deg
                header['DEC'] = self.coord.dec.deg
                header['CAT-RA'] = self.coord.ra.deg
                header['CAT-DEC'] = self.coord.dec.deg
                header['MJD-OBS'] = (np.mean([im.header['MJD-OBS'] for im in img_data_list[flt]]), 'Average MJD')
                header['DATE-OBS'] = 'T'.join(Time(header['MJD-OBS'], format='mjd', scale='utc').iso.split())
                header['DAY-OBS'] = header['DATE-OBS'].split('T')[0].split('T')[0].replace('-', '')
                header['TELESCOP'] = '2.5m'
                header['INSTRUME'] = 'SDSS'
                header['OBJECT'] = self.obj
                header['FILTER'] = flt+'p' #note, this works only with griz
                header['AIRMASS'] = 1
                header['WCSERR'] = 0
                header['PIXSCALE'] = (RESOLUTION[inst].value, 'arcsec/pixel')
                header['BUNIT'] = ('counts', 'ADUs')

                out_name = f'{OUT_FOLDER}/{survey}/{inst}_{self.obj}_{flt}'
                self.templates_paths[f'{inst}_{survey}_{flt}'] = out_name+'.fits'
                pf.writeto(out_name+'.fits', array, header, overwrite=True)

                #hotpants needs the variance, so we need the write out 1/footprint
                header['BUNIT'] = ('counts2', 'ADUs^2')
                pf.writeto(out_name+'.var.fits', 1/footprint, header, overwrite=True)

    

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


def generate_FOV_grid(_center, _fov, _step, circle=True, wcs = None):
    '''
    The function generates a grid of points, each distant :_step: to each other.
    The construction follows a recursive hexagonal ring, that becomes bigger and
    bigger until it encompass the whole FOV. It takes advantage of the class 
    `~astropy.visualization.wcsaxes.SphericalCircle`, which creates a circle that
    is formed of all the points that are within a certain angle of the central
    coordinates on a sphere. Each hexagonal ring is built on top of these circles.

    If :_step: is bigger than half of the fov size, then the grid would only have
    the central point. This can cause the FOV to not be fully covered by the 
    template tiles. To fix this, in this case the grid will be just the central
    point plus an octagon with a diagonal equal to half of the diagonal of the
    FOV. The octagon is chosen to account for any rotation of the FOV.

    If :circle: is True (default), the grid will trimmed to a circle with radius
    equal to :_fov:
    
    If :wcs: is provided, the points outside of the provided WCS are removed.

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :param _center: Coordinates of the center point of the grid
    :type _center: `~astropy.units.Quantity` ['angle']

    :param _fov: size of the FOV to tile with the grid. This assumes a square
                 FOV. This quantity effectively translates to the length of the
                 apothem of outermost hexagon, which will be equal to half of
                 the diagonal of the square FOV.
    :type _fov: `~astropy.units.Quantity` ['angle']

    :param _step: distance of each point in the grid. If a value is not
                  provided, it will default to half :_fov:
    :type _step: `~astropy.units.Quantity` ['angle']

    :param circle: If True, trims the grid to a circle with radius :_fov:
    :type circle: bool
    :default: True

    :param wcs: WCS information. If provided the points outside of it will
                be removed
    :type wcs: `~astropy.wcs.WCS` or None
    :default wcs: None

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :return: List of points of the grid. Each element is a astropy SkyCoord
    :rtype: list

    '''

    from astropy.visualization.wcsaxes import SphericalCircle

    points = [[_center, 0*u.deg]]

    # If the the step is too big, we will take an octagonal ring
    if _step > _fov/2:
        vertexes = SphericalCircle(_center, radius = _fov*np.sqrt(2)/2, resolution=8).get_verts()
        # :SphericalCircle: will produce a Polygon, so the first and the last
        # point will coincide. Therefore we need to remove the last vertex.
        for v in vertexes[:-1]:
            p = SkyCoord(*v, unit=u.deg)
            points.append([p, p.separation(_center)])

    else:
        # A square FOV of side A, its diagonal is A*sqrt(2). Therefore an
        # hexagon will have to have a diagonal of A*sqrt(2)/sqrt(3) in order to contain
        # it all, given any rotation. We add an extra order just to be sure

        for ring in range(1, math.ceil(_fov/_step*np.sqrt(2)/np.sqrt(3))+1):
            current_hexagon = []
            vertexes = SphericalCircle(_center, radius = _step*ring, resolution=6).get_verts()
            for v in vertexes:
                p = SkyCoord(*v, unit=(u.deg,u.deg))
                current_hexagon.append([p, p.separation(_center)])
            #Filling the sides of the hexagon with points
            for side in range(6):
                # We need to split the side on the nth hexagon
                # in n segments
                for i in range(1, ring):
                    a = i / ring            # rescale 0 < i < n --> 0 < a < 1
                    new_point = (1 - a) * current_hexagon[side][0].cartesian + a * current_hexagon[side+1][0].cartesian # interpolate cartesian coordinates
                    new_point = SkyCoord(new_point)
                    current_hexagon.append([new_point, new_point.separation(_center)])
            # :SphericalCircle: will produce a Polygon, so the first and the last
            # point will coincide. Therefore we need to remove the last vertex of
            # the current hexagon.
            current_hexagon.pop(6)
            points = points + current_hexagon

    if circle:
        #adding 5 arcsec to give an extra buffer
        points = [p for p in points if p[1] <= (_fov/2*np.sqrt(2))+5*u.arcsec]

    if wcs:
        points = [p for p in points if p[0].contained_by(wcs)]
    
    # Sorting the list according to the distance, with the most distant points first
    points.sort(key = lambda x: x[1], reverse=True)

    #Putting the center back at the beginning
    points = [points.pop()] + points

    #Returning just the list of coordinates, without the distance
    return [p[0] for p in points]


def reproject_and_combine(_imglist, _weightlist, _instrument):
    '''
    The name of the function is quite self explanatory. 
    It uses the `reproject` package to match the resolution 
    of :_instrument:, and to merge all the tiles
    into the final template. It produces also a weight
    image.
    
    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :param _imglist: List of hdus with the images to merge
    :type _imglist: list

    :param _weightlist: List of hdus with the weights of the images
                        to merge
    :type _weightlist: list

    :param _instrument: Name of the instrument. The string needs to
                        be a valid key in :RESOLUTION:
    :type _instrument: string


    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :return: Image data, Weight data, Header with optimal WCS
    :rtype: list
    '''
    from reproject.mosaicking import find_optimal_celestial_wcs
    from reproject import reproject_interp
    from reproject.mosaicking import reproject_and_coadd

    wcs_out, shape_out = find_optimal_celestial_wcs(_imglist, resolution=RESOLUTION[_instrument])
    array, footprint = reproject_and_coadd(_imglist, wcs_out, shape_out=shape_out, reproject_function=reproject_interp, input_weights = _weightlist, combine_function='mean')
    return array, footprint, wcs_out.to_header()


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
