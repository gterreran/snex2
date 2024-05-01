# -*- coding: utf-8 -*-
"""survey_queries.py

"""

from astropy.table import Table
from astropy.io.votable import parse
from astropy import units as u
import astropy.io.fits as pf
from astropy.time import Time
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
import astropy.units as u
import urllib.request, io
from numpy.core.defchararray import startswith
from pyvo.dal import sia
import os
import numpy as np
from astropy.wcs import WCS

#Folder definitions
OUT_FOLDER = os.getenv('GWTEMPLATES_FOLDER')
#OUT_FOLDER = '/supernova/data/extdata/GW_templates/O4'

#PS1 urls
PS1_TAB = 'https://ps1images.stsci.edu/cgi-bin/ps1filenames.py'
PS1_CUT = 'https://ps1images.stsci.edu/cgi-bin/fitscut.cgi'

#List of available services at:
#https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/DataAccessInterfaces/SimpleImageAccessSIA/SimpleImageAccessSIA.html
#COLLECTION = 'coadd_all'
COLLECTION = 'ls_dr9'
COLLECTION = 'coadd/decaps_dr2'
DECAM_SERVICE = f'https://datalab.noirlab.edu/sia/{COLLECTION}'

#Skymapper urls
SKY_CUT = 'https://api.skymapper.nci.org.au/public/siap/dr4/query'

FOV = 26.5*u.arcmin #Sinistro FOV is 26.5'x26.5'
RESOLUTION = 0.389*u.arcsec 

class survey_request:
    '''
    This class will query the surveys looking for the optical templates.
    Each function is associated with a survey. First it will check if
    the coordinates are in the footprint of that survey. Each function
    will return True or False accordingly. In addition, if the coordinates
    are indeed in the footprint, then the function will store the hdu of
    the template in the 'hdu' attribute.

    '''
    def __init__(self, _obj, _coord, _filters):
        self.obj = _obj
        self.coord = _coord
        self.filters = _filters
        self.templates_paths = {}

    
    def set_survey(_survey):
        def decorate(function):
            def initialize_survey(self):
                for flt in self.filters:
                    self.templates_paths[f'{_survey}_{flt}'] = '--'
                return function(self, _survey)
            return initialize_survey
        return decorate


#######   ---PS1---   #######
    @set_survey('PS1')
    def search_for_PS1(self, survey):
        '''
        PS1 needs the size in pixels, considering 0.25 arcsec/pixel.

        '''

        tableurl = f'{PS1_TAB}?ra={self.coord.ra.deg}&dec={self.coord.dec.deg}&size={FOV.to(u.arcsec).value/0.25}&format=fits&filters={self.filters}&type=stack,stack.wt,stack.mask'
        table = Table.read(tableurl, format='ascii')
        for tab in table:
            flt = tab['filter']
            filename = tab['filename']
            url = f'{PS1_CUT}?ra={self.coord.ra.deg}&dec={self.coord.dec.deg}&size=2500&format=fits&red={filename}'
            #out_name = f'{OUT_FOLDER}/{survey}/{self.obj}_{flt}.fits'
            out_name = f'{OUT_FOLDER}/{survey}/{self.obj}_{flt}.{".".join(filename.split(".")[10:])}'
            urllib.request.urlretrieve(url, out_name)
            self.templates_paths[f'{survey}_{flt}'] = out_name
            with pf.open(out_name, mode='update') as hdu:
                date_obs = Time(hdu[0].header['MJD-OBS'], format='mjd', scale='utc').iso.split()
                hdu[0].header['DATE-OBS'] = 'T'.join(date_obs)
                hdu[0].header['DAY-OBS'] = date_obs[0]
                hdu[0].header['FILTER'] = flt+'p' #note, this works only with griz
                hdu[0].header['TELESCOP'] = 'PS1'
                hdu[0].header['INSTRUME'] = 'GPC1'
                hdu[0].header['OBJECT'] = self.obj
                hdu[0].header['WCSERR'] = 0
                hdu[0].header['RA'] = self.coord.ra.deg
                hdu[0].header['DEC'] = self.coord.dec.deg

                hdu.flush()



#######   ---SDSS---   #######
    @set_survey('SDSS')
    def search_for_SDSS(self, survey):
        '''
        Each tile has a 10 by 13 arcminutes FOV.

        '''
        
        grid = generate_FOV_grid(self.coord, FOV)

        all_pointings = []

        for f,flt in enumerate(self.filters):

            img_data_list = []
            weight_data_list = []

            # For the first filter, we need to collect the references for all the required SDSS fields
            # We don't need to repeat the query for the following filters
            # The SDSS fields references are stored in `all_pointings`
            if f == 0:
                while True:
                    # SDSS queries radius must be less than 3.0 arcmin.
                    # initially `grid[0]` will be the center of the field.
                    # After each iteration, the `grid` array is trimmed down and its first element will be
                    # the point on the grid most distant from the center
                    xid = SDSS.query_region(grid[0], radius=3*u.arcmin, spectro=False)
                    if xid is None:
                        return
                    
                    new_pointings = [list(el) for el in np.unique(xid['run','camcol','field']) if list(el) not in all_pointings]

                    all_pointings = all_pointings + new_pointings

                    for run, camcol, field in new_pointings:
                        #fetching the fits as an astropy hdu
                        im = SDSS.get_images(run=run, camcol=camcol, field=field, band=flt, cache=True)[0]

                        sky_subtracted_original_image_hdu, weight_image_hdu = elaborate_SDSS(im, flt, run, camcol, field)
                        img_data_list.append(sky_subtracted_original_image_hdu)
                        weight_data_list.append(weight_image_hdu)
                        
                
                    for tile in img_data_list:
                        grid = [p for p in grid if not p.contained_by(WCS(tile.header))]
                    
                    if len(grid) == 0:
                        break
                
            else:
                #Note that we changed the list to all_pointings now that we have them all
                for run, camcol, field in all_pointings:
                    #fetching the fits as an astropy hdu
                    im = SDSS.get_images(run=run, camcol=camcol, field=field, band=flt, cache=True)[0]

                    sky_subtracted_original_image_hdu, weight_image_hdu = elaborate_SDSS(im, flt, run, camcol, field)
                    img_data_list.append(sky_subtracted_original_image_hdu)
                    weight_data_list.append(weight_image_hdu)

            
            array, footprint, header = reproject_and_combine(img_data_list, weight_data_list)

            header['NTILES'] = (len(img_data_list), 'Number of SDSS tiles merged')
            for i,im in enumerate(img_data_list, start=1):
                header[f'RUN{i:>03}'] = (im.header['RUN'], f'Run number of tile {i}')
                header[f'CAMCO{i:>03}'] = (im.header['CAMCOL'], f'Column in the imaging camera of tile {i}')
                header[f'FIELD{i:>03}'] = (im.header['FIELD'], f'Field number of tile {i}')
                header[f'MJD{i:>03}'] = (im.header['MJD-OBS'], f'MJD of tile {i}')


            header['MJD-OBS'] = (np.mean([im.header['MJD-OBS'] for im in img_data_list]), 'Average MJD')
            header['DATE-OBS'] = 'T'.join(Time(header['MJD-OBS'], format='mjd', scale='utc').iso.split())
            header['DAY-OBS'] = header['DATE-OBS'].split('T')[0]
            header['TELESCOP'] = '2.5m'
            header['INSTRUME'] = 'SDSS'
            header['OBJECT'] = self.obj
            header['FILTER'] = flt+'p' #note, this works only with griz
            header['AIRMASS'] = 1
            header['WCSERR'] = 0
            header['BUNIT'] = ('counts', 'ADUs')

            out_name = f'{OUT_FOLDER}/{survey}/{self.obj}_{flt}'
            self.templates_paths[f'{survey}_{flt}'] = out_name+'.fits'
            pf.writeto(out_name+'.fits', array, header, overwrite=True)

            #hotpants needs the variance, so we need the writeout 1/footprint
            header['BUNIT'] = ('counts2', 'ADUs^2')
            pf.writeto(out_name+'.var.fits', 1/footprint, header, overwrite=True)




#######   ---Skymapper---   #######
    @set_survey('Skymapper')        
    def search_for_Skymapper(self, survey):
        '''
        Skymapper has the ability to search for multiple filters in one go,
        but the returned list of file loses the information on the filter,
        so I need to send a seperate query for each filter.
        Biggest size for Skymapper is 0.17 deg

        '''
        from astropy import units as u
        for flt in self.filters:

            url = f'{SKY_CUT}?POS={self.coord.ra.deg},{self.coord.dec.deg}&SIZE=0.17,0.17&BAND={flt}&FORMAT=image/fits&INTERSECT=covers&RESPONSEFORMAT=VOTABLE'


            u = urllib.request.urlopen(url)
            s = io.BytesIO(u.read())
            votable = parse(s)

            for resource in votable.resources:
                for t in resource.tables:
                    table = t.array['get_fits'][:]

            if len(table) != 0:
                out_name = f'{OUT_FOLDER}/{survey}/{self.obj}_{flt}.fits'
                urllib.request.urlretrieve(table[0], out_name)
                self.templates_paths[f'{survey}_{flt}'] = out_name
                with pf.open(out_name, mode='update') as hdu:
                    hdu[0].header['DAY-OBS'] = hdu[0].header['DATE-OBS'].split('T')[0]
                    hdu[0].header['FILTER'] = flt+'p' #note, this works only with griz
                    hdu[0].header['OBJECT'] = self.obj
                    hdu[0].header['WCSERR'] = 0
                    hdu[0].header['RA'] = self.coord.ra.deg
                    hdu[0].header['DEC'] = self.coord.dec.deg

                    hdu.flush()
            
            else:
                print(f'Coordinates {self.coord.ra.deg} {self.coord.dec.deg} not in Skymapper with filter {flt}.')
                
            

#######   ---DECam---   #######
    @set_survey('DECam')
    def search_for_DECam(self, survey):

        connect = sia.SIAService(DECAM_SERVICE)
        table = connect.search(pos = (self.coord.ra.deg, self.coord.dec.deg), size = (FOV.to(u.deg).value, FOV.to(u.deg).value), verbosity=2).to_table()

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
    

def get_SDSS_gain_and_darkVariance(_flt, _run, _camcol):
    #retrieve gain and darkVariance from another online table
    q = SDSS.query_sql(f'SELECT DISTINCT gain_{_flt}, darkVariance_{_flt} FROM field WHERE run={_run} AND camcol={_camcol}')
    gain = q[f'gain_{_flt}'][0]
    darkVariance = q[f'darkVariance_{_flt}'][0]
    return gain, darkVariance

def elaborate_SDSS(_im, _flt, _run, _camcol, _field):
    header = _im[0].header

    #Although we do not write out each single tile, this edits in the header will
    #make it easier to keep track of each tile in the merged fits file
    header['RADESYSa'] = header['RADECSYS']
    del header['RADECSYS']
    header['MJD-OBS'] = Time(f"{header['DATE-OBS']} {header['TAIHMS']}", format='iso', scale='tai').mjd
    header['DATE-OBS'] = 'T'.join(Time(header['MJD-OBS'], format='mjd', scale='utc').iso.split())
    #For some reason the field info is not stored in the header. There is a "FRAME" key, but it's different
    header['FIELD'] = _field

    gain, darkVariance = get_SDSS_gain_and_darkVariance(_flt, _run, _camcol)

    #im[2] contains the sky
    sky_image = SDSS_sky(_im[2]) #in ADU

    #im[1] contains a row with factors nanomaggies per counts for each column
    calib_image = np.asarray([_im[1].data.tolist()]*np.shape(_im[0].data)[0])

    sky_subtracted_original_image = _im[0].data / calib_image #in ADU

    variance_image = ((_im[0].data / calib_image + sky_image) / gain  + darkVariance)# in ADU * calib_image **2

    sky_subtracted_original_image_hdu = pf.PrimaryHDU(sky_subtracted_original_image, header=header)
    #to reproject and combine, we just need the weight, not the variance
    weight_image_hdu = pf.PrimaryHDU(1/variance_image, header=header)

    return sky_subtracted_original_image_hdu, weight_image_hdu

    


def SDSS_sky(_hdu, method="linear"):
    '''
    The function converts the second entry of a hdu SDSS frame (sky) to a good
    sky array.

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :param _hdu: HDUList from SDSS database
    :type _hdu: astropy.io.fits.hdu.hdulist.HDUList

    :param method: Method of interpolation. In can be "linear" (default)
    :type method: str 
    :default: "linear"

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :return: Sky array, in ADUs.
    :rtype: numpy.ndarray

    '''

    from scipy import interpolate

    sky, x, y = _hdu.data[0]
    xold = np.arange(np.shape(sky)[0])
    yold = np.arange(np.shape(sky)[1])

    # `interp2d` is deprecated in SciPy 1.10 and will be removed in SciPy 1.14.0
    # NOTE: `interp2d` has x and y inverted but doesn't need a `meshgrid` to produce the output
    # This is the code if `interp2d` is still preferred despite of its deprecation
    # tck  = interpolate.interp2d(yold, xold, sky, kind=method)
    # return tck(x, y)

    tck = interpolate.RegularGridInterpolator((xold,yold), sky, method=method, bounds_error=False, fill_value = None)

    # need to invert x and in y to produce the meshgrid with the right shape
    x_grid, y_grid = np.meshgrid(y, x, indexing='ij')

    return tck((x_grid, y_grid))


def generate_FOV_grid(_center, _fov, _step = None, circle=True, wcs = None):
    '''
    The function generates a grid of points, each distant :_step: to each other.
    The construction follows a recursive hexagon ring, that becomes bigger and
    bigger until it encompass the whole FOV. It takes advantage of equilateral
    triangles properties, where the height to the top vertex is sqrt(3)/2 times
    the side of the triangle. 

    If :circle: is True (default), the grid will trimmed to a circle with radius
    equal to :_fov:
    
    If :wcs: is provided, the points outside of the provided WCS are removed.

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :param _center: Coordinates of the center point of the grid
    :type _center: astropy.coordinates.sky_coordinate.SkyCoord

    :param _fov: size of the FOV to tile with the grid. This assumes a square
                 FOV. This quantity effectively translates to the length of the
                 apothem of outermost hexagon, which will be equal to half of
                 the diagonal of the square FOV.
    :type _fov: astropy.units.quantity.Quantity

    :param _step: distance of each point in the grid. If a value is not
                  provided, it will default to half :_fov:
    :type _fov: astropy.units.quantity.Quantity

    :param circle: If True, trims the grid to a circle with radius :_fov:
    :type: bool
    :default: True

    :param wcs: WCS information. If provided the points outside of it will
                be removed
    :type wcs: astropy.wcs.wcs.WCS or None
    :default: None

    *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    :return: List of points of the grid. Each element is a astropy SkyCoord
    :rtype: list

    '''

    import math, itertools
    

    if _step == None:
        _step = _fov/3

    points = [[_center, 0*u.deg]]

    # Assuming a square FOV of side A, its diagonal is A*sqrt(2). Therefore an
    # hexagon will have to have a side of A*sqrt(2)/sqrt(3) in order to contain
    # it all, given any rotation.

    for ring in range(1, math.ceil(_fov/_step*np.sqrt(2)/np.sqrt(3))):
        # Moving `ring*_step` E from center to start the new hexagon
        p = SkyCoord(_center.ra+_step*ring/np.cos(_center.dec.rad), _center.dec)
        points.append([p, p.separation(_center)])

        # The following `for` loops are written in a fancy way. Since we don't need to store
        # the iterator variable, and we just need to repeat the same operation N times, we
        # can use `itertools.repeat`, which doesn't need to build and return a list, storing
        # the index as well, like range, or xrange would do. This is effectively the fastest
        # way to iterate N times in Python.
        # A more readable and common way of writing this would be
        # `for _ in range(N):`
        # And can be written like this if preferred, with negligible loss in performance. 
        # The number N of points along each hexagon side is exactly the Nth ring, so in our
        # case N=`ring`

        # The most distant points from the center are the most important, so I store the
        # distance to the center as well, to be able to sorted them later.

        #`p` it's always the previous point to move from

        # Turning 60 degrees and starting moving NW
        for _ in itertools.repeat(None, ring):
            p = SkyCoord(p.ra-_step/2/np.cos(p.dec.rad), p.dec+_step*np.sqrt(3)/2)
            points.append([p, p.separation(_center)])
        # Turning 120 degrees and starting moving W
        for _ in itertools.repeat(None, ring):
            p = SkyCoord(p.ra-_step/np.cos(p.dec.rad), p.dec)
            points.append([p, p.separation(_center)])
        # Turning 120 degrees and starting moving SW
        for _ in itertools.repeat(None, ring):
            p = SkyCoord(p.ra-_step/2/np.cos(p.dec.rad), p.dec-_step*np.sqrt(3)/2)
            points.append([p, p.separation(_center)])
        # Turning 120 degrees and starting moving SE
        for _ in itertools.repeat(None, ring):
            p = SkyCoord(p.ra+_step/2/np.cos(p.dec.rad), p.dec-_step*np.sqrt(3)/2)
            points.append([p, p.separation(_center)])
        # Turning 120 degrees and starting moving E
        for _ in itertools.repeat(None, ring):
            p = SkyCoord(p.ra+_step/np.cos(p.dec.rad), p.dec)
            points.append([p, p.separation(_center)])
        # Turning 120 degrees and starting moving NE
        for _ in itertools.repeat(None, ring):
            p = SkyCoord(p.ra+_step/2/np.cos(p.dec.rad), p.dec+_step*np.sqrt(3)/2)
            points.append([p, p.separation(_center)])

    if circle:
        points = [p for p in points if p[1] <= _fov/2*np.sqrt(2)]

    if wcs:
        points = [p for p in points if p[0].contained_by(wcs)]
    
    # Sorting the list according to the distance, with the most distant points first
    points.sort(key = lambda x: x[1], reverse=True)

    #Putting the center back at the beginning
    points = [points.pop()] + points

    #Returning just the list of coordinates, without the distance
    return [p[0] for p in points]


def reproject_and_combine(_imglist, _weightlist):#, _coord):
    from reproject.mosaicking import find_optimal_celestial_wcs
    from reproject import reproject_interp
    from reproject.mosaicking import reproject_and_coadd

    wcs_out, shape_out = find_optimal_celestial_wcs(_imglist, resolution=RESOLUTION)
    array, footprint = reproject_and_coadd(_imglist, wcs_out, shape_out=shape_out, reproject_function=reproject_interp, input_weights = _weightlist, combine_function='mean')
    return array, footprint, wcs_out.to_header()


    # import numpy as np
    # import matplotlib.pyplot as plt

    # plt.figure(figsize=(10, 8))
    # ax1 = plt.subplot(1, 2, 1, projection=wcs_out)
    # im1 = ax1.imshow(array, origin='lower', vmin=600, vmax=800)
    # ax1.set_title('Mosaic')
    # ax2 = plt.subplot(1, 2, 2, projection=wcs_out)
    # im2 = ax2.imshow(footprint, origin='lower')
    # ax2.set_title('Footprint')
    # ax1.plot([49.13674189,49.13062658,48.53523662,48.54543162]*u.deg, [ 41.39245501, 41.83609281,41.82997431, 41.38637844]*u.deg, color = 'r', lw=3, transform=ax1.get_transform("world"))
    # ax2.plot([49.13674189,49.13062658,48.53523662,48.54543162]*u.deg, [ 41.39245501, 41.83609281,41.82997431, 41.38637844]*u.deg, color = 'r', lw=3, transform=ax1.get_transform("world"))
    

    #grid = generate_FOV_grid(_coord, FOV)

    # for p in grid:
    #     ax1.plot(p.ra, p.dec, marker='o', color='r', ms=10, transform=ax1.get_transform("world"))
    #     ax2.plot(p.ra, p.dec, marker='o', color='r', ms=10, transform=ax2.get_transform("world"))

    # plt.show()

    # pf.writeto('test_combine.fits', array, wcs_out.to_header(), overwrite=True)
    # pf.writeto('test_footprint.fits', footprint, wcs_out.to_header(), overwrite=True)


# lco_fits = pf.open('/Users/giacomoterreran/lco/data/snexdata_target7856/tfn1m014-fa20-20231005-0153-e91.fits')
# wcs = WCS(lco_fits[0].header)
# footprint = wcs.calc_footprint()
# target_coord = SkyCoord(lco_fits[0].header['RA'], lco_fits[0].header['DEC'], unit=(u.hourangle, u.deg))

# import numpy as np
# import matplotlib.pyplot as plt

# im = pf.open('/Users/giacomoterreran/lco/data/data/extdata/O4/PS1/SN2023rbk_g.fits')

# plt.figure(figsize=(10, 8))
# ax1 = plt.subplot(1, 1, 1, projection=WCS(im[0]))
# im1 = ax1.imshow(im[0].data, origin='lower', vmin=600, vmax=800)
# ax1.plot([49.13674189,49.13062658,48.53523662,48.54543162]*u.deg, [ 41.39245501, 41.83609281,41.82997431, 41.38637844]*u.deg, color = 'r', lw=3, transform=ax1.get_transform("world"))


# grid = generate_FOV_grid(target_coord, FOV)

# for p in grid:
#     ax1.plot(p.ra, p.dec, marker='o', color='r', ms=10, transform=ax1.get_transform("world"))

# plt.show()
# exit()

#target_coord = SkyCoord(13.576921, -40.392550, unit=(u.deg, u.deg))
#s = survey_request('at2023pcw',target_coord , 'gri')
#s = survey_request(lco_fits[0].header['OBJECT'], target_coord, 'gri')
#s.search_for_SDSS()
#s.search_for_PS1()
#s.search_for_DECam()

galaxies = [[125.5433,-25.8880],
            [125.4069,-21.5691],
            [123.8377,-33.1189],
            [122.7733,-30.8752],
            [124.7261,-32.7621],
            [124.3686,-31.8556],
            [124.2324,-31.4309],
            [125.9628,-28.2488],
            [122.3471,-31.1373],
            [125.8190,-22.0851],
            [123.9586,-31.1844],
            [121.8064,-31.9581],
            [125.4043,-33.0136],
            [123.1610,-31.8733],
            [122.7131,-33.0211],
            [120.6195,-29.2095],
            [122.8184,-30.9066],
            [125.0188,-37.1348],
            [124.6821,-32.4162],
            [126.6102,-27.3605],
            [124.7649,-37.0665],
            [123.7104,-34.3874],
            [121.4886,-26.8309],
            [123.3382,-28.0374],
            [123.8712,-28.3101],
            [124.8092,-28.7278],
            [125.3822,-30.4291],
            [124.1169,-27.8063],
            [122.1070,-31.0502],
            [125.2081,-37.1875]]

url = 'https://alasky.cds.unistra.fr/hips-image-services/hips2fits#ra=119.51849999999999&dec=-27.298400000000004&fov=0.4&projection=AIT'

urllib.request.urlretrieve(url, 'test.html')

# for g in galaxies:
#     name = str(g[0])+'_'+str(g[1])
#     target_coord = SkyCoord(g[0], g[1], unit=(u.deg, u.deg))
#     s = survey_request(name, target_coord, 'gri')
    
#     #s.search_for_Skymapper()
#     s.search_for_DECam()


