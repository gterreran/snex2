from astropy import units as u
import os
import numpy as np

from . import PanSTARRS_search
from . import SDSS_search

# Folder definitions
OUT_FOLDER = os.getenv('GWTEMPLATES_FOLDER')
#OUT_FOLDER = '/supernova/data/extdata/GW_templates/O4'

class Instrument:
    def __init__(self, _fov, _resolution):
        self.fov = _fov
        self.resolution = _resolution

class Survey:
    def __init__(self, _telescope, _telescopeid, _instrument, _instrumentid, _skycell_size, _track_skycells_in_header):
        self.telescope = _telescope
        self.telescopeid = _telescopeid
        self.instrument = _instrument
        self.instrumentid = _instrumentid
        self.skycell_size = _skycell_size
        self.track_skycells_in_header = _track_skycells_in_header

# Specifications for different instruments
# The FOV here is thought as the diameter of a circle that would circumscribe the LCO image.
# This is effectively bigger than the actual FOV of each instrument, but it will guarantee
# that the full FOV is covered given any rotation 
LCO_INSTRUMENTS = {
    'sinistro': Instrument(_fov = 26.5*u.arcmin, _resolution = 0.389*u.arcsec),
    'muscat'  : Instrument(_fov = 9.1*u.arcmin, _resolution = 0.27*u.arcsec),
    # The FOV of 'qhy' is rectangular, so we calculate the side of a square having the same diagonal.
    'qhy'     : Instrument(_fov = np.sqrt((1.9**2 + 1.2**2)/2)*u.deg, _resolution = 0.74*u.arcsec)
}

SURVEYS ={
    'PS1': Survey(_telescope = 'PS1', _telescopeid = 41, _instrument = 'GPC1', _instrumentid = 191, _skycell_size = 0.4*u.deg, _track_skycells_in_header = PanSTARRS_search.track_skycells_in_header),
    'SDSS': Survey(_telescope = 'APO 2.5m', _telescopeid = 139, _instrument = 'SDSS', _instrumentid = 42, _skycell_size = 10*u.arcmin, _track_skycells_in_header = SDSS_search.track_skycells_in_header)
}

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
    from astropy.coordinates import SkyCoord
    import math

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