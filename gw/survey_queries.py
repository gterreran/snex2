# -*- coding: utf-8 -*-
"""survey_queries.py

"""

from astropy.table import Table
from astropy.io.votable import parse
from astropy import units as u
import urllib.request, io
from numpy.core.defchararray import startswith
from pyvo.dal import sia

#Folder definitions
OUT_FOLDER = '/Users/giacomoterreran/GW_templates/O4' #folder where images will be stored
#OUT_FOLDER = '/supernova/data/extdata/GW_templates/O4'

#PS1 urls
PS1_TAB = 'https://ps1images.stsci.edu/cgi-bin/ps1filenames.py'
PS1_CUT = 'https://ps1images.stsci.edu/cgi-bin/fitscut.cgi'

#List of available services at:
#https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/DataAccessInterfaces/SimpleImageAccessSIA/SimpleImageAccessSIA.html
COLLECTION = 'coadd_all'
DECAM_SERVICE = f'https://datalab.noirlab.edu/sia/{COLLECTION}'

#Skymapper urls
SKY_CUT = 'https://api.skymapper.nci.org.au/public/siap/dr4/query'

WIDTH_RA = 10/60 #deg
WIDTH_DEC = 10/60 #deg

class survey_request:
    '''
    This class will query the surveys looking for the optical templates.
    Each function is associated with a survey. First it will check if
    the coordinates are in the footprint of that survey. Each function
    will return True or False accordingly. In addition, if the coordinates
    are indeed in the footprint, then the function will store the hdu of
    the template in the 'hdu' attribute.

    '''
    def __init__(self, _ra, _dec, _filters):
        self.ra0 = _ra #deg
        self.dec0 = _dec #deg
        self.filters = _filters
        self.urls = {}
        self.templates_paths = {}

    def search_for_PS1_urls(self):
        '''
        PS1 needs the size in pixels, considering 0.25 arcsec/pixel.

        '''

        self.urls['PS1']={'g':'--', 'r':'--', 'i':'--'}

        tableurl = f'{PS1_TAB}?ra={self.ra0}&dec={self.dec0}&size={3600*WIDTH_RA/0.25}&format=fits&filters={self.filters}'
        table = Table.read(tableurl, format='ascii')
        for i in range(len(table)):
            flt = table[i]['filter']
            filename = table[i]['filename']
            self.urls['PS1'][flt] = f'{PS1_CUT}?ra={self.ra0}&dec={self.dec0}&size=2500&format=fits&red={filename}'


        
    def search_for_Skymapper_urls(self):
        '''
        Skymapper has the ability to search for multiple filters in one go,
        but the returned list of file loses the information on the filter,
        so I need to send a seperate quary for each filter.

        '''

        self.urls['Skymapper']={'g':'--', 'r':'--', 'i':'--'}

        for flt in self.filters:

            url = f'{SKY_CUT}?POS={self.ra0},{self.dec0}&SIZE={WIDTH_RA},{WIDTH_DEC}&BAND={flt}&FORMAT=image/fits&INTERSECT=covers&RESPONSEFORMAT=VOTABLE'

            u = urllib.request.urlopen(url)
            s = io.BytesIO(u.read())
            votable = parse(s)

            for resource in votable.resources:
                for t in resource.tables:
                    table = t.array['get_fits'][:]

            if len(table) != 0:
                self.urls['Skymapper'][flt] = table[0]


    def search_for_DECam_urls(self):

        self.urls['DECam']={'g':'--', 'r':'--', 'i':'--'}

        connect = sia.SIAService(DECAM_SERVICE)
        table = connect.search(pos = (self.ra0,self.dec0), size = (WIDTH_RA, WIDTH_DEC), verbosity=2).to_table()

        for flt in self.filters:

            sel = (table['prodtype'] == 'image') & (startswith(table['obs_bandpass'].astype(str), flt))
            
            if len(table[sel]) != 0:
                self.urls['DECam'][flt] = table[sel]['access_url'][0]
    
    def fetch_urls(self):
        '''
        It will also store the templates paths in a dictionary
        ready to be parsed by function populating the table
        '''
        for survey in self.urls:
            for flt in self.urls[survey]:
                if self.urls[survey][flt] == '--':
                    self.templates_paths[f'{survey}_{flt}'] = '--'
                else:
                    out_name = f'{OUT_FOLDER}/{survey}/{self.ra0:.4f}_{self.dec0:.4f}_{flt}.fits'
                    urllib.request.urlretrieve(self.urls[survey][flt], out_name)
                    self.templates_paths[f'{survey}_{flt}'] = out_name
