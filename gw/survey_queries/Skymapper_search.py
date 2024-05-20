# -*- coding: utf-8 -*-
"""Skymapper_search.py

"""

from astropy.io.votable import parse
import io

# Skymapper urls
SKY_CUT = 'https://api.skymapper.nci.org.au/public/siap/dr4/query'

@set_survey('Skymapper')        
def search_for_Skymapper(self, survey):
    '''
    Skymapper has the ability to search for multiple filters in one go,
    but the returned list of file loses the information on the filter,
    so I need to send a separate query for each filter.
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