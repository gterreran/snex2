# -*- coding: utf-8 -*-
"""survey_queries.py

"""

from astropy.time import Time
from astropy.wcs import WCS
import astropy.io.fits as pf
import numpy as np

from . import PanSTARRS_search, SDSS_search, DECam_search, generate_FOV_grid, LCO_INSTRUMENTS, SURVEYS, OUT_FOLDER

class template_query:
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
                to be valid, the name will be checked if it is in the
                dictionary :LCO_INSTRUMENTS:, defined above.
        :type _inst: string or list

        '''
        self.obj = _obj
        self.coord = _coord
        self.filters = _filters

        # Making the :instruments: list and checking the instruments exist.
        if isinstance(_inst, str):
            if _inst == 'all':
                self.instruments = list(LCO_INSTRUMENTS.keys())
            else:
                self.instruments = [_inst]
        else:
            self.instruments = _inst
        for i in self.instruments:
            if i not in LCO_INSTRUMENTS.keys():
                raise ValueError(f'Instrument `{i}` does not appear to be currently supported. If this is a new instrument, you will need to included it in `LCO_INSTRUMENTS`.')
        
        # Sorting the instruments in decreasing FOV size.
        self.instruments.sort(key = lambda x: LCO_INSTRUMENTS[x].fov, reverse=True)

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
    
    @set_survey('PS1')
    def search_for_PS1(self, survey):
        data, weights = PanSTARRS_search.search_for_PS1(self, survey)
        self.make_templates(_data = data, _weights = weights, _survey = survey)
        print('This is done')

    @set_survey('SDSS')
    def search_for_SDSS(self, survey):
        data, weights = SDSS_search.search_for_SDSS(self, survey)
        self.make_templates(_data = data, _weights = weights, _survey = survey)
    
    @set_survey('DECam')
    def search_for_DECam(self, survey):
        data, weights = DECam_search.search_for_DECam(self, survey)
        self.make_templates(_data = data, _weights = weights, _survey = survey)
    

    def make_templates(self, _data, _weights, _survey):

        from reproject.mosaicking import find_optimal_celestial_wcs
        from reproject import reproject_interp
        from reproject.mosaicking import reproject_and_coadd

        # We now trim the list of skycells according to each instrument FOV
        # For the first instrument this step should be redundant since we
        # selected the skycell in order to tile the FOV already.
        for inst in self.instruments:
            new_grid = generate_FOV_grid(_center = self.coord, _fov = LCO_INSTRUMENTS[inst].fov, _step=SURVEYS[_survey].skycell_size)
            for flt in _data:
                # Since enumerate() returns a generator and generators can't be reversed,
                # we need to convert it to a list first
                for skycell_idx, skycell in reversed(list(enumerate(_data[flt]))):
                    contained = False
                    for p in new_grid:
                        if p.contained_by(WCS(skycell.header)):
                            contained = True
                            break
                    # If no point of the grid is contained by the skycell,
                    # drop it from all the data lists
                    if not contained:
                        _data[flt].pop(skycell_idx)
                        _weights[flt].pop(skycell_idx)
                
                print(f'Reprojecting and combining skycells of filer {flt} for {inst}.')

                wcs_out, shape_out = find_optimal_celestial_wcs(_data[flt], resolution=LCO_INSTRUMENTS[inst].resolution)
                array, footprint = reproject_and_coadd(_data[flt], wcs_out, shape_out=shape_out, reproject_function=reproject_interp, input_weights = _weights[flt], combine_function='mean')
                header = wcs_out.to_header()

                header = SURVEYS[_survey].track_skycells_in_header(_header = header, _data = _data[flt])

                # Populating the the header with the keys needed by the pipeline
                header['RA'] = self.coord.ra.deg
                header['DEC'] = self.coord.dec.deg
                header['CAT-RA'] = self.coord.ra.deg
                header['CAT-DEC'] = self.coord.dec.deg
                header['MJD-OBS'] = (np.mean([im.header['MJD-OBS'] for im in _data[flt]]), 'Average MJD')
                header['DATE-OBS'] = 'T'.join(Time(header['MJD-OBS'], format='mjd', scale='utc').iso.split())
                header['DAY-OBS'] = header['DATE-OBS'].split('T')[0].replace('-', '')
                header['TELESCOP'] = SURVEYS[_survey].telescope
                header['INSTRUME'] = SURVEYS[_survey].instrument
                header['OBJECT'] = self.obj
                header['FILTER'] = flt+'p' # NOTE: this works only with griz
                header['AIRMASS'] = 1
                header['WCSERR'] = 0
                header['PIXSCALE'] = (LCO_INSTRUMENTS[inst].resolution.value, 'arcsec/pixel')
                header['BUNIT'] = ('counts', 'ADUs')

                # Writing out the final combined template
                out_name = f'{OUT_FOLDER}/{_survey}/{inst}_{self.obj}_{flt}'
                self.templates_paths[f'{inst}_{_survey}_{flt}'] = out_name+'.fits'
                pf.writeto(out_name+'.fits', array, header, overwrite=True)

                # Hotpants needs the variance, so we need to write out 1/footprint
                header['BUNIT'] = ('counts2', 'ADUs^2')
                pf.writeto(out_name+'.var.fits', 1/footprint, header, overwrite=True)


