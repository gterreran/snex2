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
# COLLECTION = 'coadd/decaps_dr2'
# COLLECTION = 'delve_dr2'
DECAM_SERVICE = f'https://datalab.noirlab.edu/tap/{COLLECTION}'

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

    # <td>Survey:</td>
    # <td><select name="survey" id="survey" class="survey-option" style="width:50%">
    #     <option value="ivoa_nsa.siav1">NOIRLab Astro Data Archive</option>
    #     <option value="ivoa_decaps_dr1.siav1">DECaPS DR1 image tiles</option>
    #     <option value="ivoa_des_dr1.siav1">DES DR1 image tiles</option>
    #     <option value="ivoa_des_dr2.siav1">DES DR2 image tiles</option>
    #     <option value="ivoa_des_sva1.siav1">DES SVA1 image tiles</option>
    #     <option value="ivoa_des_y1.siav1">DES Year 1 image tiles</option>
    #     <option value="ivoa_des_y2.siav1">DES Year 2 image tiles</option>
    #     <option value="ivoa_des_y3.siav1">DES Year 3 image tiles</option>
    #     <option value="ivoa_des_y4.siav1">DES Year 4 image tiles</option>
    #     <option value="ivoa_gogreen_dr1.siav1">GOGREEN DR1 image tiles</option>
    #     <option value="ivoa_ls_dr8.siav1">Legacy Survey DR8 image tiles</option>
    #     <option value="ivoa_ls_dr9.siav1">Legacy Survey DR9 image tiles</option>
    #     <option value="ivoa_nsc_dr2.siav1">NOIRLab Source Catalog DR2 image tiles</option>
    #     <option value="ivoa_sdss_dr9.siav1">SDSS DR9 image tiles</option>
    #     <option value="ivoa_smash_dr1.siav1">SMASH DR1 image tiles</option>
    #     <option value="ivoa_smash_dr2.siav1">SMASH DR2 image tiles</option>
    #     <option value="ivoa_splus_dr1.siav1">S-Plus DR1 image tiles</option>
    #     <option value="ivoa_splus_edr.siav1">S-Plus EDR image tiles</option>
    #     <option value="ivoa_stripe82.siav1">Stripe82 image tiles</option>
    #     <option value="ivoa_calibrated.calibrated_all">Instrument calibrated data only</option>
    #     <option value="ivoa_coadd.coadd_all">Image stacks/coadds only</option>
    #     <option value="ivoa_raw.raw_all">Raw data only</option>
    #     </select>
    # </td>



    'https://datalab.noirlab.edu/query/query?sql=" + survey_temp + "+q3c_radial_query(s_ra,+s_dec,+" + items[_i_ra] + ",+" + items[_i_dec] + ",+2)+)+SELECT+access_url,instrument_name,obs_bandpass,exptime,prodtype,proctype,date_obs+FROM+region+WHERE+(abs(spat_hilimit1+-+spat_lolimit1)+<+90.0+AND+("+items[_i_ra]+"+BETWEEN+spat_lolimit1+AND+spat_hilimit1)+AND+("+items[_i_dec]+"+BETWEEN+spat_lolimit2+AND+spat_hilimit2))+OR+(abs(spat_hilimit1+-+spat_lolimit1)+>+90.0+AND+(+("+items[_i_ra]+"+BETWEEN+0.0+AND+spat_lolimit1)+OR+("+items[_i_ra]+"+BETWEEN+spat_hilimit1+AND+360.0))+AND+("+items[_i_dec]+"+BETWEEN+spat_lolimit2+AND+spat_hilimit2))+ORDER+BY+date_obs+DESC&ofmt=csv&out=none&async=false";'


    connect = sia.SIAService(DECAM_SERVICE)
    table = connect.search(pos = (self.coord.ra.deg, self.coord.dec.deg), size = (LCO_INSTRUMENTS['sinistro'].fov.to(u.deg).value, LCO_INSTRUMENTS['sinistro'].fov.to(u.deg).value), verbosity=2).to_table()

    print('something')
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