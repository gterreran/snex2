# -*- coding: utf-8 -*-
"""run_template_search.py

Script to search and (if exist) download template
images from available optical surveys. The pixel scale, and the size
of the resulting template will match the characteristics of the LCO
instrument specified. The objects will also be inserted in the
`o4_galaxy` table in the SNEx1 database, as well as the path to the
downloaded files. Note that the root folder of the download is
specified in the :__init__.py: of :survey_queries:.

The :dramatiq: decorator makes it possible to run it asynchronously.

"""

from .survey_queries.query import template_query
from custom_code.hooks import _return_session, _load_table
import logging
from django.conf import settings
from tom_targets.models import Target
import dramatiq
import datetime

TABLE_NAME = 'o4_galaxies'
DEFAULT_FILTERS = 'gri'
DEFAULT_SURVEYS = ['PS1','SDSS']
DEFAULT_INSTRUMENTS = 'all'

logger = logging.getLogger(__name__)

@dramatiq.actor(max_retries=0)
def search_templates_and_update_snex1(_targets_list, _wrapped_session, _filters=DEFAULT_FILTERS, _surveys=DEFAULT_SURVEYS, _instruments=DEFAULT_INSTRUMENTS):

    if _wrapped_session:
        db_session = _wrapped_session

    else:
        db_session = _return_session(settings.SNEX1_DB_URL)

    gw_table = _load_table(TABLE_NAME, db_address=settings.SNEX1_DB_URL)
    photlco = _load_table('photlco', db_address=settings.SNEX1_DB_URL)

    for target_id, event_id in _targets_list:

        existing_target = db_session.query(gw_table).filter(gw_table.targetid==target_id)
        if existing_target.count() > 0:
            if any([t.event_id == event_id for t in existing_target]):
                logger.info('Already ingested target {} into {} table for event {}'.format(target_id, TABLE_NAME, event_id))

            else:
                logger.info('Found existing target {} in {} table, copying it'.format(target_id, TABLE_NAME))
                existing_target_row = existing_target.first()
                existing_table = existing_target_row.__table__ #Somehow this is different than the o4_galaxies object
                
                non_pk_columns = [k for k in existing_table.columns.keys() if k not in existing_table.primary_key.columns.keys()]
                data = {c: getattr(existing_target_row, c) for c in non_pk_columns}
                data['event_id'] = event_id

                db_session.add(gw_table(**data))

        else:
            t = Target.objects.get(id=target_id)

            s = template_query(t.name, t.ra, t.dec, _filters, _instruments)

            if 'PS1' in _surveys:
                s.search_for_PS1()
            if 'SDSS' in _surveys:
                s.search_for_SDSS()
            
            db_session.add(gw_table(targetid = target_id, event_id = event_id, ra0=t.ra, dec0=t.dec, **s.templates_paths))

            for template in s.templates_paths:
                db_session.add(photlco(**photlco_dict(target_id, s.templates_paths[template])))


    if not _wrapped_session:
        try:
            db_session.commit()
        except:
            db_session.rollback()
        finally:
            db_session.close()

    else:
        db_session.flush()

    logger.info('Finished ingesting target {} into {}'.format(target_id, TABLE_NAME))

def photlco_dict(_targetid, _template_path):
    
    from survey_queries import SURVEY

    templeate_header = pf.getheader(_template_path)
    path_list = _template_path.spit('/')

    now = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M:%S.%f")

    phot_lco_dic = {
        'targetid': _targetid,
        'objname': templeate_header['OBJECT'],
        'dayobs': templeate_header['DAY-OBS'],
        'dateobs': templeate_header['DATE-OBS'].split('T')[0],
        'ut': templeate_header['DATE-OBS'].split('T')[1],
        'mjd': templeate_header['MJD-OBS'],
        'exptime':1,
        'filter': templeate_header['FILTER'],
        'telescopeid': SURVEY[templeate_header['SURVEY']].telescopeid,
        'instrumentid': SURVEY[templeate_header['SURVEY']].instrumentid,
        'telescope': templeate_header['TELESCOP'],
        'instrument': templeate_header['INSTRUME'],
        'mag': 9999,
        'dmag': 9999,
        'airmass': 1,
        'wcs': 0,
        'psf': 'X',
        'apmag': 9999,
        'psfx': 9999,
        'psfy': 9999,
        'psfmag': 9999,
        'psfdmag': 9999,
        'z1': 9999,
        'z2': 9999,
        'zn': 9999,
        'c1': 9999,
        'c2': 9999,
        'dz1': 9999,
        'dz2': 9999,
        'dc1': 9999,
        'dc2': 9999,
        'quality': 127,
        'zcat': 'X',
        'abscat': 'X',
        'fwhm': 9999,
        'magtype': 1,
        'ra0': templeate_header['RA'],
        'dec0': templeate_header['DEC'],
        'tracknumber': 0,
        'filename': path_list[-1],
        'filetype': 4,
        'filepath': '/'.join(path_list[:-1])+'/',
        'groupidcode': 2199023255552,
        'datecreated': now,
        'lastmodified': now,
        'apflux': 9999,
        'dapflux': 9999,
        'dapmag': 9999,
        'lastunpacked': now
    }

    return phot_lco_dic
