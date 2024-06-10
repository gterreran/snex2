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
from ..custom_code.hooks import _return_session, _load_table
import logging
from django.conf import settings
from tom_targets.models import Target
import dramatiq


TABLE_NAME = 'o4_galaxies'

logger = logging.getLogger(__name__)

@dramatiq.actor(max_retries=0)
def search_templates_and_update_snex1(_targets_list, _filters, _surveys, _instruments, _wrapped_session):

    if _wrapped_session:
        db_session = _wrapped_session

    else:
        db_session = _return_session(settings.SNEX1_DB_URL)

    gw_table = _load_table(TABLE_NAME, db_address=settings.SNEX1_DB_URL)

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

            s = template_query(t.objname, t.ra, t.dec, _filters, _instruments)

            if 'PS1' in _surveys:
                s.search_for_PS1()
            if 'SDSS' in _surveys:
                s.search_for_SDSS()
            
            db_session.add(gw_table(targetid = target_id, event_id = event_id, ra0=t.ra, dec0=t.dec, **s.templates_paths))

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