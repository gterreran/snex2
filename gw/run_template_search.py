# -*- coding: utf-8 -*-
"""run_template_search.py

Script to search and (if exist) download template
images from available optical surveys. The pixel scale, and the size
of the resulting template will match the characteristics of the LCO
specified.

The script is based on the :survey_queries: package. This is just
a wrapper to basically run the :query: class for a list of objects.

"""

from survey_queries.query import template_query


def run_template_search(_targets_list, _filters, _surveys, _instruments):


    out_dict ={}
    for target_id, event_id, name, ra, dec in _targets_list:
        s = template_query(target_id, event_id, name, ra, dec, _filters, _instruments)
        if 'PS1' in _surveys:
            s.search_for_PS1()
        if 'SDSS' in _surveys:
            s.search_for_SDSS()
        
        out_dict[target_id] = s
        
    return out_dict


