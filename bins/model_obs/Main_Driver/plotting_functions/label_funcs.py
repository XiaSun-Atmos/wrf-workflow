#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A set of utilities for generating labeling strings for the plots
###############################################################################
@author: dmurray
"""
###############################################################################
from format_funcs import format_lat_lon_alt

"""
Create a label for the averaging minutes

Param:   avgmins    - minutes being averaged over

Return:  time related string
"""
def avgmins_to_label(avgmins=60):
    if avgmins == 60:
       avg_label='1Hr'
    elif avgmins == 30:
       avg_label='1/2Hr'
    else:
       avg_label=str(avgmins)+'Min'

    return avg_label

"""
Create a label for the grid point location:

Params:     INTPoi       - interpolation type
            sfxNum       - station id suffix number
            Grefs        - grid location reference index
            MDLlat_grid  - closest grid latitude
            MDLlon_grid  - closest grid longitude
            MDLelev_grid - closest grid altitude
            MDLlat_site  - interpolated grid latitude
            MDLlon_site  - interpolated grid longitude
            MEASmdlid    - id for the closest point model
            ALTmdlid     - id for the interpolated point model

Return:   A string descriptive of the grid point location
"""
def gridpoint_to_label(INTPoi,sfxNum,Grefs,MDLlat_grid,MDLlon_grid,MDLelev_grid,MDLlat_site,MDLlon_site,MEASmdlid,ALTmdlid):
    if (INTPoi == 0):
        if (sfxNum==0):
            if (len(Grefs) == 1):
               MDLgridinfo = 'Closest gridpoint: '+format_lat_lon_alt(MDLlat_grid,MDLlon_grid,MDLelev_grid)
            else:
               MDLgridinfo = 'Closest '+MEASmdlid+' gridpoint: '+format_lat_lon_alt(MDLlat_grid,MDLlon_grid,MDLelev_grid)
        else:
            MDLgridinfo = 'Closest '+ALTmdlid+' gridpoint: '+format_lat_lon_alt(MDLlat_grid,MDLlon_grid,MDLelev_grid)
    else:
        if (sfxNum==0):
            if (len(Grefs) == 1):
                MDLgridinfo = 'Interpolated to gridpoint: '+format_lat_lon_alt(MDLlat_site,MDLlon_site,None)
            else:
                MDLgridinfo = 'Interpolated to '+MEASmdlid+' gridpoint: '+format_lat_lon_alt(MDLlat_site,MDLlon_site,None)
        else:
            MDLgridinfo = 'Interpolated to '+ALTmdlid+' gridpoint: '+format_lat_lon_alt(MDLlat_site,MDLlon_site,None)

    return MDLgridinfo

