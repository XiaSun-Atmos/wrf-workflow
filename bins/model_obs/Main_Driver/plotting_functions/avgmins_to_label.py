#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to convert averaging minutes to a label
###############################################################################
@author: dmurray
"""
###############################################################################

def avgmins_to_label(avgmins=60):
    if avgmins == 60:
       avg_label='1Hr'
    elif avgmins == 30:
       avg_label='1/2Hr'
    else:
       avg_label=str(avgmins)+'Min'

    return avg_label
