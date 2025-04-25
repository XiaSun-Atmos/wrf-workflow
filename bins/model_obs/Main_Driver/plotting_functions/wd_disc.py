#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
The Purpose of this Funciton is to Plot the Wind Direction Data so that Sharp 
Discontinuities (i.e. those between the 0/360 line) Are Not Visualized
"""

import numpy as np

###############################################################################
def wd_plt(Xtime,WD,plt_ax,line_type,line_color,line_label,marker_size=1,line_width=0):
    with np.errstate(invalid = 'ignore'):
            Dpts = np.where(abs(np.diff(WD))>180)[0]
    ###########################################################################
    if ((len(Dpts) != 0) & False) :  ## for now, ignore the differences
        #######################################################################
        # Account for Indice Convention
        Dpts = Dpts + 1;
        #######################################################################
        plt_ax.plot(Xtime[0:Dpts[0]],WD[0:Dpts[0]],line_type,color=line_color,markersize=marker_size,\
                    linewidth=line_width,label=line_label)
        plt_ax.plot(Xtime[Dpts[-1]:-1],WD[Dpts[-1]:-1],line_type,color=line_color,markersize=marker_size,\
                    linewidth=line_width,)
        #######################################################################
        for d in np.arange(0,len(Dpts)-1):
            if len(Xtime[Dpts[d]:Dpts[d+1]]) == 1:
                xwd = np.arange(Xtime[Dpts[d]]-900,Xtime[Dpts[d]]+900)
                ywd = np.ones(xwd.shape)*WD[Dpts[d]]
                plt_ax.plot(xwd,ywd,line_type,color=line_color,markersize=marker_size,linewidth=line_width);
                ###############################################################
                #describing the half-hour period where only one data is available
            else:
                plt_ax.plot(Xtime[Dpts[d]:Dpts[d+1]],WD[Dpts[d]:Dpts[d+1]],line_type,\
                            color=line_color,markersize=marker_size,linewidth=line_width)
            ###################################################################
    else:
       #plot using traditional methods
       plt_ax.plot(Xtime,WD,line_type,color=line_color,markersize=marker_size,linewidth=line_width,label=line_label)

    ###########################################################################
    return plt_ax
    
    
   
    
 
