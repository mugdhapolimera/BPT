# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:10:34 2019

@author: mugdhapolimera

Marginalize Photoionization grid along Z and q dimensions to get overall AGN
distribution

"""

import numpy as np
import pandas as pd
import os
import sys
from scipy.integrate import cumtrapz, simps

#sys.path.append('C:\Users\mugdhapolimera\github\izi\izi_utils')
os.chdir('C:\Users\mugdhapolimera\github\izi')
from izi_utils.tabulate import idl_tabulate
from izi_utils.interpolate import grid_interpolate3d


def marginalize(grid, params = ['AGNFRAC', 'LOGZ', 'LOGQ'], 
                mar_dim = "AGNFRAC"):

    #grid, ngrid, zarr, qarr, agnarr, dlogz, dlogq, dloga = grid_interpolate3d(grid0, method = 'scipy', 
    #                                                nz1 = 50, nq1 = 50, na1 = 50)
    interp = grid_interpolate3d(grid)
    spacing = []
    for key in params: 
        spacing.append(np.diff(np.unique(interp[key]))[0])
    lines = [x for x in interp.keys() if x not in params]
    int_axis = [x for x in params if x not in mar_dim]   
    lowZmar_grid = pd.DataFrame(np.unique(interp['AGNFRAC']), columns = ["AGNFRAC"])
    highZmar_grid = pd.DataFrame(np.unique(interp['AGNFRAC']), columns = ["AGNFRAC"])
    for key in lines:
        
        nlowZ = len(np.unique(interp.LOGZ[interp.LOGZ < np.log10(0.4)]))
        nhighZ = len(np.unique(interp.LOGZ[interp.LOGZ >= np.log10(0.4)]))
                
        lowZ_interp = interp[interp.LOGZ < np.log10(0.4)]
        highZ_interp = interp[interp.LOGZ > np.log10(0.4)]

        index = params.index(int_axis[0])
        lowZmarginalized_2D = np.trapz(np.array(lowZ_interp[key]).reshape(50,nlowZ,50)
                                            ,axis=2, dx=spacing[index])
        index = params.index(int_axis[1])
        lowZmarginalized_1D = np.trapz(lowZmarginalized_2D, axis=1, dx=spacing[index])
        #integral = simps(lowZmarginalized_1D, dx=spacing[0])
        #marginalized_1D /= integral
        lowZmar_grid[key] = lowZmarginalized_1D

        index = params.index(int_axis[0])
        highZmarginalized_2D = np.trapz(np.array(highZ_interp[key]).reshape(50,nhighZ,50)
                                            ,axis=2, dx=spacing[index])
        index = params.index(int_axis[1])
        highZmarginalized_1D = np.trapz(highZmarginalized_2D, axis=1, dx=spacing[index])
        highZmar_grid[key] = highZmarginalized_1D
    
    return lowZmar_grid, highZmar_grid
