#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 20:09:06 2018

@author: mugpol
"""

import numpy as np
import os
from astropy.table import Table
import pandas as pd
import pickle

os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/BPT/')

filename = 'RESOLVE_all.fits'
infile =Table.read(filename, format = 'fits')
