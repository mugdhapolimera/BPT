# -*- coding: utf-8 -*-
"""
Created on Mon May 24 20:17:28 2021

@author: mugdhapolimera
"""

#import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
import pdb
import wget
import os.path
#sample = ['rf0358', 'rf0237', 'rf0064', 'rf0003', 'rs0010', 'rf0503','rs0124',
#'rf0477', 'rs0315', 'rs1036', 'rs1152', 'rs0107', 'rs0672', 'rs1004', 'rs1112', 
#'rs0181', 'rs1111', 'rs1106', 'rf0127', 'rs1150']
gals = pd.read_csv('samplegals.csv', header=0,sep=',')

name = gals['name']
RA = gals['radeg']
DEC = gals['dedeg']

folder = 'C:\Users\mugdhapolimera\github\BPT\sample_images/'

for i in range(0,len(RA)):
    if os.path.isfile(folder+ str(name[i]) + "_dr8_images.png"):
        print('image already exists')
        continue
    else:
        wget.download("http://legacysurvey.org/viewer/jpeg-cutout?ra=" + \
                      str(RA[i]) + "&dec=" + str(DEC[i]) + \
                      "&size=256&layer=dr8&pixscale=0.262&bands=grz", \
                      out= folder + str(name[i]) + "_dr8_images.png")