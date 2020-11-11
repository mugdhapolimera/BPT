# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 09:16:33 2019

@author: mugdhapolimera

Multiwavelength AGN crossmatch
"""

import pandas as pd
import numpy as np
import os
#SNR5 master catalog
os.chdir('C:\Users\mugdhapolimera\github\SDSS_Spectra')

resolve = pd.read_csv('RESOLVE_full_blend_dext_new.csv')
resolve.index = resolve.name

df = pd.read_csv('RESOLVE_snr5_master.csv')
df.index = df.name
total = len(df)
dwarfno = np.sum(df.logmstar < 9.5)
print('Total , Dwarf')
print(total, dwarfno)

##############################################################################
## Mid-IR ##
##############################################################################
midir = pd.read_csv('mid_ir\WISE_good.csv')
midir.index = midir.name
midiragn = pd.read_csv('mid_ir\WISE_Mid_IR-AGN.csv')
midiragn.index = midiragn.name
midiragn_sn5 = [x for x in list(df.name) if x in list(midiragn.name)]
midir_sn5 = [x for x in list(df.name) if x in list(midir.name)]

no_midir= len(midir_sn5)
no_midirdwarf = np.sum(df.loc[midir_sn5].logmstar < 9.5)
no_midiragn = len(midiragn_sn5)
no_midirdwarfagn = np.sum(df.loc[midiragn_sn5].logmstar < 9.5)

print(no_midiragn, no_midir, no_midirdwarfagn, no_midirdwarf)

##############################################################################
## Optical Emission Line##
##############################################################################
opt = pd.read_csv('resolve_emlineclass_full_snr5_master.csv')
opt.index = opt.galname
optagn = opt[opt.sftoagn | opt.defagn | opt.composite]
optagn_sn5 = [x for x in list(df.name) if x in list(optagn.galname)]
opt_sn5 = [x for x in list(df.name) if x in list(opt.galname)]

no_opt= len(opt_sn5)
no_optdwarf = np.sum(df.loc[opt_sn5].logmstar < 9.5)
no_optagn = len(optagn_sn5)
no_optdwarfagn = np.sum(df.loc[optagn_sn5].logmstar < 9.5)

print(no_optagn, no_opt, no_optdwarfagn, no_optdwarf)

##############################################################################
## Optical Broad Line##
##############################################################################
broad = pd.read_csv(r'C:\Users\mugdhapolimera\github\xray\catalog_matching\BroadlineAGN_RESOLVE.csv')
broadagn = broad['eco+res_name']
broadagn_sn5 = [x for x in list(df.name) if x in list(broadagn)]
broad_sn5 = list(df.name)

no_broad= len(broad_sn5)
no_broaddwarf = np.sum(df.loc[broad_sn5].logmstar < 9.5)
no_broadagn = len(broadagn_sn5)
no_broaddwarfagn = np.sum(df.loc[broadagn_sn5].logmstar < 9.5)

print(no_broadagn, no_broad, no_broaddwarfagn, no_broaddwarf)

##############################################################################
## Veron-Cetty 2003 Catalog##
##############################################################################
cat = pd.read_csv(r'C:\Users\mugdhapolimera\github\xray\catalog_matching\VeronAgnMatched.csv')
catagn = []
for i in range(len(cat)):
    name = cat['eco+res_name'][i]
    if name[0] == 'E':
        name = resolve.index.values[np.where(resolve.econame == name)]
    if len(name) > 1:
        catagn.append(name)
#print(catagn)
catagn_sn5 = [x for x in list(df.name) if x in list(catagn)]
cat_sn5 = list(df.name)

no_cat= len(cat_sn5)
no_catdwarf = np.sum(df.loc[cat_sn5].logmstar < 9.5)
no_catagn = len(catagn_sn5)
no_catdwarfagn = np.sum(df.loc[catagn_sn5].logmstar < 9.5)

print(no_catagn, no_cat, no_catdwarfagn, no_catdwarf)

##############################################################################
## Half Million Quasars 2015 Catalog##
##############################################################################
hmq = pd.read_csv(r'C:\Users\mugdhapolimera\github\xray\catalog_matching\HMQAgnMatched.csv')
hmqagn = []
for i in range(len(hmq)):
    name = hmq['eco+res_name'][i]
    if name[0] == 'E':
        name = resolve.index.values[np.where(resolve.econame == name)]
    if len(name) > 1:
        hmqagn.append(name)
#print(hmqagn)
hmqagn_sn5 = [x for x in list(df.name) if x in list(hmqagn)]
hmq_sn5 = list(df.name)

no_hmq= len(hmq_sn5)
no_hmqdwarf = np.sum(df.loc[hmq_sn5].logmstar < 9.5)
no_hmqagn = len(hmqagn_sn5)
no_hmqdwarfagn = np.sum(df.loc[hmqagn_sn5].logmstar < 9.5)

print(no_hmqagn, no_hmq, no_hmqdwarfagn, no_hmqdwarf)

##############################################################################
## X-ray - Newton-XMM3 ##
##############################################################################
xray_full = pd.read_csv(r'C:\Users\mugdhapolimera\github\xray\ECO+RESOLVE_xray_new.csv')
xray_full.index = xray_full.name
xray = pd.DataFrame({})
for i in range(len(xray_full)):
    flag = 0
    name = [xray_full['name'][i]]
    if name[0][0] == 'E':
        flag = 1
        oldname = name[0]
        name = resolve.index.values[np.where(resolve.econame == name[0])]
    if flag == 0:
        oldname = name[0]
    if len(name) > 0:
        #print(oldname, name[0])
        xray = xray.append(xray_full.loc[oldname])
        xray['name'][np.where(xray.name == oldname)[0][0]] = name[0]
xrayagn = xray[xray.xrayagn == 1.0]
xrayagn.index = xrayagn.name
xrayagn_sn5 = [x for x in list(df.name) if x in list(xrayagn.name)]
xray_sn5 = [x for x in list(df.name) if x in list(xray.name)]

no_xray= len(xray_sn5)
no_xraydwarf = np.sum(df.loc[xray_sn5].logmstar < 9.5)
no_xrayagn = len(xrayagn_sn5)
no_xraydwarfagn = np.sum(df.loc[xrayagn_sn5].logmstar < 9.5)

print(no_xrayagn, no_xray, no_xraydwarfagn, no_xraydwarf)

chandra = pd.read_csv(r'C:\Users\mugdhapolimera\github\xray\Chandra_crossmatch.csv',
                      comment = '#')
chandraname = [x.strip() for x in np.unique(chandra.usrid)]
chandra_sn5 = [x for x in list(df.name) if x in chandraname]