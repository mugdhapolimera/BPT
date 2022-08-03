# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 17:25:02 2020

@author: mugdhapolimera
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 15:25:30 2018

@author: mugdhapolimera

This code explores the properties of galaxies categorized using BPT plots.
"""

import numpy as np
import os
import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import ScalarFormatter
#from scipy.stats import norm

#path = os.path.dirname(os.path.realpath(__file__))
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
he2_flag = 0
full = 0
resolve = 1
eco = 0
sdsscat = 'jhu'#'nsa'
if he2_flag:
    flags = pd.read_csv('eco+resolve_emlineclass_filter_he2.csv')
else:
    flags = pd.read_csv('eco+resolve_emlineclass_filter.csv')
if full : 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_snr5_dext_nsa.csv'
    flags = pd.read_csv('eco+resolve_emlineclass_dext_snr5_nsa.csv')
if resolve: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5_dext_'+sdsscat+'.csv'
    flags = pd.read_csv('resolve_emlineclass_dext_snr5_'+sdsscat+'.csv')
    conf = pd.read_csv('RESOLVE_snr5_master_conf_new.csv')
if eco: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_filter_new.pkl'
    flags = pd.read_csv('eco_emlineclass_filter_new.csv')

full_df = pd.read_csv(inputfile)
full_df.index = full_df.name
df = full_df.loc[list(flags['galname'])]
if 'NAME' not in df.keys():
    df['NAME'] = df['name']

keys = ['defstarform', 'defagn', 'composite', 'agntosf', 'sftoagn']
if he2_flag:
    keys.append('heiisel')

    
marker = {'agntosf': 'g^', 'ambigagn': 'ms', 'composite': 'ms', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.','sftoagn': 'bs', 'sftoagn1': 's', 'sftoagn2': 'm*'}

colors = {'agntosf': 'g', 'ambigagn': 'm', 'composite': 'm', 'defagn': 'r', 
          'defliner': 'y', 'defseyf': 'c', 'heiisel': 'k',
          'defstarform': 'gray', 'sftoagn': 'b', 'sftoagn2': 'b'}

labels = {'agntosf': 'AGN-to-SF', 'ambigagn': 'Ambiguous AGN', 
          'composite': 'Composite', 'defagn': 'Definite AGN', 
          'defliner': 'LINER', 'defseyf': 'Seyfert', 
          'heiisel': 'HeII-Selected AGN', 'defstarform': 'Definite SF', 
          'sftoagn': 'SFing-AGN', 'sftoagn2' : 'MP-AGN2'}

percent = {'agntosf': 0, 'composite': 0, 'defagn': 0, 'heiisel': 57,
          'defstarform': 0, 'sftoagn': 0}

for key in keys:
    percent[key] = len(np.where(flags[key])[0])

df = df[df.logmstar > 0]
flags.index = flags.galname

sf= ['rs1126', 'rf0053', 'rs0262', 'rf0399']
sfingagn=['rs0010', 'rs0472', 'rf0503', 'rs1195']
he2=['rf0013', 'rf0372', 'rs0463', 'rs1214', 'rs1103']
composite=['rf0002','rs0261']
#agn=['rs0029','rs1116']
agn=['rs1111','rf0127']

jwst = sf+sfingagn+he2+composite+agn
odd = ['rs0582']
#G/S vs. M_star with histogram
'''rect_histy = [left_h, bottom, 0.2, height]
axScatter = plt.axes(rect_scatter)
axHisty = plt.axes(rect_histy, xscale = 'log')
nullfmt = NullFormatter() # no labels
axHisty.yaxis.set_major_formatter(nullfmt)
'''
plt.figure('Physical Properties')
#plt.suptitle(sdsscat+' catalog')
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]
axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx, yscale = 'log')
nullfmt = NullFormatter() # no labels
axHistx.xaxis.set_major_formatter(nullfmt)

# the scatter plot:
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        if key == 'defstarform':
            axScatter.plot(sel.logmstar,
                           10**sel.logmgas/10**sel.logmstar, 
                           marker[key], label = labels[key], markersize = 8, 
                           alpha = 0.3, mew = 0)
        
        #elif key == 'sftoagn':
        #    axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
        #                   label = labels[key], markersize = 15, mew = 0)
                
        elif key == 'agntosf':
            axScatter.plot(sel.logmstar,
                           10**sel.logmgas/10**sel.logmstar, 
                           marker[key], label = labels[key], markersize = 10, 
                           mew = 1, mec = 'y')            
        elif key == 'heiisel':
            axScatter.plot(sel.logmstar,
                           10**sel.logmgas/10**sel.logmstar, 
                           marker[key], label = labels[key], markersize = 12,  
                           mfc = 'None', mew = 2)            
        else:
            axScatter.plot(sel.logmstar,
                           10**sel.logmgas/10**sel.logmstar, 
                           marker[key], markersize = 12, mew = 0, 
                           label = labels[key])
        axScatter.plot(np.linspace(7.5,12.5), 
                       np.ones(len(np.linspace(7.5,12.5))), 'k-.')
        #axScatter.text(11.0, 0.005, 'RESOLVE', fontsize=14, color='k')
        axScatter.text(10.5, 1.1, r'1:1 G/S Ratio', fontsize=14, color='k')
        
        axScatter.text(9.52, -2.2, r'Gas Richness', fontsize=14, color='k')
        axScatter.text(9.52, -2.4, r'Threshold Mass', fontsize=14, color='k')

        axScatter.plot(9.5*np.ones(len(np.linspace(10**-2.5,10**1.5))), 
                       np.linspace(10**-2.5,10**1.5), 'k-.')
        axScatter.set_ylim(10**-2.5,10**1.5)
        axScatter.set_xlim(7.5,12.5) 
        axScatter.set_yscale("log")
        axScatter.yaxis.set_major_formatter(ScalarFormatter())
        axScatter.legend(loc='upper right',numpoints = 1, fontsize = 12) # bbox_to_anchor=(1,1),
        #if he2_flag:
            #axScatter.legend(loc=2, bbox_to_anchor=(1,1.15),numpoints = 1)

        axScatter.set_xlabel(r'$\rm \log(M_{stellar}/M_{\odot})$', fontsize = 22)
        axScatter.set_ylabel(r'$\rm M_{gas}/M_{stellar}$',fontsize = 22)
axScatter.plot(df.logmstar.loc[jwst], 10**df.logmgas.loc[jwst]/10**df.logmstar.loc[jwst], 
                           '*', color = 'lime',markersize = 12, mew = 0, 
                           label = labels[key])
axScatter.plot(df.logmstar.loc[odd], 10**df.logmgas.loc[odd]/10**df.logmstar.loc[odd], 
                           '*', color = 'black',markersize = 12, mew = 0, 
                           label = labels[key])
seyfert = list(flags.galname[flags.defstarform & (q.Estimate-np.log10(3e10) >= -3.0)])
for i, txt in enumerate(seyfert):
    axScatter.annotate(txt, (df.logmstar.loc[seyfert[i]], 
                             10**df.logmgas.loc[seyfert[i]]/10**df.logmstar.loc[seyfert[i]]), fontsize = 10)
#for i, txt in enumerate(jwst):
#    axScatter.annotate(txt, (df.logmstar.loc[jwst[i]], 
#                             10**df.logmgas.loc[jwst[i]]/10**df.logmstar.loc[jwst[i]]))

if not he2_flag:
    bins = np.arange(7.5,11.5,0.25)
    axHistx.hist(df.logmstar, histtype = 'stepfilled', alpha = 0.1,
             bins= bins, linewidth = 5, label = 'All Galaxies')
    mstars = df.iloc[np.where(flags['defstarform'])[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins, alpha = 0.3,
                 linewidth = 5, label = labels['defstarform'], 
                    color = colors['defstarform'])
    mstars = df.iloc[np.where((flags['defagn']) | (flags['defseyf']) | 
                    (flags['defliner']))[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins, hatch = '/',
                 linewidth = 5, label = 'Definite AGN', 
                 color = colors['defagn'])
    mstars = df.iloc[np.where((flags['sftoagn']))[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins, hatch = '\\',
                 linewidth = 5, label = 'SFing-AGN', 
                 color = colors['sftoagn'])
    axHistx.legend(loc='upper right', fontsize = 12) #bbox_to_anchor=(1,1.15), 
    
axHistx.set_ylabel('Number')
axHistx.set_xlim(axScatter.get_xlim())
axHistx.yaxis.set_major_formatter(ScalarFormatter())
 

lowsfagn = df[flags.sftoagn & (10**df.logmgas/10**df.logmstar < 1.0)]

ssfr_st = 10**(np.log10(df.sfr_nuv_wise) - df.logmstar)
ssfr_lt = 1/(1+(1/(df.meanfsmgr)))
fsmgr_st = 100*(10**6)*(ssfr_st)/(1-ssfr_st)
fsmgr_lt = df.meanfsmgr
plt.figure()
#plt.suptitle(sdsscat+' catalog')

for key in keys:
    sel = np.array(df.loc[flags[key]].name)
    if key == 'defstarform':
        plt.plot(fsmgr_st.loc[sel], fsmgr_lt.loc[sel],
            'k.', alpha = 0.3, label = labels[key])
    else:    
        plt.plot(fsmgr_st.loc[sel], fsmgr_lt.loc[sel],
            marker[key], ms = 10, label = labels[key])
xaxis = np.arange(np.min(fsmgr_st)-0.05, np.max(fsmgr_st)+0.05,0.01)
yaxis = np.ones(len(xaxis))
plt.plot(xaxis, yaxis, 'k--', lw = 3)
plt.xlim(np.min(fsmgr_st), np.max(fsmgr_st))
plt.ylim(np.min(fsmgr_lt), np.max(fsmgr_lt))
plt.text(0.0005, 1.25, r'Stellar Mass Doubled in last Gyr', 
             fontsize=14, color='k')
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 15, loc = 'lower right')
plt.plot(fsmgr_st.loc[jwst], fsmgr_lt.loc[jwst],
            '*',color = 'lime', ms = 10, label = labels[key])
plt.plot(fsmgr_st.loc[odd], fsmgr_lt.loc[odd],
            '*',color = 'black', ms = 10, label = labels[key])
plt.xlabel(r'Short Term SFH $\left(\frac{M_*(<100 Myr)}{M_*(>100 Myr)}\right)$', 
           fontsize = 22)
plt.ylabel(r'Long Term SFH $\left(\frac{M_*(<1 Gyr)}{M_*(>1 Gyr)}\right)$', 
           fontsize = 22)
dwarf =  seyfert#df.name[df.logmstar < 9.5]
#for i, txt in enumerate(df.name):
#    plt.annotate(txt, (fsmgr_st.loc[df.name[i]], 
#                             fsmgr_lt.loc[df.name[i]]))
for i, txt in enumerate(dwarf):
    plt.annotate(txt, (fsmgr_st.loc[dwarf[i]], 
                             fsmgr_lt.loc[dwarf[i]]), fontsize = 10)

plt.figure()
plt.suptitle('AGNFRAC = 0 only')

ax1 = plt.subplot()
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
Z = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_snr5_dext_jhu_csf_intrinsic_open_coincident_LOGZ.txt", 
                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
Z.index = Z.name
q = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_snr5_dext_jhu_csf_intrinsic_open_coincident_LOGQ.txt", 
                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
q.index = q.name

agnfrac = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_snr5_dext_jhu_csf_intrinsic_open_coincident_agn_AGNFRAC.txt", 
                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
agnfrac.index = agnfrac.name

Z_pp04 = Z.loc[df.name]['Estimate']+8.76
#Z_pp04[bad] = -99
for key in keys:
    if key == 'defstarform':
        ax1.plot(df.logmstar.loc[flags[key]],Z_pp04.loc[flags[key]],
                 'k.', alpha = 0.3, label = labels[key])
    else:
        ax1.plot(df.logmstar.loc[flags[key]], Z_pp04.loc[flags[key]],
         marker[key], markersize = 10, label = labels[key])
ax1.plot(df.logmstar.loc[jwst], Z_pp04.loc[jwst],
         '*', color = 'lime', markersize = 10, label = labels[key])
ax1.plot(df.logmstar.loc[odd], Z_pp04.loc[odd],
         '*', color = 'black', markersize = 10, label = labels[key])

yaxis = np.linspace(np.min(Z_pp04)-0.1,np.max(Z_pp04)+0.1)
xaxis = np.linspace(7.5, 11.5)
plt.plot(9.5*np.ones(len(yaxis)),yaxis, 'k-.', linewidth = 3)
Z04 = np.log10(0.4)+8.76
plt.plot(xaxis, Z04*np.ones(len(xaxis)), 'k--', linewidth = 2)
plt.plot(xaxis, 8.76*np.ones(len(xaxis)), 'k--', linewidth = 2)
#plt.ylim(min(yaxis),max(yaxis))
#plt.ylim(7.8,9.2)
#plt.xlim(7.5,11.5)
plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 22)
plt.ylabel(r'Z (12 + log(O/H))', fontsize = 22)
#plt.legend(fontsize = 12)
ax2 = ax1.twinx()
#ax2.set_xticks([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1.0])
yticks = np.array([7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0,9.2])
ax2.set_yticks(np.linspace(0,1,len(yticks)))
#ax1.set_yticks(yticks)#np.arange(7.8,9.2,0.2))
float_formatter = lambda x: "%.2f" % x
#xticks = np.array([7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2])
N2 = 1.754*yticks - 15.614
N2_label = ["%.2f" % z for z in N2]
ax2.set_yticklabels(N2_label)
ax2.set_ylabel(r'[NII]/H$\alpha$', fontsize = 22)
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\M-Z_Tremonti04.txt')
#ax1.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
#ax1.plot(m+10,z)

