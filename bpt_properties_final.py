# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 15:25:30 2018

@author: mugdhapolimera

This code explores the properties of galaxies categorized using BPT plots.
"""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

path = os.path.dirname(os.path.realpath(__file__))
he2_flag = 0
if he2_flag:
    flags = pd.read_csv('resolve_emlineclass_filtered_he2.csv')
else:
    flags = pd.read_csv('resolve_emlineclass_filtered.csv')

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_filtered.pkl'
full_df = pd.read_pickle(inputfile)
df = full_df.loc[flags.galname]

results_izi = 'C:/Users/mugdhapolimera/github/izi/results/IZI_Z_filtered.txt'
full_Z = np.genfromtxt(results_izi, dtype = None, names= True)
Z = full_Z[list(x for x in range(len(full_Z)) if full_Z['Name'][x] in list(flags.galname))]
        

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_izi_filtered/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

results_nb2 = pd.read_csv("results_k13_deredden/RESOLVE_param_estimates.csv")
Z_index2 = (results_nb2['Parameter'] == 'LOGZ')
full_Z_nb2 = results_nb2[Z_index2]
full_Z_nb2.index = full_Z_nb2['Galaxy Name']
Z_nb2 = full_Z_nb2.loc[flags.galname]

results_izi = 'C:/Users/mugdhapolimera/github/izi/results/IZI_Z_prior2.txt'
full_Z1 = np.genfromtxt(results_izi, dtype = None, names= True)
Z1 = full_Z1[list(x for x in range(len(full_Z1)) if full_Z1['Name'][x] in list(flags.galname))]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}

colors = {'agntosf': 'g', 'ambigagn': 'r', 'composite': 'b', 'defagn': 'r', 
          'defliner': 'y', 'defseyf': 'c', 'heiisel': 'k',
          'defstarform': 'k', 'sftoagn': 'm'}

labels = {'agntosf': 'AGN -> SF', 'ambigagn': 'Ambiguous AGN', 
          'composite': 'Composite', 'defagn': 'Definite AGN', 
          'defliner': 'LINER', 'defseyf': 'Seyfert', 
          'heiisel': 'HeII-Selected AGN', 'defstarform': 'Definite SF', 
          'sftoagn': 'MP-AGN'}

percent = {'agntosf': 0, 'composite': 0, 'defagn': 0, 'heiisel': 57,
          'defstarform': 0, 'sftoagn': 0}

for key in keys:
    percent[key] = len(np.where(flags[key])[0])

print percent
bins = np.arange(0.775,1.2,0.025)
plt.figure()
plt.hist(df.logmgas/df.logmstar, bins = bins, alpha = 0.1, 
         histtype = 'stepfilled', label = 'All Galaxies')
for key in keys:
        mgas = df.iloc[np.where(flags[key])[0]].logmgas
        mstars = df.iloc[np.where(flags[key])[0]].logmstar
        
        plt.hist(mgas/mstars, histtype = 'step', bins = bins,
                 linewidth = 5, label = labels[key],color = colors[key])
        plt.legend()
        plt.yscale('log')
        plt.xlabel('Gas/Stellar Mass Ratio')
        plt.ylabel('Number')

plt.figure(8)
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        if key == 'defstarform':
            plt.plot(sel.logmstar,sel.logmgas, marker[key], label = labels[key], 
                     alpha = 0.3)
        elif key == 'heiisel':
            plt.plot(sel.logmstar,sel.logmgas, marker[key], 
                           label = labels[key], markersize = 8,  mfc = 'None', mew = 2)            
        else:
            plt.plot(sel.logmstar,sel.logmgas, marker[key],
                     label = labels[key])
        plt.plot(np.linspace(7.5,11.5), np.linspace(7.5,11.5), 'k-.')
        plt.legend(loc = 4, numpoints = 1)        
        if he2_flag:
            plt.legend(loc = 4, numpoints = 1)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Gas Mass')

left, width = 0.1, 0.73
bottom, height = 0.1, 0.65
bottom_h = left_h = left + 0.65 + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
#rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(9)

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx, yscale = 'log')
#axHisty = plt.axes(rect_histy, xscale = 'log')

# no labels
nullfmt = NullFormatter()
axHistx.xaxis.set_major_formatter(nullfmt)
#axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        if key == 'defstarform':
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
                           label = labels[key], markersize = 8, alpha = 0.3)
        elif key == 'heiisel':
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
                           label = labels[key], markersize = 8,  mfc = 'None', mew = 2)            
        else:
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key],
                     markersize = 8, label = labels[key])
        axScatter.plot(np.linspace(7.5,11.5), np.linspace(7.5,11.5), 'k-.')
        axScatter.set_ylim(7.5,10.5)        
        if he2_flag:
            axScatter.legend(loc = 4, numpoints = 1)
        axScatter.set_xlabel('Stellar Mass')
        axScatter.set_ylabel('Gas Mass')
'''
#if not he2_flag:

plt.figure(9)
bins = np.arange(7.5,11.5,0.25)
plt.hist(df.logmstar, histtype = 'stepfilled', alpha = 0.1,
         bins= bins, linewidth = 5, label = 'All Galaxies')
mstars = df.iloc[np.where(flags['defstarform'])[0]].logmstar
plt.hist(mstars, histtype = 'step',bins = bins, alpha = 0.3,
             linewidth = 5, label = labels['defstarform'], color = colors[key])
mstars = df.iloc[np.where((flags['defagn']) | (flags['defseyf']) | 
                (flags['defliner']))[0]].logmstar
plt.hist(mstars, histtype = 'step',bins = bins,
             linewidth = 5, label = 'Traditional AGN')
mstars = df.iloc[np.where((flags['sftoagn']))[0]].logmstar
plt.hist(mstars, histtype = 'step',bins = bins,
             linewidth = 5, label = 'Non-Traditional AGN')
plt.legend(loc=2, bbox_to_anchor=(1,1.15))
'''
if not he2_flag:
    bins = np.arange(7.5,11.5,0.25)
    axHistx.hist(df.logmstar, histtype = 'stepfilled', alpha = 0.1,
             bins= bins, linewidth = 5, label = 'All Galaxies')
    mstars = df.iloc[np.where(flags['defstarform'])[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins, alpha = 0.3,
                 linewidth = 5, label = labels['defstarform'], color = colors[key])
    mstars = df.iloc[np.where((flags['defagn']) | (flags['defseyf']) | 
                    (flags['defliner']))[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins,
                 linewidth = 5, label = 'Traditional AGN', color = 'r')
    mstars = df.iloc[np.where((flags['sftoagn']))[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins,
                 linewidth = 5, label = 'Non-Traditional AGN', color = 'm')
    axHistx.legend(loc=2, bbox_to_anchor=(1,1.15))
#axHistx.set_xlabel('Stellar Mass')
axHistx.set_ylabel('Number')

'''bins = np.arange(7.5,11.5,0.25)
axHisty.hist(df.logmgas, histtype = 'stepfilled', alpha = 0.1,
         bins = bins, label = 'All Galaxies', orientation='horizontal')
for key in keys:
        mgas = df.iloc[np.where(flags[key])[0]].logmgas
        if key == 'defstarform':
            axHisty.hist(mgas, histtype = 'step',orientation='horizontal',
                bins = bins,linewidth = 5, label = labels[key], 
                color = colors[key], alpha = 0.3)
        else:
            axHisty.hist(mgas, histtype = 'step',orientation='horizontal',
                bins = bins,linewidth = 5, label = labels[key], color = colors[key])
        #axHisty.set_xlabel('Gas Mass')
        axHisty.set_xlabel('Number')
'''
axHistx.set_xlim(axScatter.get_xlim())
#axHisty.set_ylim(axScatter.get_ylim())

inputfile = 'C:/Users/mugdhapolimera/github/izi/RESOLVE_SDSS_full.pkl'
full = pd.read_pickle(inputfile)
ra=full.radeg
dec=full.dedeg
if 'fl_insample' in full.keys():
    flinsample = full.fl_insample
else:
    flinsample = np.ones(len(full), dtype = bool)
grpcz = full.grpcz
cz = full.cz
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = full.logmgas
mstars = full.logmstar
mbary = 10**mgas + 10**mstars

inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (((flinsample | (np.log10(mbary) > 9.0)) & infall) | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
full = full[inobssample]
full = full[~np.isnan(full.h_alpha_flux_ext)]
full = full[~np.isnan(full.oiii_5007_flux_ext)]
full = full[~np.isnan(full.nii_6584_flux_ext)]
full = full[~np.isnan(full.nii_6548_flux_ext)]
full = full[~np.isnan(full.h_beta_flux_ext)]
full = full[~np.isnan(full.oi_6300_flux_ext)]
full = full[~np.isnan(full.sii_6717_flux_ext)]
full = full[~np.isnan(full.sii_6731_flux_ext)]
ndx = [x for x in full.NAME if x not in list(flags.galname)]
full = full.loc[ndx]
plt.figure(10)
mbary = np.log10(10**full.logmstar + 10**full.logmgas)
if not he2_flag:
    plt.plot(full.logmh[full.fc == 0], mbary[full.fc == 0], 'ks', alpha = 0.3, 
             markersize = 5) #other group galaxies
    plt.plot(full.logmh[full.fc == 1], mbary[full.fc == 1], 'ks', alpha = 0.3, 
             markersize = 10, label = 'No Em-Lines (Dead)') #central
    plt.legend(loc = 2, numpoints = 1)
percent_fc = {'agntosf': 0, 'composite': 0, 'defagn': 0, 'heiisel': 16,
          'defstarform': 0, 'sftoagn': 0}


for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
        percent_fc[key] = len(np.where(sel.fc ==0)[0])     
        if key == 'defstarform':        
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], markersize = 5, alpha = 0.3) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], markersize = 10, label = labels[key], alpha = 0.3) #central
            
            
        elif key == 'heiisel':
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], markersize = 5, mfc ='none', mew = 2) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], markersize = 10, mfc ='none', mew = 2,label = labels[key]) #central

        else:
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], markersize = 5) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], markersize = 10, label = labels[key]) #central

#        plt.plot(sel.logmh[sel.fc ==0], mbary[sel.fc==0], marker[key], markersize = 5) #other group galaxies
                
        if he2_flag:
            plt.legend(loc = 2, numpoints = 1)
        plt.xlabel('Group Halo Mass')
        plt.ylabel('Galaxy Baryonic Mass')
        
for key in keys:
    print key, ' : ', percent_fc[key]*100.0/percent[key]
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = labels[key])
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + Levesque Grid (no Prior)')
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)
np.savetxt('Filtered_list.txt', flags.galname, fmt ='%s')
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb2.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = labels[key])
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + Levesque Grid (with Prior)')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel_ndx = Z[list(x for x in range(len(Z)) if Z['Name'][x] in sel.NAME)]
        Z_sel = Z_sel_ndx['Z_Estimate']
        plt.plot(sel.logmstar, Z_sel, marker[key], label = labels[key])
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using IZI (Python + GPy) without Prior')
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
plt.ylim(8, 9.2)
plt.plot(m+10,z)

plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel_ndx = Z1[list(x for x in range(len(Z1)) if Z1['Name'][x] in sel.NAME)]
        Z_sel = Z_sel_ndx['Z_Estimate']
        plt.plot(sel.logmstar, Z_sel, marker[key], label = labels[key])
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using IZI (Python + GPy) with Prior')
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
plt.ylim(8, 9.2)
plt.plot(m+10,z)

plt.figure()
bins = np.arange(7.6,9.1,0.1)
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        plt.hist(Z_sel, histtype = 'step', 
                         bins = bins,linewidth = 5, label = labels[key])
        plt.legend()
        plt.xlabel('Metallicity (NebulaBayes + Levesque Grid)')
        plt.ylabel('Number')

plt.figure()
bins = np.arange(7.6,9.1,0.1)
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        Z_sel = Z_nb2.loc[sel.NAME]['Estimate']
        plt.hist(Z_sel, histtype = 'step', 
                         bins = bins,linewidth = 5, label = labels[key])
        plt.legend()
        plt.xlabel('Metallicity (NebulaBayes + Levesque Grid) + Prior')
        plt.ylabel('Number')

plt.figure()
bins = np.arange(7.6,9.1,0.1)
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel_ndx = Z[list(x for x in range(len(Z)) if Z['Name'][x] in sel.NAME)]
        Z_sel = Z_sel_ndx['Z_Estimate']
        plt.hist(Z_sel, histtype = 'step', 
                        bins = bins, linewidth = 5, label = labels[key])
        plt.legend(borderaxespad=0., loc = 2)
        plt.xlabel('Metallicity (IZI - Python + Gpy)')
        plt.ylabel('Number')

plt.figure()
bins = np.arange(7.6,9.1,0.1)
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel_ndx = Z1[list(x for x in range(len(Z1)) if Z1['Name'][x] in sel.NAME)]
        Z_sel = Z_sel_ndx['Z_Estimate']
        plt.hist(Z_sel, histtype = 'step', 
                        bins = bins, linewidth = 5, label = labels[key])
        plt.legend(borderaxespad=0., loc = 2)
        plt.xlabel('Metallicity (IZI - Python + Gpy) without Prior')
        plt.ylabel('Number')

'''plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        #mbary = np.log10(10**sel.logmstar + 10** sel.logmgas)
        bins = np.arange(-6,11,1)
        plt.hist(sel.morph[(sel.morph > -999.0) & (sel.morph < 13.0)], histtype = 'step', 
                         bins= bins,linewidth = 5, label = labels[key])
        plt.legend(loc = 2)
        plt.xlabel('Morphology')
        plt.ylabel('Number')

plt.figure()
bins = np.arange(0,1.5,0.5)
print bins
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        morphel = sel.morphel
        morphel[morphel == 'E'] = 0
        morphel[morphel == 'L'] = 1
        plt.hist(morphel, histtype = 'step', 
                       bins = bins, linewidth = 5, label = labels[key])
        plt.legend(loc = 2)
        plt.xlabel('Morphology')
        plt.ylabel('Number')
'''