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
if he2_flag:
    flags = pd.read_csv('eco+resolve_emlineclass_filter_he2.csv')
else:
    flags = pd.read_csv('eco+resolve_emlineclass_filter.csv')
if full : 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_filter_new.pkl'
    flags = pd.read_csv('eco+resolve_emlineclass_filter.csv')
if resolve: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5.csv'
    flags = pd.read_csv('resolve_emlineclass_full_snr5.csv')
if eco: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_filter_new.pkl'
    flags = pd.read_csv('eco_emlineclass_filter_new.csv')

full_df = pd.read_csv(inputfile, index_col = 'name')
df = full_df.loc[list(flags['galname'])]
if 'NAME' not in df.keys():
    df['NAME'] = df['name']
#results_izi = 'C:/Users/mugdhapolimera/github/izi/results/IZI_Z_filtered.txt'
#full_Z = np.genfromtxt(results_izi, dtype = None, names= True)
#Z = full_Z[list(x for x in range(len(full_Z)) if full_Z['Name'][x] in list(flags.galname))]
        

#os.chdir('C:/Users/mugdhapolimera/github/nebulabayes')
#results_nb = pd.read_csv("eco_izi_SEL/ECO_param_estimates.csv")
#Z_index = (results_nb['Parameter'] == 'LOGZ')
#full_Z_nb = results_nb[Z_index]
#full_Z_nb.index = full_Z_nb['Galaxy Name']
#Z_nb = full_Z_nb.loc[flags.galname]

#results_nb2 = pd.read_csv("results_k13_deredden/RESOLVE_param_estimates.csv")
#Z_index2 = (results_nb2['Parameter'] == 'LOGZ')
#full_Z_nb2 = results_nb2[Z_index2]
#full_Z_nb2.index = full_Z_nb2['Galaxy Name']
#Z_nb2 = full_Z_nb2.loc[flags.galname]

#results_izi = 'C:/Users/mugdhapolimera/github/izi/results/IZI_Z_prior2.txt'
#full_Z1 = np.genfromtxt(results_izi, dtype = None, names= True)
#Z1 = full_Z1[list(x for x in range(len(full_Z1)) if full_Z1['Name'][x] in list(flags.galname))]

keys = ['defstarform', 'defagn', 'composite', 'agntosf', 'sftoagn']
#keys = ['sftoagn1', 'sftoagn2']
#'sftoagn',
if he2_flag:
    keys.append('heiisel')
    
#flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
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
#flags = flags.iloc[list(x for x in range(len(flags)) if flags.galname[x] in list(df.NAME))]
flags.index = flags.galname

plt.figure()
for key in keys:
    sel = df[flags[key]]
    plt.plot(sel.logmstar, sel.modelu_r, marker[key]) #sel.umag-sel.rmag, marker[key])
    plt.xlabel(r'log($M_{*}/M_{\odot}$)')
    plt.ylabel('(u - r)')
    plt.ylim(0,3)

plt.figure()
for key in keys:
    sel = df.iloc[np.where(flags[key])[0]]
    if key == 'defstarform':
        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), sel.meanfsmgr, 
             marker[key], markersize = 10, alpha = 0.3,  mew = 0, 
             color = colors[key], label = labels[key])
    elif key == 'agntosf': 
        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), sel.meanfsmgr, 
             marker[key], markersize = 10, mew = 1, color = colors[key],
             mec = 'y', label = labels[key])
    else:
        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), sel.meanfsmgr, 
             marker[key], markersize = 10, mew = 0, color = colors[key],
             label = labels[key])
    plt.plot(0*np.linspace(0.001,100), np.linspace(0.001,100), 'k-.')
    plt.plot(np.linspace(-2.1,1.6), 1+0*np.linspace(-2.1,1.6), 'k-.')
    plt.text(-0.1, 10**-1.5, r'1:1 G/S Ratio', fontsize=14, color='k', 
             rotation = 'vertical')
    plt.text(-2.0, 1.5, r'Stellar Mass Doubled in last Gyr', 
             fontsize=14, color='k')
    
    plt.xlabel(r'$\rm \log (M_{gas}/M_{stellar})$', size = 22)
    plt.ylabel('Mean FSMGR', size = 22)
    plt.yscale('log')
    yticks = plt.yticks()[0]
    plt.yticks(yticks, np.around(yticks,2))
    plt.ylim(10**-3, 10**2)
    plt.xlim(-2.1,1.6)
    plt.legend(title = 'RESOLVE', loc = 'lower right', fontsize = 14)
    
print percent



'''rect_histy = [left_h, bottom, 0.2, height]
axScatter = plt.axes(rect_scatter)
axHisty = plt.axes(rect_histy, xscale = 'log')
nullfmt = NullFormatter() # no labels
axHisty.yaxis.set_major_formatter(nullfmt)
'''
plt.figure('Physical Properties')
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
        axScatter.plot(np.linspace(7.5,11.5), 
                       np.ones(len(np.linspace(7.5,11.5))), 'k-.')
        axScatter.text(11.0, 0.005, 'RESOLVE', fontsize=14, color='k')
        axScatter.text(10.5, 1.1, r'1:1 G/S Ratio', fontsize=14, color='k')
        
        axScatter.text(9.52, -2.2, r'Gas Richness', fontsize=14, color='k')
        axScatter.text(9.52, -2.4, r'Threshold Mass', fontsize=14, color='k')

        axScatter.plot(9.5*np.ones(len(np.linspace(10**-2.5,10**1.5))), 
                       np.linspace(10**-2.5,10**1.5), 'k-.')
        axScatter.set_ylim(10**-2.5,10**1.5)
        axScatter.set_xlim(7.5,11.5) 
        axScatter.set_yscale("log")
        axScatter.yaxis.set_major_formatter(ScalarFormatter())
        axScatter.legend(loc=2, bbox_to_anchor=(1,1),numpoints = 1, 
                         fontsize = 14)
        #if he2_flag:
            #axScatter.legend(loc=2, bbox_to_anchor=(1,1.15),numpoints = 1)

        axScatter.set_xlabel(r'$\rm \log(M_{stellar})$', fontsize = 22)
        axScatter.set_ylabel(r'$\rm M_{gas}/M_{stellar}$',fontsize = 22)

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
    axHistx.legend(loc=2, bbox_to_anchor=(1,1.15), 
                   fontsize = 14)
    
axHistx.set_ylabel('Number')
axHistx.set_xlim(axScatter.get_xlim())
axHistx.yaxis.set_major_formatter(ScalarFormatter())
        


'''if not he2_flag:
    bins = np.arange(0,1.2,0.05)
    axHisty.hist(np.divide(df.logmgas,df.logmstar), histtype = 'stepfilled', alpha = 0.1,
             bins= bins, linewidth = 5, label = 'All Galaxies',orientation = 'horizontal')
    mgas = np.array(df.iloc[np.where(flags['defstarform'])[0]].logmgas)
    mstars = np.array(df.iloc[np.where(flags['defstarform'])[0]].logmstar)
    axHisty.hist(mgas/mstars, histtype = 'step',bins = bins, alpha = 0.3,
                 linewidth = 5, label = labels['defstarform'], 
                    color = colors['defstarform'],orientation = 'horizontal')
    mgas = df.iloc[np.where((flags['defagn']) | (flags['defseyf']) | 
                    (flags['defliner']))[0]].logmgas
    mstars = df.iloc[np.where((flags['defagn']) | (flags['defseyf']) | 
                    (flags['defliner']))[0]].logmstar
    axHisty.hist(mgas/mstars, histtype = 'step',bins = bins,
                 linewidth = 5, label = 'Traditional AGN', 
                 color = colors['defagn'],orientation = 'horizontal')
    mgas = df.iloc[np.where((flags['sftoagn']))[0]].logmgas
    mstars = df.iloc[np.where((flags['sftoagn']))[0]].logmstar
    axHisty.hist(mgas/mstars, histtype = 'step',bins = bins,
                 linewidth = 5, label = 'Non-Traditional AGN', 
                 color = colors['sftoagn'],orientation = 'horizontal')
    axHisty.legend(loc=2, bbox_to_anchor=(1,1.15))
axHisty.set_xlabel('Number')
axHisty.set_ylim(axScatter.get_ylim())
'''



inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_blend_dext.pkl'
full = pd.read_pickle(inputfile)
ra=full.radeg
dec=full.dedeg
if 'fl_insample' in full.keys():
    flinsample = full.fl_insample
else:
    flinsample = np.ones(len(full), dtype = bool)
grpcz = full.grpcz
#cz = full.cz
#infall = (ra > 22*15.) | (ra < 3*15.)
#inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = full.logmgas
mstars = full.logmstar
mbary = 10**mgas + 10**mstars
full.NAME = full.name
#inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (((flinsample | (np.log10(mbary) > 9.0)) & infall) | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
inobssample = ((grpcz >= 3000.) & (grpcz <= 7000.)) & (np.log10(mbary) > 9.2)
full = full[inobssample]
full = full[~np.isnan(full.h_alpha_flux)]
full = full[~np.isnan(full.oiii_5007_flux)]
full = full[~np.isnan(full.nii_6584_flux)]
full = full[~np.isnan(full.nii_6548_flux)]
full = full[~np.isnan(full.h_beta_flux)]
full = full[~np.isnan(full.oi_6300_flux)]
full = full[~np.isnan(full.sii_6717_flux)]
full = full[~np.isnan(full.sii_6731_flux)]
ndx = [x for x in full.name if x not in list(flags.galname)]
full = full.loc[ndx]
mbary = np.log10(10**full.logmstar + 10**full.logmgas)
'''if not he2_flag:
    plt.plot(full.logmh[full.fc == 0], mbary[full.fc == 0], 'ks', alpha = 0.3, 
             markersize = 5) #other group galaxies
    plt.plot(full.logmh[full.fc == 1], mbary[full.fc == 1], 'ks', alpha = 0.3, 
             markersize = 10, label = 'No Em-Lines (Dead)') #central
    plt.legend(loc = 2, numpoints = 1)
'''
percent_fc = {'agntosf': 0, 'composite': 0, 'defagn': 0, 'heiisel': 16,
          'defstarform': 0, 'sftoagn': 0}

plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
        mbary = np.log10(10**sel.logmstar + 10**sel.logmgas)
        percent_fc[key] = len(np.where(sel.fc ==0)[0])     
        if key == 'defstarform':        
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], 
                     markersize = 5, alpha = 0.3) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], 
                     markersize = 7, label = labels[key], alpha = 0.3) #central
            
        elif key == 'sftoagn':        
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], 
                     markersize = 5, mew = 0) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], 
                     markersize = 7, label = labels[key], mew = 0) #central
            
        elif key == 'heiisel':
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], 
                     markersize = 5, mfc ='none', mew = 2) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key], 
                     markersize = 7, mfc ='none', mew = 2,label = labels[key]) #central

        else:
            plt.plot(sel.logmh[sel.fc == 0], mbary[sel.fc == 0], marker[key], 
                     markersize = 5, mew = 0) #other group galaxies
            plt.plot(sel.logmh[sel.fc == 1], mbary[sel.fc == 1], marker[key],
                     markersize = 7, label = labels[key],mew = 0) #central

#        plt.plot(sel.logmh[sel.fc ==0], mbary[sel.fc==0], marker[key], markersize = 5) #other group galaxies
        plt.plot(11.5*np.ones(len(np.linspace(8.5,12.0))),
                 np.linspace(8.5,12.0), 'k-.', lw = 2)
        plt.plot(12*np.ones(len(np.linspace(8.5,12.0))),
                 np.linspace(8.5,12.0), 'k-.', lw = 2)
        plt.text(11.08,11, 'Gas-Richness') 
        plt.text(11.08,10.9, 'Threshold Mass')        
        plt.text(12.02,11.6, 'Bimodality Mass') 
        if he2_flag:
            plt.legend(loc = 2, numpoints = 1)
        plt.xlabel('Group Halo Mass', fontsize = 14)
        plt.ylabel('Galaxy Baryonic Mass', fontsize = 14)
        plt.xlim(10.7,14.8)
        plt.ylim(9,11.7)
        plt.legend(fontsize = 14, loc = 'upper right')
        #bbox_to_anchor=(1,1), 
        #           loc = 2)
for key in keys:
    print key, ' : ', percent_fc[key]*100.0/percent[key]

plt.figure()
data = [df.iloc[np.where(flags[key])[0]]['logmh'] for key in keys]
hist = plt.hist(data, bins = np.arange(10.5,15,0.5), normed=1, histtype='bar', 
         color=[colors[key] for key in keys], 
         label=[labels[key] for key in keys])    
xaxis = np.linspace(0,np.amax(hist[0])+0.3)
plt.plot(11.5*np.ones(len(xaxis)), xaxis, 'k-.', lw = 2)
plt.plot(12*np.ones(len(xaxis)), xaxis, 'k-.', lw = 2)
plt.text(10.88,1.5, 'Gas-Richness', fontsize = 14) 
plt.text(10.88,1.43, 'Threshold Mass',fontsize = 14)        
plt.text(12.02,1.5, 'Bimodality Mass',fontsize = 14) 
plt.xticks(np.around(hist[1],1))
plt.legend(title = 'ECO', fontsize = 14)
plt.ylim(0,max(xaxis))
plt.ylabel('Normalized Number of Galaxies', fontsize =14 )
plt.xlabel('Group Halo Mass', fontsize = 14)

nb = pd.read_csv('C:/Users/mugdhapolimera/github/nebulabayes/resolve_bpass_full_nicholls+jenkins/RESOLVE_param_estimates.csv')
Z_nb = nb[nb["Parameter"] == 'LOGQ']
Z_nb.index = Z_nb["Galaxy Name"]
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']+8.76
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = labels[key])
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + Levesque Grid (no Prior)')
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)
masses = np.arange(min(mstars), max(mstars), 0.1)
sel = []#(mstars < 0.0)
mass = []
met = []
for i in range(len(masses) - 1):
     
    sel.extend(list(np.where((mstars >= masses[i]) & (mstars < masses[i+1]))[0]))
    #print masses[i], masses[i+1], len(np.where(sel)[0]) #len(np.where((mstars >= masses[i]) & (mstars < masses[i+1]))[0])
    #if len(np.where(sel)[0])> 25:

    if len(sel) > 50:
        Z = [x for x in Z_nb.iloc[sel]['Estimate'] if x != 7.45897]
        #plt.figure()
        #Z_dist = np.hist(Z, bins = 'fd')
        med,sig = norm.fit(Z) #np.median(Z)
        mass.append(masses[i])
        met.append(med)
        sel = []#(mstars < 0.0)
    
    #m = mstars[sel]
    
plt.plot(mass, met, 'k', linewidth = 5, label = 'This Work')

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
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
plt.ylim(8, 9.2)
plt.plot(m+10,z)

plt.figure()
#bins = np.arange(7.6,9.1,0.1)
bins = np.arange(6.9,9.1,0.1)
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

plt.figure()
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

'''
bins = np.arange(0.775,1.2,0.025)
plt.figure()
plt.hist(df.logmgas/df.logmstar, bins = 'fd', alpha = 0.1, 
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
#



#################
plt.figure('Physical Properties')
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
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
                           label = labels[key], markersize = 8, alpha = 0.3, mew = 0)
        
        #elif key == 'sftoagn':
        #    axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
        #                   label = labels[key], markersize = 15, mew = 0)
                
        elif key == 'agntosf':
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
                           label = labels[key], markersize = 10,  mew = 1, 
                            mec = 'y')            
        elif key == 'heiisel':
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key], 
                           label = labels[key], markersize = 10,  mfc = 'None',
                             mew = 2)            
        else:
            axScatter.plot(sel.logmstar,sel.logmgas, marker[key],
                     markersize = 8, mew = 0, label = labels[key])
        axScatter.plot(np.linspace(7.5,11.5), np.linspace(7.5,11.5), 'k-.')
        axScatter.set_ylim(7.5,10.5)        
        if he2_flag:
            axScatter.legend(loc = 4, numpoints = 1)
        axScatter.set_xlabel('Stellar Mass')
        axScatter.set_ylabel('Gas Mass')
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
    axHistx.hist(mstars, histtype = 'step',bins = bins,
                 linewidth = 5, label = 'Traditional AGN', 
                 color = colors['defagn'])
    mstars = df.iloc[np.where((flags['sftoagn']))[0]].logmstar
    axHistx.hist(mstars, histtype = 'step',bins = bins,
                 linewidth = 5, label = 'Non-Traditional AGN', 
                 color = colors['sftoagn1'])
    axHistx.legend(loc=2, bbox_to_anchor=(1,1.15))
axHistx.set_ylabel('Number')
axHistx.set_xlim(axScatter.get_xlim())


lowsfagn = df[flags.sftoagn & (10**df.logmgas/10**df.logmstar < 1.0)]

plt.figure()
for key in keys:
    sel = df.iloc[np.where(flags[key])[0]]
    ssfr = np.log10(sel.sfr_nuv_wise) - sel.logmstar
    if key == 'defstarform':
        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), ssfr, 
             marker[key], markersize = 10, alpha = 0.3,  mew = 0, 
             color = colors[key], label = labels[key])
    elif key == 'agntosf': 
        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), ssfr, 
             marker[key], markersize = 10, mew = 1, color = colors[key],
             mec = 'y', label = labels[key])
    else:
        plt.plot(np.log10(10**sel.logmgas/10**sel.logmstar), ssfr, 
             marker[key], markersize = 10, mew = 0, color = colors[key],
             label = labels[key])
    plt.plot(0*np.linspace(0.001,100), np.linspace(0.001,100), 'k-.')
    #plt.plot(np.linspace(-2.1,1.6), 1+0*np.linspace(-2.1,1.6), 'k-.')
    plt.text(-0.1, 10**-1.5, r'1:1 G/S Ratio', fontsize=14, color='k', 
             rotation = 'vertical')
    plt.text(-2.0, 1.5, r'Stellar Mass Doubled in last Gyr', 
             fontsize=14, color='k')
    
    plt.xlabel(r'$\rm \log (M_{gas}/M_{stellar})$', size = 22)
    plt.ylabel('sSFR', size = 22)
    #plt.yscale('log')
    #yticks = plt.yticks()[0]
    #plt.yticks(yticks, np.around(yticks,2))
    #plt.ylim(10**-3, 10**2)
    #plt.xlim(-2.1,1.6)
    plt.legend(title = 'RESOLVE', loc = 'lower right', fontsize = 14)

plt.figure()
for key in keys:
    sel = df.iloc[np.where(flags[key])[0]]
    ssfr = np.log10(sel.sfr_nuv_wise) - sel.logmstar
    if key == 'defstarform':
        plt.plot(ssfr, sel.meanfsmgr,
             marker[key], markersize = 10, alpha = 0.3,  mew = 0, 
             color = colors[key], label = labels[key])
    elif key == 'agntosf': 
        plt.plot(ssfr, sel.meanfsmgr,
             marker[key], markersize = 10, mew = 1, color = colors[key],
             mec = 'y', label = labels[key])
    else:
        plt.plot(ssfr, sel.meanfsmgr,
             marker[key], markersize = 10, mew = 0, color = colors[key],
             label = labels[key])
    #plt.plot(0*np.linspace(0.001,100), np.linspace(0.001,100), 'k-.')
    plt.plot(np.linspace(-12,-8), 1+0*np.linspace(-12,-8), 'k-.')
    #plt.text(-0.1, 10**-1.5, r'1:1 G/S Ratio', fontsize=14, color='k', 
    #         rotation = 'vertical')
    plt.text(-11.5, 1.5, r'Stellar Mass Doubled in last Gyr', 
             fontsize=14, color='k')
    
    plt.xlabel('sSFR', size = 22)
    plt.ylabel('Mean FSMGR', size = 22)
    plt.yscale('log')
    yticks = plt.yticks()[0]
    plt.yticks(yticks, np.around(yticks,2))
    #plt.ylim(10**-3, 10**2)
    plt.xlim(-12,-8)
    plt.legend(title = 'RESOLVE', loc = 'lower right', fontsize = 14)
