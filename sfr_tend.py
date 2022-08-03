# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 12:30:31 2021

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
import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
pd.set_option('display.max_columns', 10)
from scipy.io.idl import readsav
import scipy
from scipy.stats import kde
from astroML.plotting import scatter_contour
#from scipy.stats import norm
def bisectorslope(fsl,isl):
    # function to compute bisector slope from forward and inverse slopes
    # using formula in Isobe et al (1990)
    bsl1 = (1./(fsl+isl))
    bsl2 = (fsl*isl - 1. + np.sqrt((1.+fsl**2)*(1.+isl**2)))
    bsl = bsl1*bsl2
    return bsl

def rms(modelpts,datapts):
    # function to compute residuals from a model evaluated at same x-values
    resids = modelpts - datapts
    rms = np.sqrt(np.mean(resids**2))
    return rms

def biweight(modelpts,datapts):
    ctune = 9.0
    resids = modelpts - datapts
    med = np.median(resids)
    MAD = np.median(np.abs(resids-med))
    ui=(resids-med)/(ctune*MAD)
    biwt = np.sqrt(len(datapts))*np.sqrt(np.sum(((resids-med)**2*(1-ui**2)**4)))/ \
        np.abs(np.sum(((1-ui**2)*(1-5*ui**2))))
    return biwt
#path = os.path.dirname(os.path.realpath(__file__))
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
he2_flag = 0
full = 1
resolve = 0
eco = 0
sdsscat = 'jhu'
if he2_flag:
    flags = pd.read_csv('eco+resolve_emlineclass_filter_he2.csv')
else:
    flags = pd.read_csv('eco+resolve_emlineclass_filter.csv')
if full : 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_snr5_dext_jhu.csv'
    flags = pd.read_csv(r'ECO\SEL\eco+resolve_emlineclass_dext_snr5_jhu.csv')
#    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_snr5_dext_'+sdsscat+'.csv'
#    flags = pd.read_csv(r'ECO\SEL\eco+resolve_emlineclass_dext_snr5_jhu.csv')
if resolve: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5_dext_'+sdsscat+'.csv'
    flags = pd.read_csv('resolve_emlineclass_dext_snr5_'+sdsscat+'.csv')
    conf = pd.read_csv('RESOLVE_snr5_master_conf_new.csv')
if eco: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_filter_new.pkl'
    flags = pd.read_csv('eco_emlineclass_filter_new.csv')

full_df = pd.read_csv(inputfile)
full_df.index = full_df.name
seldf = full_df.loc[list(flags['galname'])]
if 'NAME' not in seldf.keys():
    seldf['NAME'] = seldf['name']

keys = ['defstarform', 'defagn', 'composite', 'agntosf', 'sftoagn']
if he2_flag:
    keys.append('heiisel')

    
marker = {'agntosf': 'g^', 'ambigagn': 'ms', 'composite': 'ms', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.','sftoagn': 'bs', 'sftoagn1': 's', 'sftoagn2': 'm*'}

colors = {'agntosf': 'g', 'ambigagn': 'm', 'composite': 'm', 'defagn': 'r', 
          'defliner': 'y', 'defseyf': 'c', 'heiisel': 'k',
          'defstarform': 'gray', 'sftoagn': 'b', 'sftoagn2': 'b'}

labels = {'agntosf': 'Low-SII AGN', 'ambigagn': 'Ambiguous AGN', 
          'composite': 'Composite', 'defagn': 'Definite AGN', 
          'defliner': 'LINER', 'defseyf': 'Seyfert', 
          'heiisel': 'HeII-Selected AGN', 'defstarform': 'Definite SF', 
          'sftoagn': 'SFing-AGN', 'sftoagn2' : 'MP-AGN2'}

percent = {'agntosf': 0, 'composite': 0, 'defagn': 0, 'heiisel': 57,
          'defstarform': 0, 'sftoagn': 0}

for key in keys:
    percent[key] = len(np.where(flags[key])[0])

seldf = seldf[seldf.logmstar < 9.5]
flags.index = flags.galname



def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 

if resolve:
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_blend_dext_new.csv'
    df = pd.read_csv(inputfile)
    print len(df)
    ra=df.radeg
    dec=df.dedeg
    flinsample = df.fl_insample
    grpcz = df.grpcz
    cz = df.cz
    infall = (ra > 22*15.) | (ra < 3*15.)
    inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
    mgas = df.logmgas
    mstars = df.logmstar
    mbary = 10**mgas + 10**mstars
    if resolve:
        inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
    (((flinsample | (np.log10(mbary) > 9.2)) & infall) | \
            ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
    df = df[inobssample]
if full:
    inputfile = 'ECO+RESOLVE_inobssample.csv'
    df = pd.read_csv(inputfile)
    df = df[df.logmstar < 9.5]

ecogals = df.resname == 'notinresolve'
#for key in keys:
#    sel = np.array(df.loc[flags[key]].name)
#    ax.plot(df.logmstar.loc[flags.index.values[sel]], 
#        u_r.loc[flags.index.values[sel]], marker = marker[key], mec = 'k', markersize = 10)#, label = 'SFing-AGN'),
#ax.plot(df.logmstar.loc[target], u_r.loc[target],'r*', 
#         markersize = 10, label = 'Targets for this proposal')
#for i, txt in enumerate(target):
#    plt.annotate(txt, (df.logmstar.loc[target][i], u_r.loc[target][i]))
#plt.plot(df.logmstar.loc[ifu], u_r.loc[ifu],'*', c = 'magenta', 
#         markersize = 10, label = 'New IFU?')
#plt.plot(df.logmstar.loc[oisel], u_r.loc[oisel],
#         'rs', mfc = 'none', mec = 'r', markersize = 15, label = 'Strong [OI]')
#ax.legend(handles=[sfagnplt,targetsplt,samiplt], 
#          labels = ['SFing-AGN', 'Targets for this proposal','SFing-AGN in SAMI'],
#          fontsize = 12, loc = 'upper left')
#plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 15)
#plt.ylabel('u-r', fontsize = 15)
#plt.savefig('gemini_u_r_nov_new.jpeg', quality = 95)
#for i, txt in enumerate(nov_target):
#    plt.annotate(txt, (df.logmstar.loc[nov_target][i], u_r.loc[nov_target][i]))
#for i, txt in enumerate(sami):
#    plt.annotate(txt, (df.logmstar.loc[sami][i], u_r.loc[sami][i]))
#for i, txt in enumerate(ifu):
#    plt.annotate(txt, (df.logmstar.loc[ifu][i], u_r.loc[ifu][i]))

#FSMGR vs SSFR
ssfr_st = 10**(np.log10(df.sfr_nuv) - df.logmstar)
#ssfr_st[ecogals] = 10**(np.log10(df.sfr_nuv[ecogals]) - df.logmstar[ecogals])

ssfr_lt = 1/(1+(1/(df.meanfsmgr)))
fsmgr_st = np.log10(100*(10**6)*(ssfr_st)/(1-ssfr_st))
fsmgr_lt = np.log10(df.meanfsmgr)

xmin = np.nanmin(fsmgr_st)
xmax = np.nanmax(fsmgr_st)+0.5
ymin = np.nanmin(fsmgr_lt)
ymax = np.nanmax(fsmgr_lt)+0.5
#nbins = 500
#
fig,ax = plt.subplots(figsize=(16,8))
#sfrgrid = np.column_stack((fsmgr_st,fsmgr_lt))
#xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
#k2 = kde.gaussian_kde(sfrgrid.T)
#sfrgrid_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
#ax.pcolormesh(xgrid, ygrid, sfrgrid_z.reshape(xgrid.shape), 
#               cmap='bone_r')
#ax.contour(xgrid, ygrid, sfrgrid_z.reshape(xgrid.shape), 3, 
#            cmap='summer', corner_mask = True, extend = 'both')
X,Y,Z = density_estimation(fsmgr_st,fsmgr_lt)
#ax.pcolormesh(X, Y, Z.reshape(X.shape), 
#               cmap='bone_r')
#ax.contour(X, Y, Z.reshape(X.shape), 3, 
#            cmap='summer', corner_mask = True, extend = 'both')

#plt.imshow(np.rot90(Z), cmap='bone_r',                                                    
#          extent=[xmin, xmax, ymin, ymax], interpolation='gaussian', aspect = 'auto')
#levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
levels = 100*np.array([0.001, 0.021, 0.14, 0.34,0.5,0.68])
#pts, conts = scatter_contour(fsmgr_st, fsmgr_lt, levels = levels,
#                threshold=10, log_counts=False, ax = ax, 
#                histogram2d_args=dict(bins=100),
#                plot_args=dict(marker=',', linestyle='none', color='black'),
#                contour_args=dict(cmap=plt.cm.bone_r))
#plt.clim(0,1.8)
plt.plot(fsmgr_st, fsmgr_lt, '.', color = 'lightgray', zorder = 0, label = 'Non-SEL galaxies')


#cont = plt.contour(X, Y, Z, levels = levels, ,cmap='summer')
inversefit = np.polyfit(fsmgr_lt, fsmgr_st,1)
invslope = 1/inversefit[0]
invconst = -1*inversefit[1]/inversefit[0]
xaxis = np.arange(xmin, xmax,0.01)
#plt.plot(xaxis, invslope*xaxis+invconst,'cyan', zorder = 5, lw = 3)

fit = np.polyfit(fsmgr_st, fsmgr_lt,1)
slope = fit[0]
const = fit[1]
xaxis = np.arange(xmin, xmax,0.01)
#plt.plot(xaxis, slope*xaxis+const,'purple')

slopebis = bisectorslope(slope,invslope)
intbis = np.mean(fsmgr_lt) - slopebis*np.mean(fsmgr_st)
#plt.plot(xaxis,slopebis*xaxis+intbis,'m--')

#find median ST fsgmr for every value of long term fsmgr
dy = 0.5
bins = np.arange(ymin,ymax,dy)
bin_center = (bins[0:-1]+bins[1:])/2
fsmgr_st_50p = np.zeros(len(bins)-1)
fsmgr_st_16p = np.zeros(len(bins)-1)
fsmgr_st_84p = np.zeros(len(bins)-1)
fsmgr_st_med = np.zeros(len(bins)-1)
for i in range(len(bins)-1):
    subset = fsmgr_st[(np.logical_and(fsmgr_lt>=bins[i], fsmgr_lt<=bins[i+1]))]
    print bins[i], len(subset)
    if len(subset>0):
        fsmgr_st_med[i] = np.nanmedian(subset)
        fsmgr_st_50p[i] = np.percentile(subset,50)
        fsmgr_st_16p[i] = np.percentile(subset,16)
        fsmgr_st_84p[i] = np.percentile(subset,84)
    else:
        fsmgr_st_med[i] = np.nan
        fsmgr_st_50p[i] = np.nan
        fsmgr_st_16p[i] = np.nan
        fsmgr_st_84p[i] = np.nan

plt.plot(fsmgr_st_50p,bin_center,'cyan', zorder = 5, lw = 2)
plt.plot(fsmgr_st_16p,bin_center,'cyan',zorder = 5, lw = 2)
plt.plot(fsmgr_st_84p,bin_center,'cyan', zorder = 5, lw = 2)

selecogals = seldf.resname == 'notinresolve'

ssfr_st = 10**(np.log10(seldf.sfr_nuv) - seldf.logmstar)
#ssfr_st[selecogals] = 10**(np.log10(seldf.sfr_nuv[selecogals]) - seldf.logmstar[selecogals])
ssfr_lt = 1/(1+(1/(seldf.meanfsmgr)))
fsmgr_st = np.log10(100*(10**6)*(ssfr_st)/(1-ssfr_st))
fsmgr_lt = np.log10(seldf.meanfsmgr)
#key = 'defstarform'
#sel = np.array(seldf.loc[flags[key]].name)
#plt.plot(fsmgr_st.loc[sel], fsmgr_lt.loc[sel],
#            'k.', label = labels[key])
#
#key = 'sftoagn'
#sel = np.array(seldf.loc[flags[key]].name)
#plt.plot(fsmgr_st.loc[sel], fsmgr_lt.loc[sel],
#            marker[key])

for key in keys:
    sel = np.array(seldf.loc[flags[key]].name)
    if key == 'sftoagn':
        zo = 4
    else:
        zo = 2
    if key == 'defstarform':
        plt.plot(fsmgr_st.loc[sel], fsmgr_lt.loc[sel],
            'k.', label = labels[key], zorder = zo)
    else:    
        plt.plot(fsmgr_st.loc[sel], fsmgr_lt.loc[sel],
            marker[key], ms = 10, label = labels[key], zorder = zo)
#xaxis = np.arange(np.min(fsmgr_st)-0.05, np.max(fsmgr_st)+0.05,0.01)
yaxis = np.zeros(len(xaxis))
plt.plot(xaxis, yaxis, 'k--', lw = 3)
#plt.xlim(np.min(fsmgr_st), np.max(fsmgr_st))
#plt.ylim(np.min(fsmgr_lt), np.max(fsmgr_lt))
plt.xlim(-4.1, -0.5)
plt.ylim(-3.5,1.5)
plt.text(-4.0, 0.15, r'Stellar Mass Doubled in last Gyr', 
             fontsize=14, color='k')
#plt.yscale('log')
#plt.xscale('log')
plt.legend(fontsize = 15, loc = 'lower right')
#plt.plot(fsmgr_st.loc[jwst], fsmgr_lt.loc[jwst],
#            '*',color = 'lime', ms = 10, label = labels[key])
plt.xlabel(r'Short Term FSMGR $\left(\frac{M_*(<100 Myr)}{M_*(>100 Myr)}\right)$', 
           fontsize = 22)
plt.ylabel(r'Long Term FSMGR $\left(\frac{M_*(<1 Gyr)}{M_*(>1 Gyr)}\right)$', 
           fontsize = 22)
dwarf = seldf.name[seldf.logmstar < 9.5]
#for i, txt in enumerate(df.name):
#    plt.annotate(txt, (fsmgr_st.loc[df.name[i]], 
#                             fsmgr_lt.loc[df.name[i]]))
#for i, txt in enumerate(dwarf):
#    plt.annotate(txt, (fsmgr_st.loc[dwarf[i]], 
#                             fsmgr_lt.loc[dwarf[i]]))


#df = df[df.logmstar < 9.5]
#seldf = seldf[seldf.logmstar < 9.5]

dfinobs = df[np.log10(10**df.logmstar+10**df.logmgas) > 9.2]
dfinobs = dfinobs[dfinobs.logmstar > 7]
plt.figure()
plt.plot(dfinobs.logmstar, dfinobs.sfr_nuv,
            '.', color = 'gray', label = labels[key])
for key in keys:
    sel = np.array(seldf.loc[flags[key]].name)
    if key == 'defstarform':
        plt.plot(seldf.logmstar.loc[sel], seldf.sfr_nuv.loc[sel],
            'k.', label = labels[key])
    else:    
        plt.plot(seldf.logmstar.loc[sel], seldf.sfr_nuv.loc[sel],
            marker[key], ms = 10, label = labels[key])
plt.yscale('log')
plt.legend(fontsize = 15, loc = 'lower right')
plt.xlabel(r'$M_*$', 
           fontsize = 22)
plt.ylabel(r'SFR $(M_\odot/yr)$', 
           fontsize = 22)
#find median ST fsgmr for every value of long term fsmgr

ymin = 7.75
ymax = 9.75
dy = 0.25
bins = np.arange(ymin,ymax,dy)
bin_center = (bins[0:-1]+bins[1:])/2
sfr_50p = np.zeros(len(bins)-1)
sfr_16p = np.zeros(len(bins)-1)
sfr_84p = np.zeros(len(bins)-1)
sfr_med = np.zeros(len(bins)-1)
for i in range(len(bins)-1):
    subset = dfinobs.sfr_nuv[(np.logical_and(dfinobs.logmstar>=bins[i], dfinobs.logmstar<=bins[i+1]))]
    print bins[i], len(subset)
    if len(subset>0):
        sfr_med[i] = np.nanmedian(subset)
        sfr_50p[i] = np.percentile(subset,50)
        sfr_16p[i] = np.percentile(subset,16)
        sfr_84p[i] = np.percentile(subset,84)
    else:
        sfr_med[i] = np.nan
        sfr_50p[i] = np.nan
        sfr_16p[i] = np.nan
        sfr_84p[i] = np.nan

plt.plot(bin_center,sfr_50p,'cyan', zorder = 5, lw = 2)
plt.plot(bin_center,sfr_16p,'cyan',zorder = 5, lw = 2)
plt.plot(bin_center,sfr_84p,'cyan', zorder = 5, lw = 2)








