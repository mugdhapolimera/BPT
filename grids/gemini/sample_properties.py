# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:30:00 2019

@author: mugdhapolimera
"""
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import matplotlib as mpl
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')

df = pd.read_csv('RESOLVE_full_snr5.csv')
df.index = df.name
flags = pd.read_csv('resolve_emlineclass_full_snr5.csv')
flags.index = flags.galname
df2 = pd.read_csv('RESOLVE_full_snr5_port.csv')
df2.index = df2.name
nsa = pd.read_csv('NSA_RESOLVE_smcdext.csv')
nsa.index = nsa.resname
#target = ['rf0547', 'rf0338', 'rf0306', 'rf0477']
target = targetnames
nov_target = ['rs0323','rs0933','rs0472', 'rs1055', 'rs0124'] 
#rs1172 - low met
#'rs0930', , 'rs1105', 'rs1021' - upper locus of met
sami = ['rs0010','rs0775']
#'rs0404-not in sample','rs0167',
df.loc['rf0306','nii_6584_flux'] = nsa.loc['rf0306'].nii_6584_flux
df.loc['rf0306','h_alpha_flux'] = nsa.loc['rf0306'].h_alpha_flux
df.loc['rs0323','nii_6584_flux'] = nsa.loc['rs0323'].nii_6584_flux
df.loc['rs0323','h_alpha_flux'] = nsa.loc['rs0323'].h_alpha_flux
flags.loc['rf0547'] = flags.loc['rf0338']
flags.loc['rf0306'] = flags.loc['rf0338']
flags.loc['rf0702'] = flags.loc['rf0338']
flags.loc['rf0105'] = flags.loc['rf0338']
nonjhusfagn = [x for x in targetnames if x not in list(flags.index.values[flags.sftoagn])]
flags.sftoagn.loc[nonjhusfagn] = True
flags.defstarform.loc[nonjhusfagn] = False
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
Z_pp04[bad] = -99
plt.figure()
ax1 = plt.subplot()
plt.plot(df.logmstar,Z_pp04 ,'k.', alpha = 0.3, label = 'All galaxies')
plt.plot(df.logmstar.loc[flags.sftoagn], Z_pp04.loc[flags.sftoagn],
         'bs', mec = 'k', markersize = 10, label = 'SFing-AGN')
plt.plot(df.logmstar.loc[nov_target], Z_pp04.loc[nov_target],'r*', 
         markersize = 10, label = 'Targets for this proposal')
plt.plot(df.logmstar.loc[sami], Z_pp04.loc[sami],'*', c = 'lime', 
         markersize = 10, label = 'SFing-AGN in SAMI')
#plt.plot(df.logmstar.loc[target], Z_pp04.loc[target],'r*', 
#         markersize = 10, label = 'Targets for this proposal')
#for i, txt in enumerate(target):
#    ax1.annotate(txt, (df.logmstar.loc[target][i], Z_pp04.loc[target][i]))
oisel = df.name.loc[flags.sftoagn][(df.oi_6300_flux.loc[flags.sftoagn]/df.oi_6300_flux_err.loc[flags.sftoagn] > 10)]
#plt.plot(df.logmstar.loc[oisel], Z_pp04.loc[oisel],
#         'rs', mfc = 'none', mec = 'r', markersize = 15, label = 'Strong [OI]')
yaxis = np.linspace(np.min(Z_pp04[~bad])-0.1,np.max(Z_pp04[~bad])+0.1)
plt.plot(9.5*np.ones(len(yaxis)),yaxis, 'k-.', linewidth = 3)
#plt.ylim(min(yaxis),max(yaxis))
plt.ylim(7.8,9.0)
plt.xlim(7.5,11.0)
plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 15)
plt.ylabel(r'Z (12 + log(O/H))', fontsize = 15)
plt.legend(fontsize = 12)
ax2 = ax1.twinx()
#ax2.set_xticks([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1.0])
yticks = np.array([7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0])
ax2.set_yticks(np.linspace(0,1,len(yticks)))
ax1.set_yticks(np.arange(7.8,9.0,0.2))
float_formatter = lambda x: "%.2f" % x
#xticks = np.array([7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2])
N2 = 1.754*yticks - 15.614
N2_label = ["%.2f" % z for z in N2]
ax2.set_yticklabels(N2_label)
ax2.set_ylabel(r'[NII]/H$\alpha$', fontsize = 15)
#for i, txt in enumerate(nov_target):
#    ax1.annotate(txt, (df.logmstar.loc[nov_target][i], Z_pp04.loc[nov_target][i]))
#for i, txt in enumerate(sami):
#    ax1.annotate(txt, (df.logmstar.loc[sami][i], Z_pp04.loc[sami][i]))


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

ssfr_st = np.log10(df.sfr_nuv_wise) - df.logmstar
ssfr_lt = 1/(1+(1/(df.meanfsmgr)))
g_s = np.log10(10**df.logmgas/10**df.logmstar)
cmap = plt.get_cmap('gray',6)
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
#for i in range(1,5):
#    cmaplist[i] = cmaplist[i-1]
cmaplist[5] = (0.9,0.9,0.9,1.0)
# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)
#map = truncate_colormap(old_cmap,minval = 0.2,maxval = 0.8,n=6)
fig,ax = plt.subplots()
cax = ax.scatter(10**ssfr_st, ssfr_lt, marker = 'o', s=8, 
                 c = g_s, cmap = cmap, label = 'All galaxies')
fig.colorbar(cax,extend = 'min')
sfagn = [x for x in flags.sftoagn.index.values \
         if flags.loc[x].sftoagn & (x not in target)]   
ax.plot(10**ssfr_st.loc[flags.sftoagn], ssfr_lt.loc[flags.sftoagn],
        's', ms = 12,mec = 'b', mfc= 'none',
        mew = 1, label = 'SFing-AGN')
ax.plot(10**ssfr_st.loc['rs0010'], ssfr_lt.loc['rs0010'],'*', 
         ms = 12,mec = 'lime',mfc= 'none',mew = 2, 
         label = 'SFing-AGN in SAMI')
ax.plot(10**(ssfr_st.loc['rs0775']-0.03), ssfr_lt.loc['rs0775']+0.03,'*', 
         ms = 12,mec = 'lime',mfc= 'none',mew = 2)
ax.plot(10**ssfr_st.loc[nov_target], ssfr_lt.loc[nov_target],'*', 
         ms = 12,mec = 'r',mfc= 'none',mew = 2, 
         label = 'Targets for this proposal')
#for i, txt in enumerate(nov_target):
#    ax.annotate(txt, (10**ssfr_st.loc[nov_target][i], \
#                      ssfr_lt.loc[nov_target][i]))
#for i, txt in enumerate(sami):
#    ax.annotate(txt, (10**ssfr_st.loc[sami][i], \
#                      ssfr_lt.loc[sami][i]))
#ax.plot(10**ssfr_st.loc[target], ssfr_lt.loc[target],'*', 
#         ms = 10,mec = 'k',mfc= 'none',mew = 2, 
#         label = 'Targets for this proposal')
#for i, txt in enumerate(target):
#    ax.annotate(txt, (10**ssfr_st.loc[target][i], \
#                      ssfr_lt.loc[target][i]))
galname = 'rs0323'
#ax.plot(10**ssfr_st.loc[galname], ssfr_lt.loc[galname],
#         'rs', mfc = 'none', mec = 'r',mew = 5, markersize = 15, label = 'rs0010')
#plt.plot(ssfr_lt/10**df.logmstar/1e9,df.meanfsmgr,'.')
ax.set_xlim(np.min(10**ssfr_st), np.max(10**ssfr_st))
ax.set_ylim(np.min(ssfr_lt)-0.00001, np.max(ssfr_lt)+1.1)
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 12)
ax.set_xlabel('Short Term sSFR', fontsize = 15)
ax.set_ylabel('Long Term sSFR (x10$^9$)', fontsize = 15)

#u_r = df['modelu_rcorr']#df['umag'] - df['rmag']
#fig,ax = plt.subplots()
#ax.plot(df.logmstar,u_r,'k.', alpha = 0.3, label = 'All galaxies')
#ax.plot(df.logmstar.loc[flags.sftoagn], u_r.loc[flags.sftoagn],
#         'bs', mec = 'k', markersize = 10, label = 'SFing-AGN')
##ax.plot(df.logmstar.loc[target], u_r.loc[target],'r*', 
##         markersize = 10, label = 'Targets for this proposal')
#ax.plot(df.logmstar.loc[nov_target], u_r.loc[nov_target],'r*', 
#         markersize = 10, label = 'Targets for this proposal')
##plt.plot(df.logmstar.loc[oisel], u_r.loc[oisel],
##         'rs', mfc = 'none', mec = 'r', markersize = 15, label = 'Strong [OI]')
#ax.legend(fontsize = 12)
#ax.set_xlabel('log(M*)', fontsize = 15)
#ax.set_ylabel('u-r', fontsize = 15)
#for i, txt in enumerate(nov_target):
#    ax.annotate(txt, (df.logmstar.loc[nov_target][i], u_r.loc[nov_target][i]))

#plt.figure()
#plt.plot(df.logmstar,ssfr_st ,'k.', alpha = 0.3, label = 'All galaxies')
#plt.plot(df.logmstar.loc[flags.sftoagn], ssfr_st.loc[flags.sftoagn],
#         'bs', mec = 'k', markersize = 10, label = 'SFing-AGN')
#plt.plot(df.logmstar.loc[target], ssfr_st.loc[target],'r*', 
#         markersize = 10, label = 'Targets for this proposal')
#plt.plot(df.logmstar.loc[oisel], ssfr_st.loc[oisel],
#         'rs', mfc = 'none', mec = 'r', markersize = 15, label = 'Strong [OI]')
#plt.legend(fontsize = 12)
#plt.xlabel('log(M*)', fontsize = 15)
#plt.ylabel('Short Term sSFR', fontsize = 15)
#
#plt.figure()
#plt.plot(g_s,ssfr_st ,'k.', alpha = 0.3, label = 'All galaxies')
#plt.plot(g_s.loc[flags.sftoagn], ssfr_st.loc[flags.sftoagn],
#         'bs', mec = 'k', markersize = 10, label = 'SFing-AGN')
#plt.plot(g_s.loc[target], ssfr_st.loc[target],'r*', 
#         markersize = 10, label = 'Targets for this proposal')
#plt.plot(g_s.loc[oisel], ssfr_st.loc[oisel],
#         'rs', mfc = 'none', mec = 'r', markersize = 15, label = 'Strong [OI]')
#plt.legend(fontsize = 12)
#ax.set_xlabel('log(G/S)', fontsize = 15)
#ax.set_ylabel('Short Term sSFR', fontsize = 15)

import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
pd.set_option('display.max_columns', 10)
import numpy as np
import matplotlib.mlab as mlab
from scipy.io.idl import readsav
import pdb
from matplotlib.patches import Ellipse
from matplotlib.ticker import ScalarFormatter
import time
#import vapeplot
#from astroML.plotting import scatter_contour
import scipy
from collections import OrderedDict
from matplotlib import cm
import math
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
import matplotlib.image as mpimg

def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 

xmin = 7.5
xmax = 11.5
ymin = 0
ymax = 3
df = resolve #pd.read_pickle(inputfile)
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
inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
(((flinsample | (np.log10(mbary) > 9.0)) & infall) | \
        ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
df = df[inobssample]
u_r = df['modelu_rcorr']

X,Y,Z = density_estimation(df.logmstar,u_r)
f20 = plt.figure()#(figsize=(8,8))
plt.imshow(np.rot90(Z), cmap='bone_r',                                                    
          extent=[xmin, xmax, ymin, ymax], interpolation='gaussian')

plt.plot(df.logmstar.loc[flags.index.values[flags.sftoagn]], 
        u_r.loc[flags.index.values[flags.sftoagn]],
         'bs', mec = 'k', markersize = 10, label = 'SFing-AGN')
#ax.plot(df.logmstar.loc[target], u_r.loc[target],'r*', 
#         markersize = 10, label = 'Targets for this proposal')
#for i, txt in enumerate(target):
#    plt.annotate(txt, (df.logmstar.loc[target][i], u_r.loc[target][i]))
plt.plot(df.logmstar.loc[nov_target], u_r.loc[nov_target],'r*', 
         markersize = 10, label = 'Targets for this proposal')
plt.plot(df.logmstar.loc[sami], u_r.loc[sami],'*', c = 'lime', 
         markersize = 10, label = 'SFing-AGN in SAMI')
#plt.plot(df.logmstar.loc[oisel], u_r.loc[oisel],
#         'rs', mfc = 'none', mec = 'r', markersize = 15, label = 'Strong [OI]')
plt.legend(fontsize = 12, loc = 'upper left')
plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 15)
plt.ylabel('u-r', fontsize = 15)
plt.clim(0,1.8)
plt.contour(X, Y, Z, cmap='summer')
plt.savefig('gemini_u_r_nov_new.jpeg', quality = 95)
#for i, txt in enumerate(nov_target):
#    plt.annotate(txt, (df.logmstar.loc[nov_target][i], u_r.loc[nov_target][i]))
#for i, txt in enumerate(sami):
#    plt.annotate(txt, (df.logmstar.loc[sami][i], u_r.loc[sami][i]))



