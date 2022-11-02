# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 10:12:42 2022

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from scipy.stats import kde
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colors as mpcolors
import scipy
from scipy.io.idl import readsav
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'axes.linewidth': 2})
matplotlib.rcParams.update({'lines.linewidth': 2})

def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 
def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
  	new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap


plot = 0
#os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
#sdsscat = 'jhu'
#s06file = 'eco+resolve_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
#bptfile = 'eco+resolve_emlineclass_dext_snr5_'+sdsscat+'.csv'
#midirfile = 'mid_ir/ECO+RESOLVE_WISE_good_syserr.csv'
survey = 'ECO+RESOLVE'
#survey = 'RESOLVE'
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
sdsscat = 'nsa'
s06file = survey.lower()+'_s06emlineclass_dext_hasnr5_'+sdsscat+'.csv'
optfile = survey.lower()+'_emlineclass_dext_snr5_'+sdsscat+'.csv'
bptfile = survey.lower()+'_emlineclass_bpt1snr5_'+sdsscat+'.csv'
irrandfile = 'mid_ir/RESOLVE_WISE_good_randerr.csv'
irsysfile = 'mid_ir/RESOLVE_WISE_good_syserr.csv'

s06 = pd.read_csv(s06file)
s06.index = s06.galname
bpt = pd.read_csv(bptfile)
bpt.index = bpt.galname
opt = pd.read_csv(optfile)
opt.index = opt.galname
irrand = pd.read_csv(irrandfile)
irrand.index = irrand.name
irsys = pd.read_csv(irsysfile)
irsys.index = irsys.name

resdfname = survey+"_barysample.csv"
resdf = pd.read_csv(resdfname)
resdf.index = resdf.name

bpt = bpt.loc[set(bpt.galname) & set(resdf.name)]
opt = opt.loc[set(opt.galname) & set(resdf.name)]
s06 = s06.loc[set(s06.galname) & set(resdf.name)]



print(survey+' parent survey: ', len(resdf))
print(survey+' BPT-only sample: ', len(bpt))
print(survey+' S06-only sample: ', len(s06))
print(survey+' SEL sample: ', len(opt))
print(survey+' IR rand-err sample: ', len(irrand))
print(survey+' IR sys-err sample: ', len(irsys))

print(survey+' parent dwarf survey: ', np.sum(resdf.logmstar < 9.5))
print(survey+' BPT-only dwarf sample: ', np.sum(resdf.logmstar.loc[bpt.galname] < 9.5))
print(survey+' S06-only dwarf sample: ', np.sum(resdf.logmstar.loc[s06.galname] < 9.5))
print(survey+' SEL dwarf sample: ', np.sum(resdf.logmstar.loc[opt.galname] < 9.5))
print(survey+' IR rand-err dwarf sample: ', np.sum(resdf.logmstar.loc[irrand.name] < 9.5))
print(survey+' IR sys-err dwarf sample: ', np.sum(resdf.logmstar.loc[irsys.name] < 9.5))


optsf = set(opt.galname[opt.defstarform])
optsfagn = set(opt.galname[opt.sftoagn])
bptsf = set(bpt.galname[bpt.defstarform])
opt_tradagn = set(opt.galname[~(opt.defstarform | opt.sftoagn)])
bpt_tradagn = set(bpt.galname[~bpt.defstarform])

opt_bpt = set(np.unique(list(opt.galname) + list(bpt.galname)))
opt_bpt_sf = (optsf| bptsf) ^ optsfagn
opt_bpt_tradagn = opt_tradagn | bpt_tradagn
opt_bpt_agn = opt_bpt - opt_bpt_sf

print(survey+' updated SEL sample: ', len(opt_bpt))
print(survey+' updated SEL SF: ', len(opt_bpt_sf))
print(survey+' updated SEL Traditional AGN: ', len(opt_bpt_tradagn))

dwarf = set(resdf[resdf.logmstar < 9.5].name)
opt_bpt_dwarf = opt_bpt & dwarf
opt_bpt_sf_dwarf = opt_bpt_sf & dwarf
opt_sfagn_dwarf = optsfagn & dwarf
opt_bpt_tradagn_dwarf = opt_bpt_tradagn & dwarf

print(survey+' updated SEL dwarf sample: ', len(opt_bpt_dwarf))
print(survey+' updated SEL dwarf SF: ', len(opt_bpt_sf_dwarf))
print(survey+' updated SEL dwarf SF-AGN: ', len(opt_sfagn_dwarf))
print(survey+' updated SEL dwarf Traditional AGN: ', len(opt_bpt_tradagn_dwarf))

s06gals = set(s06.galname)
s06sf = set(s06.galname[s06.defstarform])
s06trad = set(s06.galname[s06.defagn])
s06bonus = set(s06.galname[s06.composite])

print(survey+' S06 sample: ', len(s06))
print(survey+' S06 SF: ', np.sum(s06.defstarform))
print(survey+' S06 Traditional AGN: ', np.sum(s06.defagn))
print(survey+' S06 Bonus AGN: ', np.sum(s06.composite))

s06_dwarf = s06gals & dwarf
s06sf = s06sf & dwarf
s06trad_dwarf = s06trad & dwarf
s06bonus_dwarf = s06bonus & dwarf

print(survey+' S06 dwarf sample: ', np.sum(resdf.loc[s06.galname].logmstar < 9.5))
print(survey+' S06 dwarf SF: ', np.sum(resdf.loc[s06.galname[s06.defstarform]].logmstar < 9.5))
print(survey+' S06 dwarf bonus AGN: ', np.sum(resdf.loc[s06.galname[s06.composite]].logmstar < 9.5))
print(survey+' S06 dwarf Traditional AGN: ', np.sum(resdf.loc[s06.galname[s06.defagn]].logmstar < 9.5))

optdwarfagn = s06bonus_dwarf | opt_bpt_tradagn_dwarf | opt_sfagn_dwarf
optagn = s06bonus | opt_bpt_tradagn | optsfagn

#irgals = set(irsys.name)
#iragn = set(irsys.name[irrand.agnflag]) 
irgals = set(irrand.name)
iragn = set(irrand.name[irrand.agnflag]) 

irdwarf = irgals & dwarf
irdwarfagn = iragn & dwarf 

print(survey+' IR rand-err sample: ', len(irrand))
print(survey+' IR rand-err AGN: ', np.sum(irrand.agnflag))
print(survey+' IR rand-err dwarf sample: ', np.sum(resdf.logmstar.loc[irrand.name] < 9.5))
print(survey+' IR rand-err dwarf AGN: ', np.sum(resdf.logmstar.loc[irrand.name[irrand.agnflag]] < 9.5))

print(survey+' IR sys-err sample: ', len(irsys))
print(survey+' IR sys-err AGN: ', np.sum(irsys.agnflag))
print(survey+' IR sys-err dwarf sample: ', np.sum(resdf.logmstar.loc[irsys.name] < 9.5))
print(survey+' IR sys-err dwarf AGN: ', np.sum(resdf.logmstar.loc[irsys.name[irrand.agnflag]] < 9.5))

dwarfagn = irdwarfagn | s06bonus_dwarf | opt_bpt_tradagn_dwarf | opt_sfagn_dwarf
agn = iragn | s06bonus | opt_bpt_tradagn | optsfagn

#bpt['agn'] = ~(bpt['defstarform'])
#bpt['name'] = bpt['galname']
#s06['agn'] = ~(s06['defstarform'])
#s06['name'] = s06['galname']
#midir['agn'] = midir['agnflag']

###############################################################################
#stdagn = np.intersect1d(bpt.galname[~(bpt['sftoagn'])], bpt.galname[bpt['agn']])
#commonstdagn = s06['agn'] & ~bpt['agn']
#s06bonusagn = np.array(s06.galname[~commonstdagn & s06['agn']])
#sfagn = np.array(bpt.galname[bpt['sftoagn']])
#sf = np.array(list(bpt.galname[bpt['defstarform']]) + list(s06.galname[s06['defstarform']]))
#stdagn = np.array(list(bpt.galname[bpt['composite'] | bpt['defagn'] | bpt['agntosf']]) + \
#                        list(s06.galname[s06['defagn']]))
#s06bonusagn = np.array(s06.galname[s06['composite']])
#sfagn = np.array(bpt.galname[bpt['sftoagn']])
#midiragn = np.array(midir.name[midir['agn']])


###############################################################################

#for x in resdf.keys():
#    resdf[x] = np.array(resdf[x]).byteswap().newbyteorder() 
plot = 0

if plot:
    xmin = 7.5
    xmax = 11.5
    ymin = 0.45
    ymax = 2.75
    u_r = resdf['modelu_r']
    lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    sf_colors_map = truncate_colormap(cm.gray,n=11)
    nbins = 100
    msize = 8
    fig,ax = plt.subplots()#(figsize=(8,8))
    
    X,Y,Z = density_estimation(resdf.logmstar,u_r)
    xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    resdfcontour = np.column_stack((resdf.logmstar, u_r))
    k = kde.gaussian_kde(resdfcontour.T)
    resdf_contour_z = k(np.vstack([xgrid.flatten(), ygrid.flatten()]))
    resdf_colors_map = truncate_colormap(cm.gray_r)
    ax.pcolormesh(xgrid, ygrid, resdf_contour_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=resdf_colors_map) #plt.cm.gray_r)
    #ax.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 0)
    
    sel = resdf.loc[bpt.galname]
    ax.plot(sel.logmstar, sel.modelu_rcorr, '1', color = 'red', markersize = msize, 
           mew = 2, label = 'Traditional/S06 search sample', zorder = 3)
    
    sel = resdf.loc[opt.galname]
    ax.plot(sel.logmstar, sel.modelu_rcorr, 's', color = 'blue', markersize = msize+2, 
            mew = 2, mfc = 'none', label = 'SF-AGN search sample (P22)', zorder = 2)
    plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 20)
    plt.ylabel('u-r', fontsize = 20)
    plt.legend(fontsize = 18)
    
    fig,ax = plt.subplots()#(figsize=(8,8))
    ax.pcolormesh(xgrid, ygrid, resdf_contour_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=resdf_colors_map) #plt.cm.gray_r)
    
    sel = resdf.loc[irrand.name]
    ax.plot(sel.logmstar, sel.modelu_rcorr, 'o', color = 'green', markersize = msize, 
            mew = 2, mfc = 'none', label = 'Mid-IR search sample', zorder = 14)
    #sel = resdf.loc[irsys.name]
    #ax.plot(sel.logmstar, sel.modelu_rcorr, 'o', color = 'orange', markersize = msize+2, 
    #        mew=2, mfc = 'none', label = 'Mid-IR systematic-error-limited sample', zorder = 5)
    #ax.axhline(y = 1.4)
    
    plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 20)
    plt.ylabel('u-r', fontsize = 20)
    plt.legend(fontsize = 18)
    #ax.contour(X, Y, Z, cmap='summer')
    #
    
    
    
    xmin = 7.5
    xmax = 11.5
    ymin = 0.45
    ymax = 2.75
    u_r = resdf['modelu_r']
    lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    sf_colors_map = truncate_colormap(cm.gray,n=11)
    nbins = 100
    msize = 12
    fig,ax = plt.subplots()#(figsize=(8,8))
    
    X,Y,Z = density_estimation(resdf.logmstar,u_r)
    xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    resdfcontour = np.column_stack((resdf.logmstar, u_r))
    k = kde.gaussian_kde(resdfcontour.T)
    resdf_contour_z = k(np.vstack([xgrid.flatten(), ygrid.flatten()]))
    resdf_colors_map = truncate_colormap(cm.gray_r)
    ax.pcolormesh(xgrid, ygrid, resdf_contour_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=resdf_colors_map) #plt.cm.gray_r)
    #ax.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 0)
    
    sel = resdf.loc[opt_bpt_tradagn]
    ax.plot(sel.logmstar, sel.modelu_r, '1', color = 'red', markersize = msize, 
           mew = 2, label = 'Traditional AGN', zorder = 3)
    
    sel = resdf.loc[s06bonus]
    ax.plot(sel.logmstar, sel.modelu_r, '1', color = 'lime', markersize = msize, 
           mew = 2, label = 'S06 Bonus AGN', zorder = 3)
    
    sel = resdf.loc[optsfagn]
    ax.plot(sel.logmstar, sel.modelu_r, 's', color = 'blue', markersize = msize, 
            mew = 2, mfc = 'none', label = 'SF-AGN', zorder = 2)
    
    sel = resdf.loc[iragn]
    ax.plot(sel.logmstar, sel.modelu_r, 'p', color = 'orange', markersize = msize, 
            mew = 2, mfc = 'none', label = 'Mid-IR AGN', zorder = 10)
    
    
    #sel = resdf.loc[dwarfagn]
    #ax.plot(sel.logmstar, sel.modelu_rcorr, 'o', color = 'k', markersize = msize+2, 
    #       mew = 2, mfc = 'none', label = 'Dwarf AGN', zorder = 5)
    
    #sel = resdf.loc[irsys.name]
    #ax.plot(sel.logmstar, sel.modelu_rcorr, 'o', color = 'orange', markersize = msize+2, 
    #        mew=2, mfc = 'none', label = 'Mid-IR systematic-error-limited sample', zorder = 5)
    
    ax.axvline(x = 9.5, c = 'k', linestyle = '--')
    
    plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 20)
    plt.ylabel('u-r', fontsize = 20)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.legend(fontsize = 15)
    #ax.contour(X, Y, Z, cmap='summer')
    #
    
    xmin = 9.18
    xmax = 12.0
    ymin = -3
    ymax = 2
    resdf['g_s'] = resdf.logmgas - resdf.logmstar
    resdf['mbary'] = np.log10(10**resdf.logmstar + 10**resdf.logmgas)

    fig,ax = plt.subplots()#(figsize=(8,8))
    X,Y,Z = density_estimation(resdf.mbary,resdf['g_s'])
    xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    resdfcontour = np.column_stack((resdf.mbary, resdf['g_s']))
    k = kde.gaussian_kde(resdfcontour.T)
    resdf_contour_z = k(np.vstack([xgrid.flatten(), ygrid.flatten()]))
    resdf_colors_map = truncate_colormap(cm.gray_r)
    ax.pcolormesh(xgrid, ygrid, resdf_contour_z.reshape(xgrid.shape), 
                       shading='gouraud', cmap=resdf_colors_map) #plt.cm.gray_r)
    #ax.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 0)
    
    sel = resdf.loc[opt_bpt_tradagn]
    ax.plot(sel.mbary, sel.g_s, '1', color = 'red', markersize = msize, 
           mew = 2, label = 'Traditional AGN', zorder = 3)
    
    sel = resdf.loc[s06bonus]
    ax.plot(sel.mbary, sel.g_s, '1', color = 'lime', markersize = msize, 
           mew = 2, label = 'S06 Bonus AGN', zorder = 3)
    
    sel = resdf.loc[optsfagn]
    ax.plot(sel.mbary, sel.g_s, 's', color = 'blue', markersize = msize, 
            mew = 2, mfc = 'none', label = 'SF-AGN', zorder = 2)
    
    sel = resdf.loc[iragn]
    ax.plot(sel.mbary, sel.g_s, 'p', color = 'orange', markersize = msize, 
            mew = 2, mfc = 'none', label = 'Mid-IR AGN', zorder = 10)
    
    
    ax.axvline(x = 9.9, c = 'k', linestyle = '--')
    ax.axhline(y = 0, c = 'k', linestyle = '--')
    
    plt.xlabel(r'log(M${_{baryonic}}$/M${_\odot}$)', fontsize = 20)
    plt.ylabel(r'log(M${_{gas}}$/M${_{stellar}}$)', fontsize = 20)
    plt.xlim(9.18,xmax)
    plt.ylim(-2.5,ymax)
    plt.legend(fontsize = 15)

#    fig,ax = plt.subplots()#(figsize=(8,8))
#    X,Y,Z = density_estimation(resdf.logmstar,resdf['g_s'])
#    xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
#    resdfcontour = np.column_stack((resdf.logmstar, resdf['g_s']))
#    k = kde.gaussian_kde(resdfcontour.T)
#    resdf_contour_z = k(np.vstack([xgrid.flatten(), ygrid.flatten()]))
#    resdf_colors_map = truncate_colormap(cm.gray_r)
#    ax.pcolormesh(xgrid, ygrid, resdf_contour_z.reshape(xgrid.shape), 
#                       shading='gouraud', cmap=resdf_colors_map) #plt.cm.gray_r)
#    #ax.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 0)
#    
#    sel = resdf.loc[opt_bpt_tradagn]
#    ax.plot(sel.logmstar, sel.g_s, '1', color = 'red', markersize = msize, 
#           mew = 2, label = 'Traditional AGN', zorder = 3)
#    
#    sel = resdf.loc[s06bonus]
#    ax.plot(sel.logmstar, sel.g_s, '1', color = 'lime', markersize = msize, 
#           mew = 2, label = 'S06 Bonus AGN', zorder = 3)
#    
#    sel = resdf.loc[optsfagn]
#    ax.plot(sel.logmstar, sel.g_s, 's', color = 'blue', markersize = msize, 
#            mew = 2, mfc = 'none', label = 'SF-AGN', zorder = 2)
#    
#    sel = resdf.loc[iragn]
#    ax.plot(sel.logmstar, sel.g_s, 'p', color = 'orange', markersize = msize, 
#            mew = 2, mfc = 'none', label = 'Mid-IR AGN', zorder = 10)
#    
#    
#    ax.axvline(x = 9.5, c = 'k', linestyle = '--')
#    ax.axhline(y = 0, c = 'k', linestyle = '--')
#    
#    plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 20)
#    plt.ylabel(r'log(M${_{gas}}$/M${_{stellar}}$)', fontsize = 20)
#    plt.xlim(xmin,xmax)
#    plt.ylim(-2.5,ymax)
#    plt.legend(fontsize = 15)
    
    
    from scipy.stats import kde
    import matplotlib
    import matplotlib.cm as cm
    import matplotlib.colors as mpcolors
    
    def density_estimation(m1, m2):
        X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
        positions = np.vstack([X.ravel(), Y.ravel()])                                                       
        values = np.vstack([m1, m2])                                                                        
        kernel = scipy.stats.gaussian_kde(values)                                                                 
        Z = np.reshape(kernel(positions).T, X.shape)
        return X, Y, Z 
    
    plt.figure()
    ymin, ymax = (7.8, 9.2)
    xmin,xmax = (7.5,11.5)
    
    ax1 = plt.subplot()
    
    N2 = np.log10(resdf.nii_6584_flux/resdf.h_alpha_flux)
    bad = ((N2<-2.5) | (N2>-0.3))
    Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
    Z_pp04[bad] = np.nan
    
    #TODO -- check N2O3 or N2HA or Nebula Bayes with SEL or only BPT-lines
    
#ADD RIGHT AXES
    
    #Z = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_spectra/ECO+RESOLVE_snr5_jhu_csf_M5grid_final_LOGZ.txt", 
    #                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
    #Z.index = Z.name
    #Z_pp04 = Z.loc[resdf.name]['Estimate']+8.76-0.11
    
    
    yaxis = np.linspace(np.min(Z_pp04)-0.5,np.max(Z_pp04)+0.5)
    xaxis = np.linspace(7, 11.5)
    plt.plot(9.5*np.ones(len(yaxis)),yaxis, 'k-.', linewidth = 3)
    Z04 = np.log10(0.4)+8.76
    plt.plot(xaxis, Z04*np.ones(len(xaxis)), 'k--', linewidth = 2)
    plt.plot(xaxis, 8.76*np.ones(len(xaxis)), 'k--', linewidth = 2)
    ax2 = ax1.twinx()
    yticks = np.array([7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0,9.2])
    ax2.set_yticks(np.linspace(0,1,len(yticks)))
    ax1.set_yticks(yticks)#np.arange(7.8,9.2,0.2))
    float_formatter = lambda x: "%.2f" % x
    xticks = np.array([7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2])
    N2 = 1.754*yticks - 15.614
    N2_label = ["%.2f" % z for z in N2]
    ax2.set_yticklabels(N2_label)
    ax2.set_ylabel(r'[NII]/H$\alpha$')#, fontsize = 22)
    
    #plt.ylim(min(yaxis),max(yaxis))
    ax1.set_ylim(7.8,9.2)
    plt.xlim(7.5,11.5)
    ax1.set_xlabel(r'log(M${\rm_{stellar}}$/M${_\odot}$)')#, fontsize = 22)
    ax1.set_ylabel(r'Z (12 + log(O/H))')#, fontsize = 22)
    real = np.isfinite(Z_pp04)
    realdf = resdf[real]
    real_Z_pp04 = Z_pp04[real]
    X,Y,Z = density_estimation(realdf.logmstar.loc[opt_bpt_sf & set(realdf.index)],
                              real_Z_pp04.loc[opt_bpt_sf & set(realdf.index)])
    Z = Z/Z.max()
    
    nbins = 100
    #stdagn_contour = np.column_stack((realdf.logmstar.loc[np.intersect1d(stdagn, realdf.index)], 
    #                                  real_Z_pp04.loc[np.intersect1d(stdagn, realdf.index)]))
    stdagn_contour = np.column_stack((realdf.logmstar, real_Z_pp04))
    #xmin_agn = stdagn_contour[:,0].min(); xmax_agn = stdagn_contour[:,0].max()
    #ymin_agn = stdagn_contour[:,1].min(); ymax_agn = stdagn_contour[:,1].max()
    xmin_agn = xmin; xmax_agn = xmax
    ymin_agn = ymin; ymax_agn = ymax
    xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                    ymin_agn:ymax_agn:nbins*1j]
    k = kde.gaussian_kde(stdagn_contour.T)
    agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
    agn_colors_map = truncate_colormap(cm.gray_r)
    ax1.pcolormesh(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 
                       shading='gouraud', cmap=agn_colors_map) #plt.cm.gray_r)
    
    sf_colors_map = truncate_colormap(cm.gray,n=11)
    nbins = 20
    lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    #ax1.contour(X, Y, Z, levels = lvls, colors = 'k', zorder = 1)
    
    ax1.plot(resdf.logmstar.loc[irrand.name[irrand['agnflag']]], 
             Z_pp04[irrand.name[irrand['agnflag']]], 
             'p', markersize = 10, label = 'Mid-IR AGN', mfc = 'none',
             color = 'orange', alpha = 1.0, mec = 'orange', mew = 2, zorder = 10)
    ax1.plot(resdf.logmstar.loc[s06bonus], Z_pp04.loc[s06bonus], '1', 
             markersize = 10, label = 'S06 Bonus AGN', 
             color = 'lime', alpha = 1.0, mew = 2, zorder = 4)
    ax1.plot(resdf.logmstar.loc[opt_bpt_tradagn], Z_pp04.loc[opt_bpt_tradagn], '1', markersize = 10, label = 'Traditional AGN', 
             color = 'red', alpha = 1.0, mec = 'red', mew = 2, mfc = 'none', zorder = 2)
    ax1.plot(resdf.logmstar.loc[optsfagn], Z_pp04.loc[optsfagn], 's', markersize = 10, label = 'SF-AGN', 
             color = 'blue', alpha = 1.0, mec = 'blue', mew = 2, mfc = 'none', zorder = 3)
    ax1.legend(loc='lower right', fontsize = 18)
    
    
    
    ###############################################################################
    #Substitute RESOLVE SFR_NUV_WISE to SEL-DF SFR_NUV from database
    rescat = readsav("resolvecatalog_031622.dat",python_dict=True)
    rescatphot = readsav("resolvecatalogphot_021622.dat",python_dict=True)
    match, rescatndx, resselndx = np.intersect1d(rescat['name'], resdf.name, return_indices=True)
    rescatsel_sfr = rescat['sfr_nuv_wise'][rescatndx]
    resdf.loc[resdf.name.iloc[resselndx],'sfr_nuv_wise'] = rescatsel_sfr
    
    ecocat = pd.read_csv("sfr_nuv_wisew_ecoSEL_newsnr.txt")
    ecocat['sfr_nuv_wisew'][(ecocat['mw4w'] < 0)] = ecocat['sfr_nuv'][(ecocat['mw4w'] < 0)]
    ecocat['sfr_nuv_wisew'][(ecocat['mw3w'] < 0)] = ecocat['sfr_nuv'][(ecocat['mw3w'] < 0)]
    match, ecocatndx, ecoselndx = np.intersect1d(ecocat['econame'], resdf.name, return_indices=True)
    ecocatsel_sfr = np.array(ecocat['sfr_nuv_wisew'][ecocatndx])
    ecocatsel_sfr_err = np.array(ecocat['sfr_nuv_wisew_err'][ecocatndx])
    resdf.loc[resdf.name.iloc[ecoselndx],'sfr_nuv_wise'] = ecocatsel_sfr
    
    
    ssfr_st = 10**(np.log10(resdf.sfr_nuv_wise) - resdf.logmstar)
    ssfr_lt = 1/(1+(1/(resdf.meanfsmgr)))
    #fsmgr_st = 100*(10**6)*(ssfr_st)/(0.1*1e9*(1-ssfr_st)) #using ssfr instead of sfr
    fsmgr_st = np.log10(100*(1e6)*(resdf.sfr_nuv_wise)/((10**resdf.logmstar - (100*1e6*resdf.sfr_nuv_wise)) *0.1*1e9))
    fsmgr_lt = np.log10(resdf.meanfsmgr)
    
    plt.figure()
    #xaxis = np.linspace(10**-12.3, 1.2*np.max(fsmgr_st),num=50)
    xaxis = np.linspace(-12.3, np.max(fsmgr_st),num=50)
    yaxis = np.zeros(len(xaxis))
    
    xmin = np.min(xaxis); xmax = np.max(xaxis)
    ymin = np.min(fsmgr_lt); ymax = np.max(fsmgr_lt)
    #ymin = 0.0001; ymax = 15
    sf = opt_bpt_sf | set(irrand.name[~irrand.agnflag])
    sf_contour = np.column_stack((fsmgr_st.loc[sf], fsmgr_lt.loc[sf]))
    xmin_sf = xmin; xmax_sf = xmax
    ymin_sf = ymin; ymax_sf = ymax
    xgrid_sf, ygrid_sf = np.mgrid[xmin_sf:xmax_sf:nbins*1j, 
                                    ymin_sf:ymax_sf:nbins*1j]
    k = kde.gaussian_kde(sf_contour.T)
    sf_contour_z = k(np.vstack([xgrid_sf.flatten(), ygrid_sf.flatten()]))
    nbins = 20
    lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    #plt.contour(xgrid_sf, ygrid_sf, sf_contour_z.reshape(xgrid_sf.shape), 
    #            colors = 'k')
    sf_colors_map = truncate_colormap(cm.gray_r)
    plt.pcolormesh(xgrid_sf, ygrid_sf, sf_contour_z.reshape(xgrid_sf.shape), 
                       shading='gouraud', cmap=sf_colors_map) 
    
    nbins = 30
    #stdagn_contour = np.column_stack((fsmgr_st.loc[optsfagn], 
    #                                  fsmgr_lt.loc[optsfagn]))
    #xmin_agn = xmin; xmax_agn = xmax
    #ymin_agn = ymin; ymax_agn = ymax
    #xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
    #                                ymin_agn:ymax_agn:nbins*1j]
    #k = kde.gaussian_kde(stdagn_contour.T)
    #agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
    #agn_colors_map = truncate_colormap(cm.gray_r, n=5)
    plt.plot(fsmgr_st.loc[opt_bpt_tradagn], fsmgr_lt.loc[opt_bpt_tradagn], 
                   '1', color = 'red', label = 'Traditional AGN', 
                   markersize = 10, mew = 2, alpha = 1, zorder = 4)
    plt.plot(fsmgr_st.loc[s06bonus], fsmgr_lt.loc[s06bonus], 
                   '1', color = 'lime', label = 'S06 AGN', markersize = 10, mew = 2,
                   alpha = 1, zorder = 4)
    plt.plot(fsmgr_st.loc[optsfagn], fsmgr_lt.loc[optsfagn], 's', color = 'blue', 
             mec = 'blue', mfc = 'none', label = 'SF-AGN', 
             markersize = 10, mew = 2, alpha = 1, zorder = 5)
    plt.plot(fsmgr_st.loc[iragn], fsmgr_lt.loc[iragn], 
                   'p', color = 'orange', label = 'Mid-IR AGN', markersize = 10, 
                   mew = 3, mfc = 'none', zorder = 10)
    plt.plot(xaxis, yaxis, 'k--', lw = 3)
    plt.xlim(np.min(xaxis), np.max(fsmgr_st))
    plt.ylim(np.min(fsmgr_lt), np.max(fsmgr_lt))
    plt.text(0.0005, 1.25, r'Stellar Mass Doubled in last Gyr', 
                 fontsize=14, color='k')
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend(fontsize = 15, loc = 'lower right')
    plt.tick_params(which = 'minor', length=4, width=2)
    
    #plt.plot(fsmgr_st.loc[jwst], fsmgr_lt.loc[jwst],
    #            '*',color = 'lime', ms = 10, label = labels[key])
    plt.xlabel(r'Short Term SFH $\left(\frac{M_*(<100 Myr)}{M_*(>100 Myr) \times 0.1 Gyr}\right)$')
    plt.ylabel(r'Long Term SFH $\left(\frac{M_*(<1 Gyr)}{M_*(>1 Gyr) \times 1 Gyr}\right)$')
    
    
    plt.figure()
    
    mbary = np.log10(10**resdf.logmstar + 10**resdf.logmgas)
    #centrals
    plt.plot(resdf.logmh[resdf.fc == 1].loc[opt_bpt_tradagn], 
             mbary[resdf.fc == 1].loc[opt_bpt_tradagn], 
             '1', color = 'red', mew = 2, mfc = 'none', markersize = 10,
             label = 'Traditional AGN') #other group galaxies
    plt.plot(resdf.logmh[resdf.fc == 1].loc[s06bonus], 
             mbary[resdf.fc == 1].loc[s06bonus], 
             '1', color = 'lime', mew = 2, mfc = 'none', markersize = 10,label = 'S06 Bonus AGN') #other group galaxies
    plt.plot(resdf.logmh[resdf.fc == 1].loc[optsfagn], 
             mbary[resdf.fc == 1].loc[optsfagn], 
             's', color = 'blue', mew = 2, mfc = 'none', markersize = 10,label = 'SF-AGN') #other group galaxies
    plt.plot(resdf.logmh[resdf.fc == 1].loc[iragn], 
             mbary[resdf.fc == 1].loc[iragn], 
             'p', color = 'orange', mew = 2, mfc = 'none', markersize = 10,label = 'Mid-IR AGN') #other group galaxies
    plt.legend(loc = 'lower right', numpoints = 1, fontsize = 15)
#satellites
    plt.plot(resdf.logmh[resdf.fc == 0].loc[opt_bpt_tradagn], 
             mbary[resdf.fc == 0].loc[opt_bpt_tradagn], 
             '1', color = 'red', mew = 2, mfc = 'none', markersize = 5) #other group galaxies
    plt.plot(resdf.logmh[resdf.fc == 0].loc[s06bonus], 
             mbary[resdf.fc == 0].loc[s06bonus], 
             '1', color = 'lime', mew = 2, mfc = 'none',markersize = 5) #other group galaxies
    plt.plot(resdf.logmh[resdf.fc == 0].loc[optsfagn], 
             mbary[resdf.fc == 0].loc[optsfagn], 
             's', color = 'blue', mew = 2, mfc = 'none', markersize = 5) #other group galaxies
    plt.plot(resdf.logmh[resdf.fc == 0].loc[iragn], 
             mbary[resdf.fc == 0].loc[iragn], 
             'p', color = 'orange', mew = 2, mfc = 'none', markersize = 5, 
             zorder = 10) #other group galaxies

    plt.plot(resdf.logmh[resdf.fc == 1].loc[sf], 
             mbary[resdf.fc == 1].loc[sf], 
             '.', color = 'gray', mec= 'none', markersize = 10, alpha = 0.3,
             zorder = 0) 
    plt.plot(resdf.logmh[resdf.fc == 0].loc[sf], 
             mbary[resdf.fc == 0].loc[sf], 
             '.', color = 'gray', mec = 'none', markersize = 5, alpha = 0.3,
             zorder = 0) 

    plt.xlabel(r'log(Galaxy M$\rm_{halo}$/M$_{\odot}$)')
    plt.ylabel(r'log(Galaxy M$\rm_{bary}$/M$_{\odot}$)')
    plt.xlim(10.8,14.54)
    plt.ylim(9.15,11.5)
    fullmbary = np.log10(10**resdf.logmstar + 10**resdf.logmgas)
    yaxis = np.arange(np.min(fullmbary)-0.5, np.max(fullmbary)+0.5,0.1)
    xaxis = np.ones(len(yaxis))
    plt.plot(11.5*xaxis, yaxis, 'k--', lw = 3)
    plt.plot(12.0*xaxis, yaxis, 'k-.', lw = 3)

