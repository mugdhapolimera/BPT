# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 11:25:51 2020

@author: mugdhapolimera
"""

# -*- coding: utf-8 -*-
"""
Created on Sun May 26 15:38:42 2019

@author: mugdhapolimera
"""

#This program makes a Line-Ratio diagram (also known as a BPT plot or Kewley diagram)
#with labels using data from the RESOLVE survey to classify galaxies as LINERs,
#Seyferts, Composites, or AGNs on the basis of their flux ratio for distinction.

#Original code from Ashley Bittner 08/03/2017
#Edited version from Margie Bruff 01/07/2018
#Updated by Carlynn Ferguson 03/30/2018

#suggested use of python debugger to understand the code more thoroughly
#KDE plotting 
#https://python-graph-gallery.com/86-avoid-overlapping-in-scatterplot-with-2d-density/
#import pdb

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import sys
from scipy.stats import kde
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as colors
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
from scipy.optimize import curve_fit

def fit_bptline(x,a,b,c,d,e):
    y = (a + b*x + c*(x**2))*(np.tanh(d*x)) + e
    return y
#define demarcation function: log_NII_HA vs. log_OIII_HB
def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2halocus(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.1 + (0.61 / (log_NII_HA + 0.08))
def n2hamain(x): #main line for NII/H-alpha from equation 5, Kewley 2006
#    return 1.19 + (0.61 / (log_NII_HA - 0.47)) #Ke01
#    return 0.57 + (0.13 / (log_NII_HA - 0.003))
    log_NII_HA= x.copy()
#    OvHb = y.copy()
#    log_NII_HA[log_NII_HA > -0.4] = -0.4
    return (-30.787 + 1.1358*log_NII_HA + 0.27297*log_NII_HA**2)*np.tanh(5.7409*log_NII_HA) - 31.093
#    return fit_bptline(NvH,*fn[0])
def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
    return 1.30 + (0.72 / (log_SII_HA - 0.32))
#    return 0.58 + (0.04 / (log_SII_HA +0.012))
def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
    return 0.76 + 1.89*log_SII_HA
def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
    return 1.33 + (0.73 / (log_OI_HA + 0.59))
#    return 0.61 + (0.056 / (log_OI_HA + 0.40))
def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
    return 1.3 + 1.18*log_OI_HA
def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
    return -0.59
def n2hacrit(log_NII_HA): #boundary for OI/H-alpha
    return -0.077
def s2hacrit(log_SII_HA): #boundary for OI/H-alpha
    return -0.033

def he2hbmain(log_NII_HA):
    return -1.22+1.0/(8.92*log_NII_HA+1.32)

def he2hbmainclass(log_NII_HA):
    refn2ha = np.linspace(-3.0, -0.15)
    main = he2hbmain(refn2ha)
    line = np.poly1d(np.polyfit(refn2ha, main, 15))
    return line(log_NII_HA) 

def he2hblimit(log_NII_HA):
        return -1.07+1.0/(8.92*log_NII_HA-0.95)

def he2hblimitclass(log_NII_HA):
    refn2ha = np.linspace(-3.0, -0.15)
    limit = he2hblimit(refn2ha)
    line = np.poly1d(np.polyfit(refn2ha, limit, 15))
    return line(log_NII_HA) 

def ratioerror(num,num_err,den, den_err):
    err_num2 = (num_err/(num*np.log(10)))**2
    err_den2 = (den_err/(den*np.log(10)))**2
    return np.sqrt(err_num2 + err_den2)

def s06_bpt(inputfile, outputfile, eco, resolve, full, sdsscat, save, ax):
    df = pd.read_csv(inputfile)
    #define alternate catalog names
    if 'name' in df.keys():
        df['NAME'] = df['name']
    #if (sdsscat == 'nsa') & ('resname' in df.keys()):
    #    df['name'] = df['resname']
    if 'CATID' in df.keys():
        df['NAME'] = df['CATID']
    name = df['name']
    df['NAME'] = df['name']
    df.index = df.name
    
    if eco or full: 
        resname = df['resname'] #for eco
        resname = resname != 'notinresolve'
    if resolve:
        econame = df['NAME']#df['econame'] #for resolve
        econame = df['NAME']#econame != 'notineco'

    # note that the ratio uses only the stronger line, but for S/N reasons we add
    # the weaker and multiply by 3/4 since Chris Richardson says the canonical
    # line ratio is 3:1 (this needs to be updated with a more precise number)
    if sdsscat == 'port':
        nii = df['Flux_NII_6583']
        nii_sum = (df['Flux_NII_6583']+ df['Flux_NII_6547'])*3./4
        nii_sum_err = (np.sqrt(df['Flux_NII_6547_Err']**2 + df['Flux_NII_6583_Err']**2))*3./4
        oiii = df['Flux_OIII_5006']
        oiii_err = df['Flux_OIII_5006_Err']
        h_alpha = df['Flux_Ha_6562']
        h_alpha_err = df['Flux_Ha_6562_Err']
        h_beta = df['Flux_Hb_4861']
        h_beta_err = df['Flux_Hb_4861_Err']
        oi = df['Flux_OI_6300']
        oi_err = df['Flux_OI_6300_Err']
        sii_sum = df['Flux_SII_6716'] + df['Flux_SII_6730']
        sii_sum_err = np.sqrt(df['Flux_SII_6716_Err']**2 + df['Flux_SII_6730_Err']**2)
    
    elif sdsscat == 'jhu' or sdsscat == 'nsa' or sdsscat == 'master':
        
        nii = df['nii_6584_flux']
        if 'nii_6548_flux' in df.keys():
            nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
            nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
        else:        
            nii_sum = df['nii_6584_flux']
            nii_sum_err = df['nii_6584_flux_err']
        
        # note that the ratio uses only the stronger line, but for S/N reasons we add
        # the weaker and multiply by 3/4 since Chris Richardson says the canonical
        # line ratio is 3:1 (this needs to be updated with a more precise number)
        oiii = df['oiii_5007_flux']
        oiii_err = df['oiii_5007_flux_err']
        h_alpha = df['h_alpha_flux']
        h_alpha_err = df['h_alpha_flux_err']
        h_beta = df['h_beta_flux']
        h_beta_err = df['h_beta_flux_err']
        oi = df['oi_6300_flux']
        oi_err = df['oi_6300_flux_err']
        if 'sii_6717_flux' in df.keys():
            sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
            sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
        else:
            sii_sum = df['sii_6731_flux']
            sii_sum_err = df['sii_6731_flux_err']
    
    gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & \
                (h_beta > 0) & (h_beta_err > 0) & \
                (h_alpha_err > 0) & (nii_sum_err > 0) & (oiii_err > 0))
    
    snr = ((h_alpha > 5*h_alpha_err) & (nii_sum > 5*nii_sum_err) & (oiii > 5*oiii_err) \
           & (h_beta > 5*h_beta_err))
    print(np.sum(gooddata))
    print(np.sum(snr))
    data = gooddata #& snr #use ALL galaxy data within catalog
    #if sdsscat == 'master':
    #data = df.name > 0
    #print total points shared with alternate catalog
    if full: 
        sel = (np.where(data & name)[0]) #for eco
    elif eco: 
        sel = (np.where(data & resname)[0]) #for eco
    elif resolve:
        sel = (np.where(data & econame)[0]) #for resolve
    print ''
    print 'TOTAL DATA WITH ALTERNATE CATALOG NAME: ', len(sel)
    df = df[data]
    
    
    nii = nii[data]
    nii_sum = nii_sum[data]
    nii_sum_err = nii_sum_err[data]
    oiii = oiii[data]
    oiii_err = oiii_err[data]
    oi = oi[data]
    oi_err = oi_err[data]
    sii_sum = sii_sum[data]
    sii_sum_err = sii_sum_err[data]
    h_beta = h_beta[data]
    h_beta_err = h_beta_err[data]
    h_alpha = h_alpha[data]
    h_alpha_err = h_alpha_err[data]
    
    #nii = df['nii_6584_flux']
    ##if 'nii_6548_flux' in df.keys():
    #nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
    #nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
    ##else:
    ## note that the ratio uses only the stronger line, but for S/N reasons we add
    ## the weaker and multiply by 3/4 since Chris Richardson says the canonical
    ## line ratio is 3:1 (this needs to be updated with a more precise number)
    #oiii = df['oiii_5007_flux']
    #oiii_err = df['oiii_5007_flux_err']
    #h_alpha = df['h_alpha_flux']
    #h_alpha_err = df['h_alpha_flux_err']
    #h_beta = df['h_beta_flux']
    #h_beta_err = df['h_beta_flux_err']
    #oi = df['oi_6300_flux']
    #oi_err = df['oi_6300_flux_err']
    #if 'sii_6717_flux' in df.keys():
    #    sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
    #
    #    sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
    #else:
    #    sii_sum = df['sii_6731_flux']
    #
    #    sii_sum_err = df['sii_6731_flux_err']
    
    if sdsscat =='nsa':
        df = df[data]
    
    #length of data to be used for debugging
    datalen = np.sum(data)
    
    # data ratios
    #n2ha = np.log10(nii_sum/h_alpha)
    o3hb = np.log10(oiii/h_beta) # always the y-axis
    o1ha = np.log10(oi/h_alpha)
    s2ha = np.log10(sii_sum/h_alpha)
    n2ha = np.log10(nii_sum/h_alpha)
    s2ha_err = np.array(ratioerror(sii_sum, sii_sum_err, h_alpha, h_alpha_err))
    o3hb_err = np.array(ratioerror(oiii, oiii_err, h_beta, h_beta_err))
    o1ha_err = np.array(ratioerror(oi, oi_err, h_alpha, h_alpha_err))
    n2ha_err = np.array(ratioerror(nii_sum, nii_sum_err, h_alpha, h_alpha_err))
        

    #Below are the selectors for the data to distinguish btwn: Seyferts, Composites,
    #and AGN's based on the flux ratio diagnostic as understood via Kewley 2006.
    
    #NII plot selectors
    compsel = (o3hb <= n2hacompmin(n2ha)) & (o3hb >= n2hamain(n2ha)) & (n2ha < 0.)
    sfsel = (o3hb < n2hamain(n2ha)) #& (n2ha < 0.) & ~(o3hb > n2hamain(n2ha)) #~(o3hb > n2hamain(n2ha)) & ~compsel1
    agnsel= (o3hb > n2hacompmin(n2ha)) | (n2ha > 0.)
    #plt.hist(o1ha_err[o1ha_err < 1e5], bins = 'fd')
    #SII plot selectors
    #sfsel2 = (o3hb <= s2hamain(s2ha)) & ~compsel1
    #seyfsel2 = ((o3hb > s2hamain(s2ha)) & (o3hb >= s2halinseyf(s2ha)))
    #linersel2 = ((o3hb > s2hamain(s2ha)) & (o3hb < s2halinseyf(s2ha)))
    #agnsel2 = (o3hb > s2hamain(s2ha)) & ~compsel1
    #
    ##OI plot selectors
    #sfsel3 = (o3hb <= o1hamain(o1ha)) & (o1ha < -0.7) & ~compsel1
    #seyfsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb >= o1halinseyf(o1ha))
    #linersel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb < o1halinseyf(o1ha))
    #agnsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & ~compsel1
    
    
    #REFERENCE for cumulative plot selectors
    #seyfselr = seyfsel2 & seyfsel3
    #linerselr = linersel2 & linersel3
    
    #cumulative plot selectors
    #compsel = compsel1  #composite galaxies
    #seyfsel = agnsel1 #& seyfselr #Seyfert AGN galaxies
    #linersel = agnsel1 #& linerselr #LINER AGN galaxies
    #ambigsel1 = sfsel1 & (agnsel2 | agnsel3) #SF in first plot, AGN in subsequent plot
    #ambigsel2 = np.array(agnsel1) & (np.array(sfsel2) | np.array(sfsel3)) #AGN in first plot, SF in subsequent plot
    #ambagnsel = agnsel1 & ~seyfselr & ~linerselr & ~(sfsel2 | sfsel3) #Ambiguous AGN
    #
    #sftoagn1 = sfsel1 & agnsel2
    #sftoagn2 = sfsel1 & agnsel3
    
    #Save the BPT flags to a CSV file
    #emlineclass = sfsel ^ compsel ^ seyfsel ^ linersel ^ ambigsel1 ^ ambigsel2 ^ ambagnsel
    #defagn = seyfsel | linersel | ambagnsel
    #if not he2_flag:    
    #    flags = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
    #                          'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
    #                          'sftoagn':ambigsel1, 'agntosf':ambigsel2, 'defagn': defagn,
    #                          'sftoagn1':sftoagn1, 'sftoagn2': sftoagn2})
    #else:
    flags = pd.DataFrame({'galname':np.array(df.name), 'defstarform':sfsel, 'composite':compsel, 
                              'defagn':agnsel})
            
    flags.to_csv(outputfile ,index=False)
    
    #checking that plotted points are within the total data range
    print ''
    sfselpts = (len(np.where(sfsel)[0]))
    #seyfselpts = (len(np.where(seyfsel)[0]))
    #linerselpts = (len(np.where(linersel)[0]))
    compselpts = (len(np.where(compsel)[0]))
    agnselpts = (len(np.where(agnsel)[0]))
    #ambigsel1pts = (len(np.where(ambigsel1)[0]))
    #ambigsel2pts = (len(np.where(ambigsel2)[0]))
    totalselpts = sfselpts+compselpts+agnselpts
    #\seyfselpts+linerselpts+ambigsel1pts+ambigsel2pts
    sfpercent = float(sfselpts)/float(datalen)*100
    #seyfpercent = float(seyfselpts)/float(datalen)*100
    #linerpercent = float(linerselpts)/float(datalen)*100
    comppercent = float(compselpts)/float(datalen)*100
    agnpercent = float(agnselpts)/float(datalen)*100
    #ambig1percent = float(ambigsel1pts)/float(datalen)*100
    #ambig2percent = float(ambigsel2pts)/float(datalen)*100
    fulldf = df
    fulldf.index = fulldf.name
    if 'logmstar' in df.keys():
        dwarf = (df.logmstar[data] < 9.5)
        giant = (df.logmstar[data] > 9.5)
    else:
        dwarf = (fulldf.loc[df.name].logmstar < 9.5)
        giant = (fulldf.loc[df.name].logmstar > 9.5)
    agn = (agnsel|compsel) #(ambigsel1|seyfsel|linersel|ambagnsel|compsel)
    dwarfagn = dwarf & agn
    giantagn = giant & agn
    
    print ("DATA POINTS: "),datalen
    print ("TOTAL PLOTTED POINTS: "), totalselpts
    print ("TOTAL PLOTTED POINTS OMITTED: "), datalen-totalselpts
    print "* IF NUMBERS ABOVE ARE AT ALL NEGATIVE THEN THERE IS OVERPLOTTING"
    print ("Definite Star Forming: "),sfselpts,("("),round(sfpercent, 2),("%"),(")")
    print ("Composite: "),compselpts, ("("),round(comppercent, 2),("%"),(")")
    #print ("SF --> AGN: "), ambigsel1pts, ("("),round(ambig1percent, 2),("%"),(")")
    #print ("AGN --> SF: "), ambigsel2pts, ("("),round(ambig2percent, 2),("%"),(")")
    print ("Ambiguous AGN: "),agnselpts, ("("),round(agnpercent, 2),("%"),(")")
    #print ("Seyfert: "),seyfselpts, ("("),round(seyfpercent, 2),("%"),(")")
    #print ("LINER: "),linerselpts, ("("),round(linerpercent, 2),("%"),(")")
    #print ("TOTAL KNOWN AGN: "),linerselpts+seyfselpts+agnselpts, ("("), \
    print ("TOTAL KNOWN AGN: "),agnselpts, ("("), \
    #round(linerpercent+seyfpercent+agnpercent, 2), ("% )")
    round(agnpercent, 2), ("% )")
    #print ("POSSIBLE TOTAL AGN: "),linerselpts+seyfselpts+agnselpts+ambigsel1pts+ambigsel2pts,("("),\
    print ("POSSIBLE TOTAL AGN: "),+agnselpts+compselpts,("("),\
    round(comppercent+agnpercent, 2), ("% )")
    #round(linerpercent+seyfpercent+agnpercent+ambig1percent+ambig2percent, 2), ("% )")
    #print ("Percent Omitted: "), round((100-(sfpercent+seyfpercent+linerpercent+comppercent+agnpercent+ambig1percent+ambig2percent)), 2), ("%")
    print ''
    
    print ("Dwarf AGN: "), np.sum(dwarfagn)
    print ("AGN % in Dwarf Galaxies: "), round(100.0*np.sum(dwarfagn)/float(np.sum(dwarf)),2), ("%")
    #print ("AGN in Giant Galaxies: "), 100*round(np.sum(giantagn)/float(np.sum(giant)),2), ("%")
    #print ("SF-AGN in dwarfs: "), np.sum(ambigsel1 & dwarf)
    #print ("Number of Dwarfs:"), np.sum(dwarf)
    ###PLOTS###
    #reference points in x-direction for demarcation lines on plots
    def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
      	new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
      	return new_cmap
    sf_colors_map = truncate_colormap(cm.gray_r)

    refn2ha = np.linspace(-3.0, 0.35)
#    refoiha = np.linspace(-2.5, -0.4)
#    refsiiha = np.linspace(-2, 0.3,100)
    
    #lowsfagn = ['rf0376', 'rf0503', 'rs0063', 'rs0626', 'rs1195', 'rs1292']
    #NII/OIII plot
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    xmin = refn2ha.min(); xmax = refn2ha.max()
    ymin = -1.0; ymax = 1.0
    nbins = 50

    definite = np.column_stack((n2ha[sfsel], o3hb[sfsel]))
    xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    k2 = kde.gaussian_kde(definite.T)
    definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
    ax.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
                   shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
    ax.set_xlim(-1.5,0.32)
    ax.set_ylim(-1.0,1.0)
    main1, = ax.plot(refn2ha, n2hamain(refn2ha), 'g--', lw = 5)#, label = 'Ke01 Maximum Starburst Line')
    composite, = ax.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
                          '-.', color = 'cyan', lw = 5)#, label = 'Ka03 Composite Line')
#    sfsel1, = ax.plot(n2ha[sfsel], o3hb[sfsel], 'k.', alpha = 0.1, 
#                       markersize = 5)#, label = 'Definite Star Forming')
    compdata1, = ax.plot(n2ha[compsel], o3hb[compsel], '1', color = 'lime', 
                          markersize = 8, mew = 2, alpha = 0.5, label = 'S06 Bonus AGN')
    agnsel1, = ax.plot(n2ha[agnsel], o3hb[agnsel], 'r1',
                       markersize = 8, mew = 2, alpha = 0.5, label = 'Conventional AGN')
    dwarfagn1, = ax.plot(n2ha[dwarfagn], o3hb[dwarfagn], 'kv', mfc = 'none',
                       mec = 'k', mew = 2,  markersize = 12, label = 'Dwarf AGN')
            
    ax.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
    ax.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
    ax.legend(loc = 'lower left', fontsize = 15)


    
    #SII/OIII plot
#    plt.figure('SII Scatter Plot')
#    ax2 = plt.subplot(111)
#    main2, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
#    ax2.set_xlim(-1.5, 0.5)
#    ax2.set_ylim(-1.0,1.0)
#    ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#    ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#    sfsel1, = ax2.plot(s2ha[sfsel], o3hb[sfsel], 'k.', alpha = 0.1, 
#                       markersize = 5)#, label = 'Definite Star Forming')
#    compdata1, = ax2.plot(s2ha[compsel], o3hb[compsel], 'm1', 
#                          markersize = 8, mew = 2, alpha = 0.5, label = 'S06 Composite')
#    agnsel1, = ax2.plot(s2ha[agnsel], o3hb[agnsel], 'r1',
#                       markersize = 8, mew = 2, alpha = 0.5, label = 'S06 Conventional AGN')
#    dwarfagn1, = ax2.plot(s2ha[dwarfagn], o3hb[dwarfagn], 'kv', mfc = 'none',
#                       mec = 'k', mew = 2,  markersize = 12, label = 'S06 Dwarf AGN')
#    #
    ##OI/OIII plot
#    plt.figure('OI Scatter Plot')
#    ax3 = plt.subplot(111)
#    main3, = ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
#                      'k', label = 'Ke01 Maximum Starburst Line')
#    comp3, = ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
#                      'k-.', label = 'Ka03 Composite Line')
#    liner2, = ax3.plot(refoiha[refoiha > -1.1], o1halinseyf(refoiha[refoiha > -1.1]),
#                       'k--', label = 'Ke06 Liner/Seyfert Division Line')
#    ax3.set_xlim(-2.0, -0.4)
#    ax3.set_ylim(-1.0,1.0)
#    ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
#    ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#    sfdata3, = ax3.plot(o1ha[sfsel], o3hb[sfsel], 'k.', alpha = 0.1, 
#                        markersize = 5, label = 'SF')
#    compdata3, = ax3.plot(o1ha[compsel], o3hb[compsel], 'ms',
#                          markersize = 8, mew = 0, label = 'Composite')
##    plt.legend(bbox_to_anchor=(1.25, 1),loc=2, borderaxespad=0., 
##               numpoints = 1, fontsize = 14)
#    
#    
#    Z = np.array([0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.5,2.5])
#    OvHb = (0.023606 - 0.667627*Z)*np.tanh(-3.412213 + 5.743451*Z) + 0.712143
#    NvH = (-1.0577 - 0.055221*Z)*np.tanh(2.00404 - 3.82832*Z) - 1.55079
#    ax.plot(NvH,OvHb,'r')
#    ax.scatter(NvH,OvHb,color = 'b')
#    #ndx = np.where(Z == 1.0)
#    #plt.scatter(NvH[ndx],OvHb[ndx])
#    
#    fn = curve_fit(fit_bptline, NvH, OvHb)
#    xaxis = np.arange(-1.6,-0.25,0.01)
#    ax.plot(refn2ha,fit_bptline(refn2ha,*fn[0]),'g') 
#    
#    OvHa = (-0.83751 + 0.110241*Z)*np.tanh(2.35279 - 3.97006*Z) - 2.11304
#    SvH = (-0.86928 + 0.052481*Z)*np.tanh(2.66503 - 4.44255*Z) - 1.251617
#    ax2.plot(SvH,OvHb,'r')
#    ax2.scatter(SvH,OvHb,color = 'b')
#
#    ax3.plot(OvHa,OvHb,'r')
#    ax3.scatter(OvHa,OvHb,color = 'b')

    return (ax)