# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 12:58:09 2022

@author: mugdhapolimera
"""

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

def ratioerror(num,num_err,den, den_err):
    err_num2 = (num_err/(num*np.log(10)))**2
    err_den2 = (den_err/(den*np.log(10)))**2
    return np.sqrt(err_num2 + err_den2)

def s06_bptagn(inputfile, outputfile, seloutputfile, midirfile, 
               eco, resolve, full, sdsscat, save, ax):
    sel = pd.read_csv(seloutputfile)
    sel.index = sel.galname
    selagn = np.array(sel[~sel.defstarform].galname)

    midir = pd.read_csv(midirfile)
    midir.index = midir.name
    midiragn = np.array(midir.name.iloc[np.where(midir.agnflag == True)[0]])
    
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
    
    if sdsscat =='nsa':
        df = df[data]
    
    #length of data to be used for debugging
    datalen = np.sum(data)
    
    # data ratios
    o3hb = np.log10(oiii/h_beta) # always the y-axis
    n2ha = np.log10(nii_sum/h_alpha)
    n2ha_err = np.array(ratioerror(nii_sum, nii_sum_err, h_alpha, h_alpha_err))
        

    #Below are the selectors for the data to distinguish btwn: Seyferts, Composites,
    #and AGN's based on the flux ratio diagnostic as understood via Stasinska 2006.
    
    #NII plot selectors
    compsel = (o3hb <= n2hacompmin(n2ha)) & (o3hb >= n2hamain(n2ha)) & (n2ha < 0.)
    sfsel = (o3hb < n2hamain(n2ha)) #& (n2ha < 0.) & ~(o3hb > n2hamain(n2ha)) #~(o3hb > n2hamain(n2ha)) & ~compsel1
    agnsel= (o3hb > n2hacompmin(n2ha)) | (n2ha > 0.)
    flags = pd.DataFrame({'galname':np.array(df.name), 'defstarform':sfsel, 'composite':compsel, 
                              'defagn':agnsel})
            
    flags.to_csv(outputfile ,index=False)
    
    #checking that plotted points are within the total data range
    print ''
    sfselpts = (len(np.where(sfsel)[0]))
    compselpts = (len(np.where(compsel)[0]))
    agnselpts = (len(np.where(agnsel)[0]))
    totalselpts = sfselpts+compselpts+agnselpts
    sfpercent = float(sfselpts)/float(datalen)*100
    comppercent = float(compselpts)/float(datalen)*100
    agnpercent = float(agnselpts)/float(datalen)*100
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
    print ("Ambiguous AGN: "),agnselpts, ("("),round(agnpercent, 2),("%"),(")")
    print ("TOTAL KNOWN AGN: "),agnselpts, ("("), \
    round(agnpercent, 2), ("% )")
    print ("POSSIBLE TOTAL AGN: "),+agnselpts+compselpts,("("),\
    round(comppercent+agnpercent, 2), ("% )")
    print ''
    
    print ("Dwarf AGN: "), np.sum(dwarfagn)
    print ("AGN % in Dwarf Galaxies: "), round(100.0*np.sum(dwarfagn)/float(np.sum(dwarf)),2), ("%")
    #print ("AGN in Giant Galaxies: "), 100*round(np.sum(giantagn)/float(np.sum(giant)),2), ("%")
    #print ("SF-AGN in dwarfs: "), np.sum(ambigsel1 & dwarf)
    #print ("Number of Dwarfs:"), np.sum(dwarf)

    refn2ha = np.linspace(-3.0, 0.35)
    refoiha = np.linspace(-2.5, -0.4)
    refsiiha = np.linspace(-2, 0.3,100)
    
    ax.set_xlim(-1.5,0.32)
    ax.set_ylim(-1.0,1.0)
    main1, = ax.plot(refn2ha, n2hamain(refn2ha), 'k')#, label = 'Ke01 Maximum Starburst Line')
    composite, = ax.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
                          'k-.')#, label = 'Ka03 Composite Line')
    sfsel1, = ax.plot(n2ha[sfsel], o3hb[sfsel], 'k.', alpha = 0.1, 
                       markersize = 5)#, label = 'Definite Star Forming')
    dwarfagn1, = ax.plot(n2ha[dwarfagn], o3hb[dwarfagn], 'kv', mfc = 'none',
                       mec = 'k', mew = 2,  markersize = 12, label = 'S06 Dwarf AGN',
                       zorder = 10)
    seldwarfagn = np.intersect1d(selagn, df.name[dwarf])
    seldwarfagn1, = ax.plot(n2ha.loc[seldwarfagn], o3hb.loc[seldwarfagn], 'ks', mfc = 'none',
                       mec = 'k', mew = 2,  markersize = 12, label = 'SEL Dwarf AGN',
                       zorder = 11)
    midirdwarfagn = np.intersect1d(midiragn, df.name[dwarf])
    midirdwarfagn1, = ax.plot(n2ha.loc[midirdwarfagn], o3hb.loc[midirdwarfagn], 'kp', mfc = 'none',
                       mec = 'k', mew = 2,  markersize = 12, label = 'Mid-IR Dwarf AGN',
                       zorder = 12)
    compdata1, = ax.plot(n2ha[compsel], o3hb[compsel], '1', color = 'lime',
                          markersize = 8, mew = 2, alpha = 0.5, label = 'S06 Composite')
    agnsel1, = ax.plot(n2ha[agnsel], o3hb[agnsel], 'r1',
                       markersize = 8, mew = 2, alpha = 0.5, label = 'S06 Traditional AGN')
    ax.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
    ax.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
    ax.legend(loc = 'lower left', fontsize = 15)
    
    return (ax)