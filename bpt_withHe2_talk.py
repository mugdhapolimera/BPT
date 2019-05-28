# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 16:05:06 2018

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
#from astroML.plotting import scatter_contour
from scipy.stats import kde
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as colors
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
#display catalog being used
print ''
#pdb.set_trace()

#read in data
he2_flag = 0
save = 1
resolve = 1
eco = 0
full = 0
if sys.platform == 'linux2':
        os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/SDSS_spectra/')

else:
    os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')


if full: 
    inputfile = 'ECO+RESOLVE_filter_new.pkl'
    print 'ECO RESULTS'
    if he2_flag:
        outputfile = 'eco+resolve_emlineclass_filter_he2.csv'
    else: 
        outputfile = 'eco+resolve_emlineclass_filter.csv'
elif eco: 
    inputfile = 'ECO_filter_new.pkl'
    print 'ECO RESULTS'
    if he2_flag:
        outputfile = 'eco_emlineclass_filter_he2.csv'
    else: 
        outputfile = 'eco_emlineclass_filter_new.csv'

else:
    inputfile = 'RESOLVE_filter_new.pkl'
    print 'RESOLVE RESULTS'
    if he2_flag:
        outputfile = 'resolve_emlineclass_filter_he2_new.csv'
    else: 
        outputfile = 'resolve_emlineclass_filter_new.csv'
df = pd.read_pickle(inputfile)

if ('heii_4685_flux_port' in df.keys()):
    df = df[~np.isnan(df.heii_4685_flux_port)]
    heii = df['heii_4685_flux_port']
    heii_err = df['heii_4685_flux_port_err']

else:
    df = df[~np.isnan(df.Flux_HeII_4685)]
    heii = df['Flux_HeII_4685']
    heii_err = df['Flux_HeII_4685_Err']

#define alternate catalog names
if 'name' in df.keys():
    df['NAME'] = df['name']

name = df['NAME']
if eco: 
    resname = df['resname'] #for eco
    resname = resname != 'notinresolve'
if resolve:
    econame = df['econame'] #for resolve
    econame = econame != 'notineco'

#define demarcation function: log_NII_HA vs. log_OIII_HB
def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
    return 1.19 + (0.61 / (log_NII_HA - 0.47))
def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
    return 1.30 + (0.72 / (log_SII_HA - 0.32))
def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
    return 0.76 + 1.89*log_SII_HA
def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
    return 1.33 + (0.73 / (log_OI_HA + 0.59))
def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
    return 1.3 + 1.18*log_OI_HA
def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
    return -0.59

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

#create line ratios [NII]/H-alpha and [OIII]/H-beta
nii = df['nii_6584_flux']
nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4

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
sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)


#need to check Kewley paper and figure out if ratio used in nii_sum applies to sii_sum as well

#Filter Data: all non-negative SEL fluxes and errors; Hbeta >3sigma
gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & (oi > 0) &
            (sii_sum > 0) & (h_beta > 0) & (h_beta > 3*h_beta_err) &
            (h_alpha_err > 0) & (nii_sum_err > 0) & (oiii_err > 0) & 
            (oi_err > 0) & (sii_sum_err > 0))

he2data = (heii/heii_err >=5) & (heii_err > 0)

if he2_flag:
    data = gooddata & he2data
else:
    data = gooddata #use ALL galaxy data within catalog

#print data type being used
print ''
print 'ALL DATA'


#print total points shared with alternate catalog

if full: 
    sel = (np.where(data)[0]) #for eco
elif eco: 
    sel = (np.where(data & resname)[0]) #for eco
elif resolve:
    sel = (np.where(data & econame)[0]) #for resolve
print ''
print 'TOTAL DATA WITH ALTERNATE CATALOG NAME: ', len(sel)

nii = nii[data]
nii_sum = nii_sum[data]
oiii = oiii[data]
oi = oi[data]
sii_sum = sii_sum[data]
h_beta = h_beta[data]
h_alpha = h_alpha[data]
heii = heii[data] # 3-sigma cut for HeII selection
subsetname = name[data]

#length of data to be used for debugging
datalen = np.sum(data)

# data ratios
#n2ha = np.log10(nii_sum/h_alpha)
o3hb = np.log10(oiii/h_beta) # always the y-axis
o1ha = np.log10(oi/h_alpha)
s2ha = np.log10(sii_sum/h_alpha)
he2hb = np.log10(heii/h_beta)
n2ha = np.log10(nii/h_alpha)
#Below are the selectors for the data to distinguish btwn: Seyferts, Composites,
#and AGN's based on the flux ratio diagnostic as understood via Kewley 2006.

#NII plot selectors
compsel1 = (o3hb >= n2hacompmin(n2ha)) & (o3hb <= n2hamain(n2ha))
sfsel1 = (o3hb < n2hacompmin(n2ha)) & (n2ha < 0.) & ~(o3hb > n2hamain(n2ha)) #~(o3hb > n2hamain(n2ha)) & ~compsel1
agnsel1= (o3hb > n2hamain(n2ha))

#SII plot selectors
sfsel2 = (o3hb <= s2hamain(s2ha)) & ~compsel1
seyfsel2 = ((o3hb > s2hamain(s2ha)) & (o3hb >= s2halinseyf(s2ha)))
linersel2 = ((o3hb > s2hamain(s2ha)) & (o3hb < s2halinseyf(s2ha)))
agnsel2 = (o3hb > s2hamain(s2ha)) & ~compsel1

#OI plot selectors
sfsel3 = (o3hb <= o1hamain(o1ha)) & (o1ha < -0.7) & ~compsel1
seyfsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb >= o1halinseyf(o1ha))
linersel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & (o3hb < o1halinseyf(o1ha))
agnsel3 = ((o3hb > o1hamain(o1ha)) | (o1ha > -0.7)) & ~compsel1

#HEII plot selectors
if he2_flag:
    refn2ha = np.linspace(-3.0, 0.35)
    main = he2hbmain(refn2ha[refn2ha < -0.15])
    limit = he2hblimit(refn2ha[refn2ha < -0.15])
#    sfsel4 = (he2hb <= main) & (he2hb <= limit)
#    agnsel4 = (he2hb > main) & (he2hb > limit)
    sfsel4 = (he2hb <= he2hbmain(n2ha)) & \
                (he2hb <= he2hblimit(n2ha))
    agnsel4 = ((he2hb > he2hbmain(n2ha)) | (n2ha >= -0.2)) & \
                ((he2hb > he2hblimit(n2ha)) | (n2ha >= 0))
#for BPT comparison

print ''
print 'CHECKING FOR OVERPLOTTING:'
overlap = (np.where(agnsel1 & compsel1)[0])
print 'NII - AGN/ COMPOSITE:', len(overlap)
overlap = (np.where(compsel1 & sfsel1)[0])
print 'NII - COMPOSITE/ STAR FORMING:', len(overlap)
overlap = (np.where(agnsel1 & sfsel1)[0])
print 'NII - AGN/ STAR FORMING:', len(overlap)
overlap = (np.where(agnsel2 & sfsel2)[0])
print 'SII - AGN/ STAR FORMING:', len(overlap)
overlap = (np.where(seyfsel2 & sfsel2)[0])
print 'SII - SEYFERT/ STAR FORMING:', len(overlap)
overlap = (np.where(seyfsel2 & linersel2)[0])
print 'SII - SEYFERT/ LINER:', len(overlap)
overlap = (np.where(linersel2 & sfsel2)[0])
print 'SII - LINER/ STAR FORMING:', len(overlap)
overlap = (np.where(agnsel3 & sfsel3)[0])
print 'OI - AGN/ STAR FORMING:', len(overlap)
overlap = (np.where(seyfsel3 & sfsel3)[0])
print 'OI - SEYFERT/ STAR FORMING:', len(overlap)
overlap = (np.where(seyfsel3 & linersel3)[0])
print 'OI - SEYFERT/ LINER:', len(overlap)
overlap = (np.where(linersel3 & sfsel3)[0])
print 'OI - LINER/ STAR FORMING:', len(overlap)
if he2_flag:
    overlap = (np.where(agnsel4 & sfsel4)[0])
    print 'HeII - AGN / STAR FORMING:', len(overlap)

####stop here for HeII additions (i.e. NOT included in calculations) ####

#REFERENCE for cumulative plot selectors
seyfselr = seyfsel2 & seyfsel3
linerselr = linersel2 & linersel3

#cumulative plot selectors
if he2_flag:
    sfsel = sfsel1 & sfsel2 & sfsel3 & sfsel4 #definite star forming galaxies
else:
    sfsel = sfsel1 & sfsel2 & sfsel3 
compsel = compsel1  #composite galaxies
seyfsel = agnsel1 & seyfselr #Seyfert AGN galaxies
linersel = agnsel1 & linerselr #LINER AGN galaxies
ambigsel1 = sfsel1 & (agnsel2 | agnsel3) #SF in first plot, AGN in subsequent plot
ambigsel2 = agnsel1 & (sfsel2 | sfsel3) #AGN in first plot, SF in subsequent plot
ambagnsel = agnsel1 & ~seyfselr & ~linerselr & ~(sfsel2 | sfsel3) #Ambiguous AGN

sftoagn1 = sfsel1 & agnsel2
sftoagn2 = sfsel1 & agnsel3

emlineclass = sfsel ^ compsel ^ seyfsel ^ linersel ^ ambigsel1 ^ ambigsel2 ^ ambagnsel
defagn = seyfsel ^ linersel ^ ambagnsel
#print np.sum(emlineclass)

if save:
    if not he2_flag:    
        dfout = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                          'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                          'sftoagn':ambigsel1, 'agntosf':ambigsel2, 'defagn': defagn,
                          'sftoagn1':sftoagn1, 'sftoagn2': sftoagn2})
        dfout.to_csv(outputfile,index=False)
    
    else:
        dfout = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                          'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                          'sftoagn':ambigsel1, 'agntosf':ambigsel2, 'defagn': defagn,
                          'heiisel':agnsel4})
        dfout.to_csv(outputfile ,index=False)

#create alternate catalog name-based agn selector, print len
if full: 
    agn = (np.where(agnsel1)[0]) #for eco
if eco: 
    agn = (np.where(agnsel1 & resname)[0]) #for eco
if resolve: 
    agn = (np.where(agnsel1 & econame)[0]) #for resolve
print ''
print 'AGN WITH ALTERNATE CATALOG NAME: ', len(agn)


#recheck for overplotting
print ''
print 'RE-CHECKING FOR OVERPLOTTING:'
overlap = (np.where(sfsel & ambagnsel)[0])
print 'DEFINITE STAR FORMING/ AMBIGUOUS AGN:', len(overlap)
overlap = (np.where(sfsel & seyfsel)[0])
print 'DEFINITE STAR FORMING/ SEYFERT:', len(overlap)
overlap = (np.where(sfsel & linersel)[0])
print 'DEFINITE STAR FORMING/ LINER:', len(overlap)
overlap = (np.where(sfsel & ambigsel1)[0])
print 'DEFINITE STAR FORMING/ AMBIGUOUS1:', len(overlap)
overlap = (np.where(sfsel & ambigsel2)[0])
print 'DEFINITE STAR FORMING/ AMBIGUOUS2:', len(overlap)
overlap = (np.where(ambagnsel & seyfsel)[0])
print 'AMBIGUOUS AGN/ SEYFERT:', len(overlap)
overlap = (np.where(ambagnsel & linersel)[0])
print 'AMBIGUOUS AGN/ LINER:', len(overlap)
overlap = (np.where(ambagnsel & ambigsel1)[0])
print 'AMBIGUOUS AGN/ AMBIGUOUS1', len(overlap)
overlap = (np.where(ambagnsel & ambigsel2)[0])
print 'AMBIGUOUS AGN/ AMBIGUOUS2', len(overlap)
overlap = (np.where(seyfsel & linersel)[0])
print 'SEYFERT/ LINER:', len(overlap)
overlap = (np.where(seyfsel & ambigsel1)[0])
print 'SEYFERT/ AMBIGUOUS1:', len(overlap)
overlap = (np.where(seyfsel & ambigsel2)[0])
print 'SEYFERT/ AMBIGUOUS2:', len(overlap)
overlap = (np.where(compsel1 & ambigsel1)[0])
print 'COMPOSITE/ AMBIGUOUS1:', len(overlap)
overlap = (np.where(compsel1 & ambigsel2)[0])
print 'COMPOSITE/ AMBIGUOUS2:', len(overlap)
overlap = (np.where(ambigsel1 & ambigsel2)[0])
print 'AMBIGUOUS1/ AMBIGUOUS2:', len(overlap)


#checking that plotted points are within the total data range
print ''
sfselpts = (len(np.where(sfsel)[0]))
seyfselpts = (len(np.where(seyfsel)[0]))
linerselpts = (len(np.where(linersel)[0]))
compselpts = (len(np.where(compsel)[0]))
agnselpts = (len(np.where(ambagnsel)[0]))
ambigsel1pts = (len(np.where(ambigsel1)[0]))
ambigsel2pts = (len(np.where(ambigsel2)[0]))
if he2_flag:
    heiiselpts = (len(np.where(agnsel4)[0]))
    totalselpts = sfselpts+seyfselpts+linerselpts+compselpts+agnselpts+ambigsel1pts+ambigsel2pts+heiiselpts
else:
    totalselpts = sfselpts+seyfselpts+linerselpts+compselpts+agnselpts+ambigsel1pts+ambigsel2pts

print ("DATA POINTS: "),datalen
print ("TOTAL PLOTTED POINTS: "), totalselpts
print ("TOTAL PLOTTED POINTS OMITTED: "), datalen-totalselpts
print "* IF NUMBERS ABOVE ARE AT ALL NEGATIVE THEN THERE IS OVERPLOTTING"
print ("Definite Star Forming: "),sfselpts
print ("Composite: "),compselpts
print ("SF --> AGN: "), ambigsel1pts
print ("AGN --> SF: "), ambigsel2pts
print ("Ambiguous AGN: "),agnselpts
print ("Seyfert: "),seyfselpts
print ("LINER: "),linerselpts
print ("TOTAL KNOWN AGN: "),linerselpts+seyfselpts+agnselpts
print ("POSSIBLE TOTAL AGN: "),linerselpts+seyfselpts+agnselpts+ambigsel1pts+ambigsel2pts

#display percent of each category
print ''
print "PERCENT OF GALAXIES IN EACH CATEGORY"
sfpercent = float(sfselpts)/float(datalen)*100
seyfpercent = float(seyfselpts)/float(datalen)*100
linerpercent = float(linerselpts)/float(datalen)*100
comppercent = float(compselpts)/float(datalen)*100
agnpercent = float(agnselpts)/float(datalen)*100
ambig1percent = float(ambigsel1pts)/float(datalen)*100
ambig2percent = float(ambigsel2pts)/float(datalen)*100
print ("Definite Star Forming:"), round(sfpercent, 2), ("%")
print("Composite: "),round(comppercent, 2), ("%")
print("SF --> AGN: "),round(ambig1percent, 2), ("%")
print("AGN --> SF: "),round(ambig2percent, 2), ("%")
print("Ambiguous AGN: "),round(agnpercent, 2), ("%")
print("Seyfert: "),round(seyfpercent, 2), ("%")
print("LINER: "),round(linerpercent, 2), ("%")
print("TOTAL KNOWN AGN: "),round(linerpercent+seyfpercent+agnpercent, 2), ("%")
print("POSSIBLE TOTAL AGN: "),round(linerpercent+seyfpercent+agnpercent+ambig1percent+ambig2percent, 2), ("%")
print ("Percent Omitted: "), round((100-(sfpercent+seyfpercent+linerpercent+comppercent+agnpercent+ambig1percent+ambig2percent)), 2), ("%")
print ''

###PLOTS###
'''fig = plt.figure(0)
fig1 = plt.figure(1)
fig2 = plt.figure(2)
fig3 = plt.figure(3)
plt.clf()
'''
#reference points in x-direction for demarcation lines on plots
refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, -0.59)
refsiiha = np.linspace(-2, 0.3)

#NII/OIII plot
import matplotlib
def makeColours( vals ):
    colours = np.zeros( (len(vals),3) )
    norm = matplotlib.colors.Normalize( vmin=vals.min(), vmax=vals.max() )

    #Can put any colormap you like here.
    colours = [matplotlib.cm.ScalarMappable(norm=norm, cmap='Blues_r').to_rgba(val) for val in vals]

    return colours
if resolve == 1:    
    fig = plt.figure('Full Scatter Plot')
    ax1 = fig.add_subplot(221)
    main1, = ax1.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'ke01 Theoretical Maximum Starburst Line')
    composite, = ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 'k-.', label = 'Ka03 Composite Line')
    ax1.set_xlim(-2,1)
    ax1.set_ylim(-1.25,1.5)
    sfsel1, = ax1.plot(n2ha[sfsel], o3hb[sfsel], 'k.', alpha = 0.3, markersize = 5)#, label = 'Definite Star Forming')
    ambig1data1, = ax1.plot(n2ha[ambigsel1], o3hb[ambigsel1],'rs', markersize = 8, mew = 0)#, label = 'SF -> AGN (SF-to-AGN)')
    compdata1, = ax1.plot(n2ha[compsel], o3hb[compsel], 'bs', markersize = 8, mew = 0)#, label = 'Composite')
    seyfsel1, = ax1.plot(n2ha[seyfsel], o3hb[seyfsel], 'mo', markersize = 8, mew = 0)#, label = 'Seyfert')
    liner1, = ax1.plot(n2ha[linersel], o3hb[linersel], 'co', markersize = 8, mew = 0)#, label = 'LINER')
    ambig1, = ax1.plot(n2ha[ambagnsel], o3hb[ambagnsel], 'yo', markersize = 8, mew = 0)#, label = 'Ambiguous AGN')
    ambig2data1, = ax1.plot(n2ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 12, mew = 0)#, label = 'AGN -> SF')
    ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
    ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
    if he2_flag:
        agndata4, = ax1.plot(n2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
    
    #SII/OIII plot
    #plt.figure('SII Scatter Plot')
    ax2 = plt.subplot(222)
    main2, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
    liner, = ax2.plot(refsiiha, s2halinseyf(refsiiha), 'k--', label = 'Liner/Seyfert Division')
    ax2.set_xlim(-2, 0.5)
    ax2.set_ylim(-1.25,1.5)
    ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
    ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
    sfdata2, = ax2.plot(s2ha[sfsel], o3hb[sfsel], 'k.', markersize = 5, alpha = 0.5, label = 'SF')
    ambig1data2, = ax2.plot(s2ha[ambigsel1], o3hb[ambigsel1], 'rs', markersize = 8, mew = 0, label = 'SF-to-AGN')
    ambig2data2, = ax2.plot(s2ha[ambigsel2], o3hb[ambigsel2], 'g^', markersize = 12, mew = 0, label = 'AGN-to-SF ')
    seyfdata2, = ax2.plot(s2ha[seyfsel], o3hb[seyfsel], 'mo', markersize = 8, mew = 0, label = 'Seyfert')
    linerdata2, = ax2.plot(s2ha[linersel], o3hb[linersel],'co', markersize = 8, mew = 0, label = 'LINER')
    agndata2, = ax2.plot(s2ha[ambagnsel], o3hb[ambagnsel],'yo', markersize = 8, mew = 0, label = 'Ambiguous AGN')
    compdata2, = ax2.plot(s2ha[compsel], o3hb[compsel], 'bs',markersize = 8, mew = 0, label = 'Composite')
    if he2_flag:
        agndata4, = ax2.plot(s2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
    
    #OI/OIII plot
    #plt.figure('OI Scatter Plot')
    ax3 = plt.subplot(223)
    main3, = ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]), 'k', label = 'Ke01 Theoretical Maximum Starburst Line')
    comp3, = ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]), 'k-.', label = 'Ka03 Composite Line')
    #main3, = ax3.plot(refoiha, o1hamain(refoiha), 'k', label = 'Main Line')
    #liner, = ax3.plot(refoiha, o1halinseyf(refoiha), 'k--', label = 'Composite Line (Plot 1)')
    liner2, = ax3.plot(refoiha, o1halinseyf(refoiha), 'k--', label = 'Ke06 Liner/Seyfert Division Line')
    ax3.set_xlim(-2.5, 0)
    ax3.set_ylim(-1.25,1.5)
    ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
    ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
    sfdata3, = ax3.plot(o1ha[sfsel], o3hb[sfsel], 'k.', alpha = 0.3, markersize = 5, mew = 0, label = 'SF')
    ambig1data3, = ax3.plot(o1ha[ambigsel1], o3hb[ambigsel1],'rs', markersize = 8, mew = 0, label = 'SF-to-AGN')#'SF-to-AGN')
    ambig2data3, = ax3.plot(o1ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 10, mew = 0, label = 'AGN-to-SF')
    seyfdata3, = ax3.plot(o1ha[seyfsel], o3hb[seyfsel], 'mo', markersize = 8, mew = 0, label = 'Seyfert')
    linerdata3, = ax3.plot(o1ha[linersel], o3hb[linersel],'co', markersize = 8, mew = 0, label = 'LINER')
    agndata3, = ax3.plot(o1ha[ambagnsel], o3hb[ambagnsel],'yo', markersize = 8, mew = 0, label = 'Ambiguous AGN')
    #agnsel3, = ax3.plot(o1ha[defagn], o3hb[defagn],'rs', alpha = 0.5,markersize = 8, mew = 0, label = 'Definite AGN')
    compdata3, = ax3.plot(o1ha[compsel], o3hb[compsel], 'bs', markersize = 8, mew = 0, label = 'Composite')
    if he2_flag:
        agndata4, = ax3.plot(o1ha[agnsel4], o3hb[agnsel4],'ks',  markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
    plt.legend(bbox_to_anchor=(1.25, 1), loc=2, borderaxespad=0., numpoints = 1)
    #plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)
    
    ##N2/HeII plot
    if he2_flag:
    
        plt.figure()
        ax4 = plt.subplot(111)
        main4, = ax4.plot(refn2ha[refn2ha < -0.15], he2hbmain(refn2ha[refn2ha < -0.15]), 'k',  label = 'Main Line')
        liner4, = ax4.plot(refn2ha[refn2ha < 0.1], he2hblimit(refn2ha[refn2ha < 0.1]), 'k--', label = 'Composite Line (Plot 1)')
#        main4, = ax4.plot(refn2ha, he2hbmain(refn2ha), 'k',  label = 'Main Line')
#        liner4, = ax4.plot(refn2ha, he2hblimit(refn2ha), 'k--', label = 'Composite Line (Plot 1)')
        ax4.set_xlim(-3., 1.)
        ax4.set_ylim(-3, 1.)
        ax4.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
        ax4.set_ylabel(r"$\rm \log([HeII]/H\beta)$", fontsize = 22)
        #data4, = ax4.plot(n2ha, he2hb, 'k.')
        sfdata4, = ax4.plot(n2ha[sfsel4], he2hb[sfsel4], 'ko',  markersize = 10, alpha = 0.5, label = 'SF')
        agndata4, = ax4.plot(n2ha[agnsel4], he2hb[agnsel4],'ks',  markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
        #alldata, = ax4.plot(n2ha, he2hb, 'kx',  markersize = 10, alpha = 0.5, label = 'SF')
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
        plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)
'''
xmin = refn2ha.min(); xmax = refn2ha.max()
ymin = -1.25; ymax = 1.5
nbins = 50
fig = plt.figure('ECO Final Plot', figsize = (5,10))
ax1 = fig.add_subplot(311)

definite = np.column_stack((n2ha[defagn|sfsel], o3hb[defagn|sfsel]))
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((n2ha[defagn], o3hb[defagn]))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 6, colors='k')
ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), shading='gouraud', cmap=plt.cm.gray_r)
main1, = ax1.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'Main Line')
composite, = ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 'k-.', label = 'Composite Line')
ax1.set_xlim(-2,0.5)
ax1.set_ylim(-1.25,1.5)
ambig1data1, = ax1.plot(n2ha[ambigsel1], o3hb[ambigsel1],'rs', alpha = 0.5, 
                    markersize = 8, mew = 0, label = 'SF-to-AGN')
#sftoagn1data, = ax1.plot(n2ha[sftoagn1], o3hb[sftoagn1],'cs', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#sftoagn2data, = ax1.plot(n2ha[sftoagn2], o3hb[sftoagn2],'ms', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
compdata1, = ax1.plot(n2ha[compsel], o3hb[compsel], 'bs', alpha = 0.5, markersize = 8, mew = 0, label = 'Composite')
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ambig2data1, = ax1.plot(n2ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 10,
                        mew = 1, mec = 'y', label = 'AGN-to-SF')
if he2_flag:
    agndata4, = ax1.plot(n2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
#plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)

xmin = refsiiha.min(); xmax = refsiiha.max()
ymin = -1.25; ymax = 1.5
nbins = 50
ax1 = fig.add_subplot(312)

definite = np.column_stack((s2ha[defagn|sfsel], o3hb[defagn|sfsel]))
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((s2ha[defagn], o3hb[defagn]))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 6, colors='k')
ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), shading='gouraud', cmap=plt.cm.gray_r)
main1, = ax1.plot(refsiiha, s2hamain(refsiiha), 'k', label = 'Main Line')
ax1.set_xlim(-1.5,xmax)
ax1.set_ylim(-1.25,1.5)
ambig1data1, = ax1.plot(s2ha[ambigsel1], o3hb[ambigsel1],'rs', alpha = 0.5, 
                    markersize = 8, mew = 0, label = 'SF-to-AGN')
#sftoagn1data, = ax1.plot(s2ha[sftoagn1], o3hb[sftoagn1],'cs', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#sftoagn2data, = ax1.plot(s2ha[sftoagn2], o3hb[sftoagn2],'ms', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
compdata1, = ax1.plot(s2ha[compsel], o3hb[compsel], 'bs', alpha = 0.5, markersize = 8, mew = 0, label = 'Composite')
ax1.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ambig2data1, = ax1.plot(s2ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 10,
                        mew = 1, mec = 'y', label = 'AGN-to-SF')
if he2_flag:
    agndata4, = ax1.plot(s2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
#plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)


#OI Plot
xmin = refoiha.min(); xmax = 0#refoiha.max()
ymin = -1.25; ymax = 1.5
nbins = 50

ax1 = fig.add_subplot(313)

definite = np.column_stack((o1ha[defagn|sfsel], o3hb[defagn|sfsel]))
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((o1ha[defagn], o3hb[defagn]))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 6, 
            colors='k')
ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
               shading='gouraud', cmap=plt.cm.gray_r)
composite_fake, = ax1.plot(refoiha, o1hamain(refoiha), 'k-.', 
                           label = 'Ka03 Composite Line')
main1, = ax1.plot(refoiha, o1hamain(refoiha), 'k', 
                  label = 'Ke01 Maximum Starburst Line')
ax1.set_xlim(xmin,0)
ax1.set_ylim(-1.25,1.5)
compdata1, = ax1.plot(o1ha[compsel], o3hb[compsel], 'bs', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Composite')
#sftoagn1data, = ax1.plot(o1ha[sftoagn1], o3hb[sftoagn1],'cs', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#sftoagn2data, = ax1.plot(o1ha[sftoagn2], o3hb[sftoagn2],'ms', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
ambig1data1, = ax1.plot(o1ha[ambigsel1], o3hb[ambigsel1],'rs', alpha = 0.5, 
                        markersize = 8, mew = 0, label = 'SF-to-AGN')
ax1.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ambig2data1, = ax1.plot(o1ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 10,
                        mew = 1, mec = 'y', label = 'AGN-to-SF')
if he2_flag:
    agndata4, = ax1.plot(o1ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8,
                         mfc ='none', mew = 2, label = 'HeII-Selected AGN')
ax1.legend(borderaxespad=0, loc = 3, numpoints = 1, prop={'size': 6})
#plt.subplots_adjust(left=0.11, bottom=0.11, right=0.11, top=0.11, hspace=0.5)

'''
def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
  	new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap
#metal_colors_map = truncate_colormap(cm.Blues, n = metals)
#u_colors_map = truncate_colormap(cm.Reds, n = ionps)

sf_colors_map = truncate_colormap(cm.gray_r)


#xray = pd.read_csv('../SDSS_spectra/xray/ECO+RESOLVE_xray.csv')
#xray.index = xray.name
xmin = refn2ha.min(); xmax = refn2ha.max()
ymin = -1.25; ymax = 1.5
nbins = 50
fig = plt.figure('NII Plot')
ax1 = fig.add_subplot(111)

definite = np.column_stack((n2ha[defagn|sfsel], o3hb[defagn|sfsel]))
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((n2ha[defagn], o3hb[defagn]))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
               shading='gouraud', cmap=sf_colors_map) #plt.cm.gray_r)
main1, = ax1.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'Main Line')
composite, = ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 
'k-.', label = 'Composite Line')
ax1.set_xlim(-1.7,0.4)
ax1.set_ylim(-1.0,1.0)
ambig1data1, = ax1.plot(n2ha[ambigsel1], o3hb[ambigsel1], 'bs',
                        alpha = 0.5, markersize = 8, mew = 0, label = 'SF-to-AGN')
#sftoagn1data, = ax1.plot(n2ha[sftoagn1], o3hb[sftoagn1],'rs', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#sftoagn2data, = ax1.plot(n2ha[sftoagn2], o3hb[sftoagn2],'m*', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#xraydata, = ax1.plot(n2ha[xray.index.values], o3hb[xray.index.values],'ks', alpha = 0.5, markersize = 8, mew = 0, label = 'X-Ray')
compdata1, = ax1.plot(n2ha[compsel], o3hb[compsel], 's', color = 'purple', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Composite')
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ambig2data1, = ax1.plot(n2ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 10,
                        mew = 1, mec = 'y', label = 'AGN-to-SF')
ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 3, 
            colors='k', corner_mask = True, extend = 'both')
if he2_flag:
    agndata4, = ax1.plot(n2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8,
                         mfc ='none', mew = 2, label = 'HeII-Selected AGN')
plt.legend(loc=3, numpoints = 1, fontsize = 14)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
#plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)
plt.show()
xmin = refsiiha.min(); xmax = refsiiha.max()
ymin = -1.25; ymax = 1.5
nbins = 50

fig = plt.figure('SII Plot')
ax1 = fig.add_subplot(111)

definite = np.column_stack((s2ha[defagn|sfsel], o3hb[defagn|sfsel]))
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((s2ha[defagn], o3hb[defagn]))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 3,
            colors='k')
ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), shading='gouraud', cmap=plt.cm.gray_r)
main1, = ax1.plot(refsiiha, s2hamain(refsiiha), 'k', label = 'Main Line')
ax1.set_xlim(-1.2,0.2)
ax1.set_ylim(-1,1)
compdata1, = ax1.plot(s2ha[compsel], o3hb[compsel], 'bs', alpha = 0.5, markersize = 8, mew = 0, label = 'Composite')
ambig1data1, = ax1.plot(s2ha[ambigsel1], o3hb[ambigsel1],'rs', alpha = 0.5, 
                    markersize = 8, mew = 0, label = 'SF-to-AGN')
#sftoagn1data, = ax1.plot(s2ha[sftoagn1], o3hb[sftoagn1],'rs', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#sftoagn2data, = ax1.plot(s2ha[sftoagn2], o3hb[sftoagn2],'m*', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#xraydata, = ax1.plot(s2ha[xray.index.values], o3hb[xray.index.values],'ks', alpha = 0.5, markersize = 8, mew = 0, label = 'X-Ray')
ax1.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ambig2data1, = ax1.plot(s2ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 10,
                        mew = 1, mec = 'y',label = 'AGN-to-SF')
if he2_flag:
    agndata4, = ax1.plot(s2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, 
                         mfc ='none', mew = 2, label = 'HeII-Selected AGN')
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
#plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)
plt.show()

#OI Plot
xmin = refoiha.min(); xmax = 0#refoiha.max()
ymin = -1.25; ymax = 1.5
nbins = 50
fig = plt.figure('OI Plot')
ax1 = fig.add_subplot(111)

definite = np.column_stack((o1ha[defagn|sfsel], o3hb[defagn|sfsel]))
xgrid, ygrid = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
k2 = kde.gaussian_kde(definite.T)
definite_z = k2(np.vstack([xgrid.flatten(), ygrid.flatten()]))
agn_contour = np.column_stack((o1ha[defagn], o3hb[defagn]))
xmin_agn = agn_contour[:,0].min(); xmax_agn = agn_contour[:,0].max()
ymin_agn = agn_contour[:,1].min(); ymax_agn = agn_contour[:,1].max()
xgrid_agn, ygrid_agn = np.mgrid[xmin_agn:xmax_agn:nbins*1j, 
                                ymin_agn:ymax_agn:nbins*1j]
k = kde.gaussian_kde(agn_contour.T)
agn_contour_z = k(np.vstack([xgrid_agn.flatten(), ygrid_agn.flatten()]))
ax1.contour(xgrid_agn, ygrid_agn, agn_contour_z.reshape(xgrid_agn.shape), 3, 
            colors='k')
ax1.pcolormesh(xgrid, ygrid, definite_z.reshape(xgrid.shape), 
               shading='gouraud', cmap=plt.cm.gray_r)
composite_fake, = ax1.plot(refoiha, o1hamain(refoiha), 'k-.', 
                           label = 'Ka03 Composite Line')
main1, = ax1.plot(refoiha, o1hamain(refoiha), 'k', 
                  label = 'Ke01 Maximum Starburst Line')
ax1.set_xlim(-2.1,-0.2)
ax1.set_ylim(-1,1)
compdata1, = ax1.plot(o1ha[compsel], o3hb[compsel], 'bs', alpha = 0.5, 
                      markersize = 8, mew = 0, label = 'Composite')
ambig1data1, = ax1.plot(o1ha[ambigsel1], o3hb[ambigsel1],'rs', alpha = 0.5,
                        markersize = 8, mew = 0, label = 'SF-to-AGN')
#sftoagn1data, = ax1.plot(o1ha[sftoagn1], o3hb[sftoagn1],'rs', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#sftoagn2data, = ax1.plot(o1ha[sftoagn2], o3hb[sftoagn2],'m*', alpha = 0.5, markersize = 8, mew = 0, label = 'SF -> AGN (SF-to-AGN)')
#xraydata, = ax1.plot(o1ha[xray.index.values], o3hb[xray.index.values],'ks', alpha = 0.5, markersize = 8, mew = 0, label = 'X-Ray')
ax1.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ambig2data1, = ax1.plot(o1ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 10,
                        mew = 1, mec = 'y', label = 'AGN-to-SF')
if he2_flag:
    agndata4, = ax1.plot(o1ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, 
                         mfc ='none', mew = 2, label = 'HeII-Selected AGN')
#plt.legend(loc=2, borderaxespad=0., numpoints = 1)
#plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, 
#                    hspace=None)
plt.show()

'''
flags = pd.read_csv(outputfile)
flags.index = flags.galname
xray_flags = flags.loc[xray.index.values]
keys = np.array(xray_flags.keys())
print 'X-Ray Detected AGN (Total: 98)'
print 'Definite SF : ',np.sum(xray_flags['defstarform'])
print 'Definite AGN: ', np.sum(xray_flags['defagn'])
print 'Composite : ',np.sum(xray_flags['composite']) 
print 'SF-to-AGN : ',np.sum(xray_flags['sftoagn']) 
print 'AGN-to-SF : ',np.sum(xray_flags['agntosf']) 
#for i in range(len(xray_flags)): 
#    print xray_flags.galname[i],keys[np.where(xray_flags.iloc[i])[0]]
'''