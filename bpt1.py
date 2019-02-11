# -*- coding: utf-8 -*-
"""
Created on Thu Dec 06 03:12:35 2018

@author: mugdhapolimera
"""

#This program makes a Line-Ratio diagram (also known as a BPT plot or Kewley diagram)
#with labels using data from the RESOLVE survey to classify galaxies as LINERs,
#Seyferts, Composites, or AGNs on the basis of their flux ratio for distinction.

#Original code from Ashley Bittner 08/03/2017
#Edited version from Margie Bruff 01/07/2018
#Updated by Carlynn Ferguson 03/30/2018

#suggested use of python debugger to understand the code more thoroughly
#pdb.set_trace()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#plt.ion()
#import pdb
#os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/BPT/')
#display catalog being used
print ''
#print 'ECO RESULTS'
print 'RESOLVE RESULTS'

#read in data
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_full.pkl'
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpt1_filter.pkl'
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_blend_dext.pkl'
df = pd.read_pickle(inputfile) #ECO catalog
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_all_dext.fits'
#inputfile = 'RESOLVE_SDSS_dext.fits'
#dat = Table.read(inputfile, format='fits')
#df = dat.to_pandas()
he2_flag = 0
save = 1

if ('heii_4685_flux_port_ext' in df.keys()):
    df = df[~np.isnan(df.heii_4685_flux_port_ext)]
    heii = df['heii_4685_flux_port']
    heii_err = df['heii_4685_flux_port_err']

else:
    df = df[~np.isnan(df.Flux_HeII_4685)]
    heii = df['Flux_HeII_4685']
    heii_err = df['Flux_HeII_4685_Err']

#define alternate catalog names
name = df['name']
resname = df['resname'] #for eco
resname = resname != 'notinresolve'

#econame = df['econame']
#econame = df['econame'] #for resolve
#econame = econame != 'notineco'

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
def he2hblimit(log_NII_HA):
        return -1.07+1.0/(8.92*log_NII_HA-0.95)


#create line ratios [NII]/H-alpha and [OIII]/H-beta
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
sii_sum_err = np.sqrt(df['sii_6717_flux']**2 + df['sii_6731_flux']**2)


#need to check Kewley paper and figure out if ratio used in nii_sum applies to sii_sum as well

#Filter Data: all non-negative SEL fluxes and errors; Hbeta >3sigma
gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & 
            (h_beta > 0) & (h_beta > 3*h_beta_err) &
            (h_alpha_err > 0) & (nii_sum_err > 0) & (oiii_err > 0))# & 
#            (oi > 0) & (sii_sum > 0) & (oi_err > 0) & (sii_sum_err > 0))

he2data = (heii/heii_err >=3) & (heii_err > 0)

if he2_flag:
    data = gooddata & he2data
else:
    data = gooddata #use ALL galaxy data within catalog

#results = pd.read_csv("C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/results_bpass_agn_snr/RESOLVE_param_estimates.csv")
#agn_index = (results['Parameter'] == 'AGNFRAC')
#results_agn = results[agn_index]
#results_agn.index = results_agn['Galaxy Name']
#Z_index = (results['Parameter'] == 'LOGZ')
#results_Z = results[Z_index]
#results_Z.index = results_agn['Galaxy Name']

#agn_sel = (results_agn['Estimate'] > 0) & (results_Z['Estimate'] != min(results_Z['Estimate']))
data = data #& agn_sel
print len(np.where(data)[0])

#print data type being used
print ''
print 'ALL DATA'
#print 'DWARF ONLY DATA'
#print 'ALL BLUE E/S0 DATA'
#print 'DWARF ONLY BLUE E/S0 DATA'
#print 'BLUE NUGGET CANDIDATES'


#print total points shared with alternate catalog
#sel = (np.where(data & resname)[0]) #for eco
sel = (np.where(data))# & econame)[0]) #for resolve
print ''
print 'TOTAL DATA WITH ALTERNATE CATALOG NAME: ', len(sel)

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
n2ha = np.log10(nii_sum/h_alpha)
o3hb = np.log10(oiii/h_beta) # always the y-axis
o1ha = np.log10(oi/h_alpha)
s2ha = np.log10(sii_sum/h_alpha)
he2hb = np.log10(heii/h_beta)

#Below are the selectors for the data to distinguish btwn: Seyferts, Composites,
#and AGN's based on the flux ratio diagnostic as understood via Kewley 2006.

#need to figure out why plot does not 'go back up' and consider adding to selectors
#to prohibit the data once the line 'goes back up'

#NII plot selectors
compsel1 = (o3hb >= n2hacompmin(n2ha)) & (o3hb <= n2hamain(n2ha))
sfsel1 = (o3hb < n2hacompmin(n2ha)) & (n2ha < 0.) & ~(o3hb > n2hamain(n2ha)) #~(o3hb > n2hamain(n2ha)) & ~compsel1
agnsel1= (o3hb >= n2hamain(n2ha))


#HEII plot selectors
if he2_flag:
    sfsel4 = (he2hb <= he2hbmain(n2ha)) & (he2hb <= he2hblimit(n2ha))
    agnsel4 = (he2hb > he2hbmain(n2ha)) 
#for BPT comparison

print ''
print 'CHECKING FOR OVERPLOTTING:'
overlap = (np.where(agnsel1 & compsel1)[0])
print 'NII - AGN/ COMPOSITE:', len(overlap)
overlap = (np.where(compsel1 & sfsel1)[0])
print 'NII - COMPOSITE/ STAR FORMING:', len(overlap)
overlap = (np.where(agnsel1 & sfsel1)[0])
print 'NII - AGN/ STAR FORMING:', len(overlap)

####stop here for HeII additions (i.e. NOT included in calculations) ####


#cumulative plot selectors
if he2_flag:
    sfsel = sfsel1 #definite star forming galaxies
else:
    sfsel = sfsel1 
compsel = compsel1  #composite galaxies
agnsel = agnsel1  #Seyfert AGN galaxies

emlineclass = sfsel ^ compsel ^ agnsel
#print np.sum(emlineclass)

if save:
    if not he2_flag:    
        dfout = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                           'defagn':agnsel})
        dfout.to_csv('eco_emlineclass_bpt1.csv',index=False)
    
    else:
        dfout = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                          'defagn':agnsel,
                          'heiisel':agnsel4})
        dfout.to_csv('resolve_emlineclass_filtered_he2.csv',index=False)

#create alternate catalog name-based agn selector, print len
agn = (np.where(agnsel1 & resname)[0]) #for eco
#agn = (np.where(agnsel1 & econame)[0]) #for resolve
print ''
print 'AGN WITH ALTERNATE CATALOG NAME: ', len(agn)



#checking that plotted points are within the total data range
print ''
sfselpts = (len(np.where(sfsel)[0]))
compselpts = (len(np.where(compsel)[0]))
agnselpts = (len(np.where(agnsel)[0]))
#if he2_flag:
#    heiiselpts = (len(np.where(agnsel4)[0]))
#    totalselpts = sfselpts+seyfselpts+linerselpts+compselpts+agnselpts+ambigsel1pts+ambigsel2pts+heiiselpts
#else:
totalselpts = sfselpts+agnselpts+compselpts

print ("DATA POINTS: "),datalen
print ("TOTAL PLOTTED POINTS: "), totalselpts
print ("TOTAL PLOTTED POINTS OMITTED: "), datalen-totalselpts
print "* IF NUMBERS ABOVE ARE AT ALL NEGATIVE THEN THERE IS OVERPLOTTING"
print ("Definite Star Forming: "),sfselpts
print ("Composite: "),compselpts
print ("TOTAL KNOWN AGN: "),agnselpts

#display percent of each category
print ''
print "PERCENT OF GALAXIES IN EACH CATEGORY"
sfpercent = float(sfselpts)/float(datalen)*100
comppercent = float(compselpts)/float(datalen)*100
agnpercent = float(agnselpts)/float(datalen)*100

print ("Definite Star Forming:"), round(sfpercent, 2), ("%")
print("Composite: "),round(comppercent, 2), ("%")
print("TOTAL KNOWN AGN: "),round(agnpercent, 2), ("%")
print ("Percent Omitted: "), round((100-(sfpercent+comppercent+agnpercent)), 2), ("%")
print ''

###PLOTS###
#reference points in x-direction for demarcation lines on plots
refn2ha = np.linspace(-3.0, 0.35)

#NII/OIII plot
plt.figure(1)
ax = plt.subplot(111)
main1, = ax.plot(refn2ha, n2hamain(refn2ha), 'k', label = 'Main Line')
composite, = ax.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 'k-.', label = 'Composite Line')
#composite, = ax.plot(refn2ha, n2hacompmin(refn2ha), 'k--')
ax.set_xlim(-2,1)
ax.set_ylim(-1.25,1.5)
ax.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#data, = ax.plot(n2ha, o3hb, 'k.')
sfdata1, = ax.plot(n2ha[sfsel], o3hb[sfsel], 'k.', markersize = 3, alpha = 0.5, label = 'SF')
agndata1, = ax.plot(n2ha[agnsel], o3hb[agnsel],'rs', markersize = 8, label = 'Definite AGN')
compdata1, = ax.plot(n2ha[compsel], o3hb[compsel], 'bs', markersize = 8, label = 'Composite')
#comp2he2, = ax.plot(n2ha[agnsel4], o3hb[agnsel4], 'y*', label = 'HeII-selected AGN')
if he2_flag:
    agndata4, = ax.plot(n2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)
