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
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter.pkl'
df = pd.read_pickle(inputfile) #ECO catalog
#inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_all_dext.fits'
#inputfile = 'RESOLVE_SDSS_dext.fits'
#dat = Table.read(inputfile, format='fits')
#df = dat.to_pandas()
he2_flag = 1
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
name = df['NAME']
#resname = df['resname'] #for eco
#resname = resname != 'notinresolve'

#econame = df['econame']
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

emlineclass = sfsel ^ compsel ^ seyfsel ^ linersel ^ ambigsel1 ^ ambigsel2 ^ ambagnsel
#print np.sum(emlineclass)

if save:
    if not he2_flag:    
        dfout = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                          'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                          'sftoagn':ambigsel1, 'agntosf':ambigsel2})
        dfout.to_csv('resolve_emlineclass_filtered.csv',index=False)
    
    else:
        dfout = pd.DataFrame({'galname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                          'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                          'sftoagn':ambigsel1, 'agntosf':ambigsel2, 'heiisel':agnsel4})
        dfout.to_csv('resolve_emlineclass_filtered_he2.csv',index=False)

#create alternate catalog name-based agn selector, print len
#agn = (np.where(agnsel1 & resname)[0]) #for eco
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
#if he2_flag:
#    heiiselpts = (len(np.where(agnsel4)[0]))
#    totalselpts = sfselpts+seyfselpts+linerselpts+compselpts+agnselpts+ambigsel1pts+ambigsel2pts+heiiselpts
#else:
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
refoiha = np.linspace(-2.5, 0.0)
refsiiha = np.linspace(-2, 0.3)

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
seyfdata1, = ax.plot(n2ha[seyfsel], o3hb[seyfsel], 'co',  markersize = 5, label = 'Seyfert')
linerdata1, = ax.plot(n2ha[linersel], o3hb[linersel],'yo',  markersize = 5, label = 'LINER')
agndata1, = ax.plot(n2ha[ambagnsel], o3hb[ambagnsel],'rs', markersize = 8, label = 'Ambiguous AGN')
compdata1, = ax.plot(n2ha[compsel], o3hb[compsel], 'bs', markersize = 8, label = 'Composite')
ambig1data1, = ax.plot(n2ha[ambigsel1], o3hb[ambigsel1],'m^', markersize = 8, label = 'Star Forming --> AGN')
ambig2data1, = ax.plot(n2ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 8, label = 'AGN --> Star Forming')
#comp2he2, = ax.plot(n2ha[agnsel4], o3hb[agnsel4], 'y*', label = 'HeII-selected AGN')
if he2_flag:
    agndata4, = ax.plot(n2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)


#SII/OIII plot
plt.figure(2)
ax2 = plt.subplot(111)
main2, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Main Line')
liner, = ax2.plot(refsiiha, s2halinseyf(refsiiha), 'k--', label = 'Liner/Seyfert Division')
ax2.set_xlim(-2, 0.5)
ax2.set_ylim(-1.25,1.5)
ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
sfdata2, = ax2.plot(s2ha[sfsel], o3hb[sfsel], 'k.', markersize = 3, alpha = 0.5, label = 'SF')
seyfdata2, = ax2.plot(s2ha[seyfsel], o3hb[seyfsel], 'co', markersize = 5, label = 'Seyfert')
linerdata2, = ax2.plot(s2ha[linersel], o3hb[linersel],'yo', markersize = 5, label = 'LINER')
agndata2, = ax2.plot(s2ha[ambagnsel], o3hb[ambagnsel],'rs', markersize = 8, label = 'Ambiguous AGN')
compdata2, = ax2.plot(s2ha[compsel], o3hb[compsel], 'bs',markersize = 8, label = 'Composite')
ambig1data2, = ax2.plot(s2ha[ambigsel1], o3hb[ambigsel1], 'm^', markersize = 8, label = 'Star Forming --> AGN')
ambig2data2, = ax2.plot(s2ha[ambigsel2], o3hb[ambigsel2], 'g^', markersize = 8, label = 'AGN --> Star Forming ')
if he2_flag:
    agndata4, = ax2.plot(s2ha[agnsel4], o3hb[agnsel4],'ks', markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)


#OI/OIII plot
plt.figure(3)
ax3 = plt.subplot(111)
main3, = ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]), 'k', label = 'Main Line')
#main3, = ax3.plot(refoiha, o1hamain(refoiha), 'k', label = 'Main Line')
#liner, = ax3.plot(refoiha, o1halinseyf(refoiha), 'k--', label = 'Composite Line (Plot 1)')
liner2, = ax3.plot(refoiha, o1halinseyf(refoiha), 'k--', label = 'Liner/Seyfert Division')
ax3.set_xlim(-2.5, 0)
ax3.set_ylim(-1.25,1.5)
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
sfdata3, = ax3.plot(o1ha[sfsel], o3hb[sfsel], 'k.', alpha = 0.5, markersize = 3,label = 'SF')
seyfdata3, = ax3.plot(o1ha[seyfsel], o3hb[seyfsel], 'co', markersize = 5, label = 'Seyfert')
linerdata3, = ax3.plot(o1ha[linersel], o3hb[linersel],'yo', markersize = 5, label = 'LINER')
agndata3, = ax3.plot(o1ha[ambagnsel], o3hb[ambagnsel],'rs', markersize = 8, label = 'Ambiguous AGN')
compdata3, = ax3.plot(o1ha[compsel], o3hb[compsel], 'bs', markersize = 8, label = 'Composite')
ambig1data3, = ax3.plot(o1ha[ambigsel1], o3hb[ambigsel1],'m^', markersize = 8, label = 'MP-AGN')#'Star Forming --> AGN')
ambig2data3, = ax3.plot(o1ha[ambigsel2], o3hb[ambigsel2],'g^', markersize = 8, label = 'AGN --> Star Forming')
if he2_flag:
    agndata4, = ax3.plot(o1ha[agnsel4], o3hb[agnsel4],'ks',  markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)
plt.show()

##N2/HeII plot
if he2_flag:

    plt.figure()
    ax4 = plt.subplot(111)
    main4, = ax4.plot(refn2ha[refn2ha < -0.15], he2hbmain(refn2ha[refn2ha < -0.15]), 'k',  label = 'Main Line')
    liner4, = ax4.plot(refn2ha[refn2ha < 0.1], he2hblimit(refn2ha[refn2ha < 0.1]), 'k--', label = 'Composite Line (Plot 1)')
    ax4.set_xlim(-3., 1.)
    ax4.set_ylim(-3, 1.)
    ax4.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
    ax4.set_ylabel(r"$\rm \log([HeII]/H\beta)$", fontsize = 22)
    #data4, = ax4.plot(n2ha, he2hb, 'k.')
    sfdata4, = ax4.plot(n2ha[sfsel4], he2hb[sfsel4], 'k.',  markersize = 3, alpha = 0.5, label = 'SF')
    agndata4, = ax4.plot(n2ha[agnsel4], he2hb[agnsel4],'ks',  markersize = 8, mfc ='none', mew = 2, label = 'HeII-Selected AGN')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ', numpoints = 1)
    plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)


'''
ra=df.radeg
dec=df.dedeg
if 'fl_insample' in df.keys():
    flinsample = df.fl_insample
else:
    flinsample = np.ones(len(df), dtype = bool)
grpcz = df.grpcz
cz = df.cz
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars

inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (((flinsample | (np.log10(mbary) > 9.0)) & infall) | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
df = df[inobssample]
df = df[~np.isnan(df.h_alpha_flux_ext)]
df = df[~np.isnan(df.oiii_5007_flux_ext)]
df = df[~np.isnan(df.nii_6584_flux_ext)]
df = df[~np.isnan(df.nii_6548_flux_ext)]
df = df[~np.isnan(df.h_beta_flux_ext)]
df = df[~np.isnan(df.oi_6300_flux_ext)]
df = df[~np.isnan(df.sii_6717_flux_ext)]
df = df[~np.isnan(df.sii_6731_flux_ext)]
'''