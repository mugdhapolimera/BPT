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
#import pdb

#display catalog being used
print ''
print 'ECO RESULTS'
#print 'RESOLVE RESULTS'

#read in data
df = pd.read_pickle('ECO.pkl') #ECO catalog
#df = pd.read_pickle('RESOLVE.pkl') #RESOLVE catalog
df = df[~np.isnan(df.h_alpha_flux_ext)]
df = df[~np.isnan(df.oiii_5007_flux_ext)]
df = df[~np.isnan(df.nii_6584_flux_ext)]
df = df[~np.isnan(df.nii_6548_flux_ext)]
df = df[~np.isnan(df.h_beta_flux_ext)]
df = df[~np.isnan(df.oi_6300_flux_ext)]
df = df[~np.isnan(df.sii_6717_flux_ext)]
df = df[~np.isnan(df.sii_6731_flux_ext)]

#define alternate catalog names
name = df['name']
resname = df['resname'] #for eco
resname = resname != 'notinresolve'
#econame = df['econame'] #for resolve
#econame = econame != 'notineco'

#define (dwarf) blue E/S0's:
logmstar = df['logmstar']
#below certain mass, completeness 8.9
urcolor = df['modelu_rcorr']
morphel = df['morphel']
blueseq = (((urcolor < 1.457) & (logmstar <= 9.1)) | ((urcolor < (0.24*logmstar - 0.7)) & ((logmstar > 9.1) & (logmstar < 10.1))) | ((urcolor < 1.7) & (logmstar >=10.1)))
dwarf = (logmstar < 9) #dwarf galaxies
bluees0 = blueseq & (morphel == 'E') #Blue E/S0's

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

#create line ratios [NII]/H-alpha and [OIII]/H-beta
nii_sum = (df['nii_6584_flux_ext'] + df['nii_6548_flux_ext'])*3./4
# note that the ratio uses only the stronger line, but for S/N reasons we add
# the weaker and multiply by 3/4 since Chris Richardson says the canonical
# line ratio is 3:1 (this needs to be updated with a more precise number)
oiii = df['oiii_5007_flux_ext']
h_alpha = df['h_alpha_flux_ext']
h_beta = df['h_beta_flux_ext']
oi = df['oi_6300_flux_ext']
sii_sum = df['sii_6717_flux_ext'] + df['sii_6731_flux_ext']
#need to check Kewley paper and figure out if ratio used in nii_sum applies to sii_sum as well

#filter out bad data and merge gooddata with blue E/S0's
gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & (oi > 0) &
            (sii_sum > 0) & (h_beta > 0))

#data = gooddata & bluees0 #use all blue E/S0 data
#data = gooddata & bluees0 & dwarf #use all DWARF blue E/S0 data
#data = gooddata & dwarf #all DWARF galaxies
data = gooddata #use ALL galaxy data within catalog

#print data type being used
print ''
print 'ALL DATA'
#print 'DWARF ONLY DATA'
#print 'ALL BLUE E/S0 DATA'
#print 'DWARF ONLY BLUE E/S0 DATA'

#print total points shared with alternate catalog
sel = (np.where(data & resname)[0]) #for eco
#sel = (np.where(data & econame)[0]) #for resolve
print ''
print 'TOTAL DATA WITH ALTERNATE CATALOG NAME: ', len(sel)

nii_sum = nii_sum[data]
oiii = oiii[data]
oi = oi[data]
sii_sum = sii_sum[data]
h_beta = h_beta[data]
h_alpha = h_alpha[data]
subsetname = name[data]

#length of data to be used for debugging
datalen = np.sum(data)

# data ratios
n2ha = np.log10(nii_sum/h_alpha)
o3hb = np.log10(oiii/h_beta) # always the y-axis
o1ha = np.log10(oi/h_alpha)
s2ha = np.log10(sii_sum/h_alpha)

#Below are the selectors for the data to distinguish btwn: Seyferts, Composites,
#and AGN's based on the flux ratio diagnostic as understood via Kewley 2006.

#need to figure out why plot does not 'go back up' and consider adding to selectors
#to prohibit the data once the line 'goes back up'

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

#REFERENCE for cumulative plot selectors
seyfselr = seyfsel2 & seyfsel3
linerselr = linersel2 & linersel3

#cumulative plot selectors
sfsel = sfsel1 & sfsel2 & sfsel3 #definite star forming galaxies
compsel = compsel1  #composite galaxies
seyfsel = agnsel1 & seyfselr #Seyfert AGN galaxies
linersel = agnsel1 & linerselr #LINER AGN galaxies
ambigsel1 = sfsel1 & (agnsel2 | agnsel3) #SF in first plot, AGN in subsequent plot
ambigsel2 = agnsel1 & (sfsel2 | sfsel3) #AGN in first plot, SF in subsequent plot
ambagnsel = agnsel1 & ~seyfselr & ~linerselr & ~(sfsel2 | sfsel3) #Ambiguous AGN

emlineclass = sfsel ^ compsel ^ seyfsel ^ linersel ^ ambigsel1 ^ ambigsel2 ^ ambagnsel
#print np.sum(emlineclass)

dfout1 = pd.DataFrame({'name':name, 'emlineclass':emlineclass},index=name)
dfout1.fillna(False,inplace=True)
dfout2 = pd.DataFrame({'subsetname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                      'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                      'sftoagn':ambigsel1, 'agntosf':ambigsel2},index=subsetname)
dfout = pd.merge(dfout1,dfout2,how="left",left_index=True,right_index=True)
dfout.to_csv('ecoemlineclass.csv',index=False)

dfout3 = pd.DataFrame({'subsetname':subsetname, 'defstarform':sfsel, 'composite':compsel, 
                      'defseyf':seyfsel, 'defliner':linersel, 'ambigagn':ambagnsel,
                      'sftoagn':ambigsel1, 'agntosf':ambigsel2})
dfout3.to_csv('ecoemlineclass_subset.csv',index=False)

#create alternate catalog name-based agn selector, print len
agn = (np.where(agnsel1 & resname)[0]) #for eco
#agn = (np.where(agnsel1 & econame)[0]) #for resolve
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
fig = plt.figure(1)
plt.clf()

#reference points in x-direction for demarcation lines on plots
refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, 0.0)
refsiiha = np.linspace(-2, 0.3)

#NII/OIII plot
ax = fig.add_subplot(131)
main1, = ax.plot(refn2ha, n2hamain(refn2ha), 'k')
#composite, = ax.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]), 'k--')
composite, = ax.plot(refn2ha, n2hacompmin(refn2ha), 'k--')
ax.set_xlim(-2,1)
ax.set_ylim(-1.25,1.5)
ax.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
data, = ax.plot(n2ha, o3hb, 'k.')
sfdata1, = ax.plot(n2ha[sfsel], o3hb[sfsel], 'b.')
agndata1, = ax.plot(n2ha[ambagnsel], o3hb[ambagnsel],'yo')
seyfdata1, = ax.plot(n2ha[seyfsel], o3hb[seyfsel], 'rs')
linerdata1, = ax.plot(n2ha[linersel], o3hb[linersel],'gs')
compdata1, = ax.plot(n2ha[compsel], o3hb[compsel], 'co')
ambig1data1, = ax.plot(n2ha[ambigsel1], o3hb[ambigsel1],'m*', alpha = 0.5)
ambig2data1, = ax.plot(n2ha[ambigsel2], o3hb[ambigsel2],'k*', alpha = 0.5)

#SII/OIII plot
ax2 = fig.add_subplot(132)
main2, = ax2.plot(refsiiha, s2hamain(refsiiha), 'k')
liner, = ax2.plot(refsiiha, s2halinseyf(refsiiha), 'k--')
ax2.set_xlim(-2, 0.5)
ax2.set_ylim(-1.25,1.5)
ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
sfdata2, = ax2.plot(s2ha[sfsel], o3hb[sfsel], 'b.')
agndata2, = ax2.plot(s2ha[ambagnsel], o3hb[ambagnsel],'yo')
seyfdata2, = ax2.plot(s2ha[seyfsel], o3hb[seyfsel], 'rs')
linerdata2, = ax2.plot(s2ha[linersel], o3hb[linersel],'gs')
compdata2, = ax2.plot(s2ha[compsel], o3hb[compsel], 'co')
ambig1data2, = ax2.plot(s2ha[ambigsel1], o3hb[ambigsel1], 'm*', alpha = 0.5)
ambig2data2, = ax2.plot(s2ha[ambigsel2], o3hb[ambigsel2], 'k*', alpha = 0.5)


#OI/OIII plot
ax3 = fig.add_subplot(133)
#main, = ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]), 'k', label = 'Main Line')
main3, = ax3.plot(refoiha, o1hamain(refoiha), 'k', label = 'Main Line')
liner, = ax3.plot(refoiha, o1halinseyf(refoiha), 'k--', label = 'Composite Line (Plot 1)')
liner2, = ax3.plot(refoiha, o1halinseyf(refoiha), 'k--', label = 'Seyfert (Plots 2 & 3)')
ax3.set_xlim(-2.5, 0)
ax3.set_ylim(-1.25,1.5)
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
sfdata3, = ax3.plot(o1ha[sfsel], o3hb[sfsel], 'b.', label = 'Definite Star Forming')
agndata3, = ax3.plot(o1ha[ambagnsel], o3hb[ambagnsel],'yo', label = 'Ambiguous AGN')
seyfdata3, = ax3.plot(o1ha[seyfsel], o3hb[seyfsel], 'rs', label = 'Seyfert')
linerdata3, = ax3.plot(o1ha[linersel], o3hb[linersel],'gs', label = 'LINER')
compdata3, = ax3.plot(o1ha[compsel], o3hb[compsel], 'co', label = 'Composite')
ambig1data3, = ax3.plot(o1ha[ambigsel1], o3hb[ambigsel1],'m*', alpha = 0.5, label = 'Star Forming --> AGN')
ambig2data3, = ax3.plot(o1ha[ambigsel2], o3hb[ambigsel2],'k*', alpha = 0.5, label = 'AGN --> Star Forming')

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title = 'Galaxy Types: ')
plt.subplots_adjust(left=0.08, bottom=0.11, right=0.72, top=0.96, wspace=0.49, hspace=None)


plt.show()
