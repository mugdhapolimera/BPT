
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 11:13:54 2019

@author: mugdhapolimera

Creating target list for SOAR observations

"""
import pandas as pd
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.io.idl import readsav
import numpy as np
saveoutput = 0

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra/')
resolve = pd.read_csv('RESOLVE_full_raw.csv', index_col = 'name')
internal = readsav('resolvecatalog.dat')
internalphot = readsav('resolvecatalogphot.dat')
telarr = readsav('telarray.dat')['telarr']
code = readsav('telarray.dat')['code_2018A']

jhuflag = pd.read_csv('resolve_emlineclass_dext_snr5_jhu.csv', index_col = 'galname')
portflag = pd.read_csv('resolve_emlineclass_dext_snr5_port.csv', index_col = 'galname')
nsaflag = pd.read_csv('resolve_emlineclass_dext_snr5_nsa.csv', index_col = 'galname')

port = pd.read_csv('RESOLVE_full_snr5_dext_port.csv', index_col = 'name')[portflag.sftoagn]
jhu = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv', index_col = 'name')[jhuflag.sftoagn]
nsa = pd.read_csv('RESOLVE_full_snr5_dext_nsa.csv', index_col = 'name').loc[nsaflag.index.values[nsaflag.sftoagn]]
 
port = port[['radeg','dedeg']]   
jhu = jhu[['radeg','dedeg']]   
nsa = nsa[['radeg','dedeg']]

allunq = np.unique(list(jhuflag.index) + list(nsaflag.index) + list(portflag.index))   
#sftoagn = df[ambigsel1 & dwarf][['radeg','dedeg']]

#sfagn = pd.read_csv('uniquesfagn.csv')
unique = np.unique(list(jhu.index) + list(nsa.index) + list(port.index))

print('List of Unique SFing-AGN from JHU, Portsmouth and NSA Catalogs')
print(len(unique))
#print (jhu.loc[[x for x in jhu.index if x in unique]])
#print (port.loc[[x for x in port.index \
#                 if (x in unique) & (x not in jhu.index)]])
#print (nsa.loc[[x for x in nsa.index \
#                if (x in unique) & (x not in jhu.index) & (x not in port.index)]])

unq = jhu.loc[[x for x in jhu.index if x in unique]]
unq = unq.append(port.loc[[x for x in port.index \
                 if (x in unique) & (x not in jhu.index)]])
unq = unq.append(nsa.loc[[x for x in nsa.index \
                if (x in unique) & (x not in jhu.index) & (x not in port.index)]])

sfagn = unq
#nov_target = unq.index.values
#c = SkyCoord(ra=sfagn.radeg*u.degree, dec=sfagn.dedeg*u.degree, frame = 'fk5')
#sfagn['h'] = c.ra.hms.h
#sfagn['m'] = c.ra.hms.m
#sfagn['s'] = c.ra.hms.s
#
#sfagn['d'] = c.dec.dms.d
#sfagn['dm'] = c.dec.dms.m
#sfagn['ds'] = c.dec.dms.s
#
##print(sfagn)
#
##target range 
##November 2019 deadline- obs from jan to march 
#ra_low = 22
#ra_high = 3
#sel = (ra_low < sfagn.h) & (sfagn.h< ra_high)

#print(sfagn[sel])
#nugget_bluesetup = list(pd.read_csv('SOAR_2021A_all_bluesetup.csv').name)
#nugget_bluesetup.append('rs1036')
#nugget_bluesetup = ['rs0018','rs0049','rs0091','rs0169','rs0202','rs0239','rs0399','rs0608','rs0610','rs0708','rs0797','rs0801','rs0896','rs0978','rs1036','rs1205','rs1231']
nugget_bluesetup = ['rs0018','rs0049','rs0091','rs0136','rs0169','rs0202',
                    'rs0237','rs0239','rs0399','rs0608','rs0610','rs0693',
                    'rs0708','rs0741','rs0797','rs0801','rs0835','rs0896',
                    'rs0942','rs0978','rs1036','rs1205','rs1231','rs1415']

#nugget_broadsetup = list(pd.read_csv(r'C:\Users\mugdhapolimera\Desktop\nuggets_SOAR2020A_broad.txt',
#                      delim_whitespace=True).name)
#melanie = ['rs0661', 'rs1510', 'rs0106', 'rs0680']
#agn_priority = ['rs0124', 'rs0421', 'rs0909', 'rs1105', 'rs1047', 'rs1038']
#Below are sfagn without deep broad obs taken after 2019
agn_priority = ['rs0150','rs0421','rs0909','rs1038','rs1047',
                'rs1143','rs1195']
agn_priority = ['rf0338', 'rf0503']
#agn_priority = list(sfagn.index)
#agn_priority = ['rs0472','rs0930','rs1038',
#                'rs1105','rs1143','rs0167','rs0323','rs0430','rs0573','rs0933',
#                'rs0944','rs0981','rs1021','rs1055']

#lowsb - 'rs0150','rs0421', 'rs0775','rs1047',
#nuggets = ['rf0056','rf0283']
#semester = internal.infall
semester = internal.inspring
blueflag = ((internal.blue) | (internal.gemdone) | \
        (internal.koalablue) | (internal.saltlsblue))
redflag = ((internal.red) | (internal.gemdone) | \
        (internal.koalared) | (internal.saltls))
s0flag = internal.inobssample & (internal.morph <= 0) & semester & \
((((blueflag == 0) | (internal.broad == 0)) & (internal.obstag == 'A')) | \
(blueflag & (internal.broad == 0)))
HR = ((internal.red) | (internal.blue) | \
        (internal.gemdone) | (internal.koalared) | \
        (internal.koalablue) | (telarr == 'KINDONEHI')) | \
        ((internal.saltfp) | (internal.saltls) | \
                (internal.saltlsblue))
redsetup = (semester & \
           ((internalphot.smoothrestumag - internalphot.smoothrestrmag < 1.8)))
bluesetup = (semester & \
           ((internalphot.smoothrestumag - internalphot.smoothrestrmag > 1.8)))

broadsetup = list(internal.name[np.where(semester & (internal.broad == 0))])

redfillers = list(internal.name[np.where((HR == 0) & (internal.broad == 1) \
                                         & redsetup)])
bluefillers = list(internal.name[np.where((HR == 0) & (internal.broad == 1) \
                                          & bluesetup )])
s0 =  list(internal.name[np.where(s0flag)])
broadfillers = list(internal.name[np.where(semester & (internal.broad == 0) & 
                                           HR == 1)])
easybroadfillers = list(internal.name[np.where(semester & (internal.broad == 0) & 
                                           (HR == 1) & (internal.obstag == 'A'))])

#ra_filter = list(internal.name[np.where((internal.ra > 1.7*15) & (internal.ra < 5*15))])
#ra_filter = list(internal.name[np.where((internal.ra > 131.25) & (internal.ra < 236.25))])
soartarg = list(internal.name[np.where( internal.inobssample & 
                                       ((telarr == 'SALT') | \
                                        (telarr == 'SALTblue') | \
                                       (telarr == 'SOARABS') | \
                                       (telarr =='SOAREM') | \
                                       (telarr == 'KINDONE') | \
                                       (telarr == 'KINDONEHI')))])
#targetnames= list(np.intersect1d(redfillers+bluefillers+broadfillers,ra_filter))
easybroadandred = list(internal.name[np.where((internal.broad == 0) &
                                              (redflag== 0) & redsetup & \
                                              (internal.obstag == 'A'))])
#springlatera = list(internal.name[np.where((internal.ra > 204.8) & (internal.inspring) & \
#                                           (((redflag== 0) & redsetup) | \
#                                      ((blueflag ==0) & bluesetup) | \
#                                      internal.broad== 0))])
springlatera = list(internal.name[np.where((internal.ra > 14*15) & \
                                           (((redflag== 0) & redsetup) | \
                                      ((blueflag ==0) & bluesetup) | \
                                      internal.broad== 0))])
falllatera = list(internal.name[np.where(((internal.ra > 20*15)) & \
                                         (((redflag== 0) & redsetup) | \
                                      ((blueflag ==0) & bluesetup) | \
                                      internal.broad== 0))])
latera = list([springlatera + falllatera])[0]
#latera = list(np.intersect1d([springlatera + falllatera],list(internal.name[np.where(semester)])))

nuggets = list(np.intersect1d(nugget_bluesetup,
                              internal.name[np.where(blueflag == 0)]))
#targetnames= np.unique(nuggets+easybroadandred+latera)
#targetnames= np.unique(nuggets+easybroadandred+redfillers+easybroadfillers+falllatera)
agnnames= list(np.intersect1d(agn_priority,list(internal.name[np.where(semester)])))
agnnames= list(np.intersect1d(agn_priority,list(jhuflag.index[jhuflag.sftoagn])))
#targetnames= np.unique(targetnames+agnnames+nuggets)
maintargetnames= np.unique(agnnames+nuggets)
maintargetnames= list(np.intersect1d(maintargetnames,soartarg))

targetnames= np.unique(maintargetnames + easybroadandred + latera + easybroadfillers)
targetnames= list(np.intersect1d(targetnames,soartarg))+['rf0746', 'rf0747']

#targetnames= list(np.intersect1d(targetnames,broadsetup))
#skip red setup for these targets because no emission lines
#skip = ['rf0424', 'rf0475', 'rf0112', 'rf0304'] 
#targetnames = np.array([x for x in targetnames if x not in skip])

#targetnames= np.unique(nuggets+agnnames)
#targetnames = np.unique(list(s0+redfillers))
#list(sfagn[sel].index.values)
#resolve = readsav('../../../SDSS_spectra/resolvecatalog.dat')
#mu_r = pd.DataFrame(dict(zip(resolve.name, resolve.ifusb)))
targets = resolve.loc[targetnames][['radeg','dedeg']]
#print targets
c = SkyCoord(ra=targets.radeg*u.degree, dec=targets.dedeg*u.degree, frame = 'fk5')
targets['h'] = c.ra.hms.h
targets['m'] = c.ra.hms.m
targets['s'] = c.ra.hms.s

targets['d'] = c.dec.dms.d
targets['dm'] = c.dec.dms.m
targets['ds'] = c.dec.dms.s
#
targets['sample'] = 'resolve'
targets = targets.sort_values(by='radeg')
#targets_unsorted = targets.copy()
#sortra1 = targets.radeg > 330.0
#sortra2 = targets.radeg < 45.0
#
#targets = targets_unsorted[sortra1]
#targets = targets.append(targets_unsorted[sortra2])
#
#print targets
ra = []
dec = []
rsys = []

for x in targets.index.values:
    ra.append(str(int(targets.loc[x].h)).zfill(2)+':'+\
              str(int(targets.loc[x].m)).zfill(2)+':'+\
              '{0:05.2f}'.format(np.abs(targets.loc[x].s)))
    if targets.loc[x].dm < 0:
        dec.append('-'+str(int(np.abs(targets.loc[x].d)))+':'+\
              str(int(np.abs(targets.loc[x].dm))).zfill(2)+':'+\
              '{0:05.2f}'.format(np.abs(targets.loc[x].ds)))
    else:
        dec.append(str(int(targets.loc[x].d)).zfill(2)+':'+\
              str(int(np.abs(targets.loc[x].dm))).zfill(2)+':'+\
              '{0:05.2f}'.format(np.abs(targets.loc[x].ds)))
    ndx = np.where(internal.name == x)        
#    if x in melanie: 
#        targets['sample'].loc[x] = 'nascent - red (usual 666 counts S/N ~ 25) '
#    if x in nuggets: 
#        targets['sample'].loc[x] = 'nuggets - ~60 min blue setup (25 S/N for continuum) + std. broad'
#    elif (x in s0) & (internal.broad[ndx] == 0): 
#        targets['sample'].loc[x] = 'S0 - broad only'
#    elif (x in s0) & (blueflag[ndx] == 0) & (internal.broad[ndx] == 0): 
#        targets['sample'].loc[x] = 'S0 - blue and broad'
#    elif (x in s0) & (blueflag[ndx] == 0) & (internal.broad[ndx] == 1): 
#        targets['sample'].loc[x] = 'S0 - blue only'
    if x in agn_priority: 
        targets['sample'].loc[x] = 'SFing-agn - broad 1-2 hrs; S/N ~ 15, 225 counts for [OI]; deep red ~1hr'      
#    elif x in list(targetlist.index.values) and x in agn_priority: 
#        targets['sample'].loc[x] = 'SFing-agn - broad 1-2 hrs [OI] S/N ~ 15; red S/N ~ 25?'      
#    elif x in nugget_bluesetup: 
#        targets['sample'].loc[x] = 'nuggets - usual blue setup'
    elif x in easybroadandred: 
        targets['sample'].loc[x] = 'easy std. broad and red setup'
    elif x in easybroadfillers: 
        targets['sample'].loc[x] = 'easy std. broad setup filler'
    elif (x in latera) or (x in springlatera): 
        targets['sample'].loc[x] = 'late RA fillers; standard obs. '
    elif x in redfillers: 
        targets['sample'].loc[x] = 'std. red setup; filler obs.'
#    elif x in broadfillers: 
#        targets['sample'].loc[x] = 'std. broad'
    elif (x not in nugget_bluesetup): #x in nugget_broadsetup: 
        targets['sample'].loc[x] = 'nuggets - usual broad setup'
    elif x in nugget_bluesetup: 
        targets['sample'].loc[x] = 'nuggets - usual broad setup; 25 S/N for blue setup continuum'
    
#    rsys.append('Vega')
    #ndx = np.where(resdata.name == x)
    #mu_r = resdata.ifusb[ndx]
#print 'name           ra        dec       ra2           dec2     PAout  PAin    Re   min  iout   iin mur90 M_r  LR HR obstag sample'
sdss = np.genfromtxt('specdr82spring.txt', dtype = None,
                     names = ['namedr8','radr8','decdr8',
                                                    'specdr8flagstring'])

if saveoutput:
    with open(r'C:\Users\mugdhapolimera\Desktop\targetlist_forus_052422_test.txt', 'w') as the_file:
        the_file.write('name           ra        dec       ra2           dec2        PAout   PAin    Re     min    iout   iin    mur90   M_r   LR  HR  obstag sample\n')
        for x in range(len(targets.index.values)):
        #for x in range(1):
        #    galname = 'rs0087'
            galname = targets.index.values[x]
            ndx = np.where(internal.name == galname)[0][0]
            Re = internalphot['radr50p'][ndx]
            b_a = internalphot['b_a'][ndx]
            b_a_alt = internalphot['b_ainnerdisk'][ndx]
            minor_axis = 2.*1.3*Re*b_a
            sini=min(np.sqrt((1-b_a**2)/(1-0.2**2)), 1.0)
            sini_alt=min(np.sqrt((1-b_a_alt**2)/(1-0.2**2)),1.0)
            inclination=np.arcsin(sini)*180./np.pi
            if np.isfinite(inclination):
                if inclination == 0:
                    inclination = 90.
            inclination_alt = np.arcsin(sini_alt)*180./np.pi
            if (np.isfinite(inclination_alt)):
                if inclination_alt == 0:
                    inclination_alt = 90.
            ind2=np.where(sdss['namedr8'] == galname)
            if sdss['specdr8flagstring'][ind2] == -9999:
                specdr8flag=0
            elif sdss['specdr8flagstring'][ind2] == 11111:
                specdr8flag=2
            else:    
                specdr8flag=1
                
            if internal.broad[ndx] : 
                LR = 1
            else: 
                LR = 0
            if ((internal.red[ndx]) or (internal.blue[ndx]) or \
            (internal.gemdone[ndx]) or (internal.koalared[ndx]) or \
            (internal.koalablue[ndx]) or (telarr[ndx] == 'KINDONEHI')) or \
            ((internal.saltfp[ndx]) or (internal.saltls[ndx]) or \
                    (internal.saltlsblue[ndx])):
                HR = 1
            else: 
                HR = 0
        #    if targets['sample'].loc[x] == 'sfing-agn':
        #        LR = 1
        #    HR = 0

            string =  str(x+1).zfill(3)+'_'+galname+\
            '   '+'{0:08.4f}'.format(targets.loc[galname]['radeg'])+\
            '   '+'{0:06.4f}'.format(targets.loc[galname]['dedeg'])+\
            '   '+ra[x]+'   '+dec[x]+\
            '   '+'{0:05.1f}'.format(internalphot['pa'][ndx])+\
            '   '+'{0:05.1f}'.format(internalphot['painnerdisk'][ndx])+\
            '   '+'{0:04.1f}'.format(Re)+\
            '   '+'{0:04.1f}'.format(minor_axis)+\
            '   '+'{0:04.1f}'.format(inclination)+\
            '   '+'{0:04.1f}'.format(inclination_alt)+\
            '   '+'{0:.1f}'.format(internalphot.mur90[ndx])+\
            '   '+'{0:.1f}'.format(resolve.loc[galname]['absrmag'])+\
            '   '+'{}'.format(LR)+\
            '   '+'{}'.format(HR)+\
            '   '+internal['obstag'][ndx]+str(specdr8flag)+\
            '   '+targets.loc[galname]['sample']+'\n'
        
            the_file.write(string)
    #    print str(x+1).zfill(3)+'_'+galname+\
    #    '   '+'{0:08.4f}'.format(targets.loc[galname]['radeg'])+\
    #    '   '+'{0:06.4f}'.format(targets.loc[galname]['dedeg'])+\
    #    '   '+ra[x]+'   '+dec[x]+\
    #    '   '+'{0:05.1f}'.format(internalphot['pa'][ndx])+\
    #    '   '+'{0:05.1f}'.format(internalphot['painnerdisk'][ndx])+\
    #    '   '+'{0:04.1f}'.format(Re)+\
    #    '   '+'{0:04.1f}'.format(minor_axis)+\
    #    '   '+'{0:04.1f}'.format(inclination)+\
    #    '   '+'{0:04.1f}'.format(inclination_alt)+\
    #    '   '+'{0:.1f}'.format(internalphot.mur90[ndx])+\
    #    '   '+'{0:.1f}'.format(resolve.loc[galname]['absrmag'])+\
    #    '   '+'{}'.format(LR)+\
    #    '   '+'{}'.format(HR)+\
    #    '   '+internal['obstag'][ndx]+str(specdr8flag)+\
    #    '   '+targets.loc[galname]['sample']

for x in range(len(targets.index.values)):
    galname = targets.index.values[x]
    ndx = np.where(internal.name == galname)[0][0]
    print str(x+1).zfill(3)+'_'+galname+\
    '   '+ra[x]+'   '+dec[x]+'   2000   '#+telarr[ndx]
#M_r = list(df.loc[nov_target]['absrmag'])                
#data = zip(nov_target,ra,dec, M_r,rsys)
#nov = pd.DataFrame(data = data, columns = ['Name','RAJ200','DecJ200','r','r_sys'])
#nov.to_csv('November_targetlist.csv')