
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 12:50:08 2019

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import os
from scipy import stats 
from astropy.stats import binom_conf_interval

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')

resolve = pd.read_csv('RESOLVE_full_raw.csv', index_col = 'name')
resjhuflag = pd.read_csv('resolve_emlineclass_dext_snr5_jhu.csv', index_col = 'galname')
resportflag = pd.read_csv('resolve_emlineclass_dext_snr5_port.csv', index_col = 'galname')
resnsaflag = pd.read_csv('resolve_emlineclass_dext_snr5_nsa.csv', index_col = 'galname')

resport = pd.read_csv('RESOLVE_full_snr5_dext_port.csv', index_col = 'name')#[portflag.sftoagn]
resjhu = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv', index_col = 'name')#[jhuflag.sftoagn]
resnsa = pd.read_csv('RESOLVE_full_snr5_dext_nsa.csv', index_col = 'name')#.loc[nsaflag.index.values[nsaflag.sftoagn]]

fulljhu = pd.read_csv("ECO+RESOLVE_snr5_dext_jhu.csv")
fullport = pd.read_csv("ECO+RESOLVE_snr5_dext_port.csv")
fullnsa = pd.read_csv("ECO+RESOLVE_snr5_dext_nsa.csv")

fullportflag = pd.read_csv("eco+resolve_emlineclass_dext_snr5_port.csv")
fullnsaflag = pd.read_csv("eco+resolve_emlineclass_dext_snr5_nsa.csv")

eco = pd.read_csv('ECO_full_raw.csv', index_col = 'name')
os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra\ECO\SEL')
ecojhuflag = pd.read_csv('eco_emlineclass_dext_snr5_jhu.csv', index_col = 'galname')
ecoportflag = pd.read_csv('eco_emlineclass_dext_snr5_port.csv', index_col = 'galname')
econsaflag = pd.read_csv('eco_emlineclass_dext_snr5_nsa.csv', index_col = 'galname')

ecoport = pd.read_csv('ECO_full_snr5_dext_port.csv', index_col = 'name')#[portflag.sftoagn]
ecojhu = pd.read_csv('ECO_full_snr5_dext_jhu.csv', index_col = 'name')#[jhuflag.sftoagn]
econsa = pd.read_csv('ECO_full_snr5_dext_nsa.csv', index_col = 'name')#.loc[nsaflag.index.values[nsaflag.sftoagn]]

fulljhuflag = pd.read_csv("eco+resolve_emlineclass_dext_snr5_jhu.csv")

percent = pd.DataFrame({"Category": ["Definite SF", "SFing-AGN", "Composite",
                                     "Seyfert", "LINER", "Ambiguous AGN", 
                                     "AGN-to-SF"],
                        "JHU RESOLVE":np.zeros(7),
                        "JHU ECO":np.zeros(7), 
                        "JHU Overall": np.zeros(7),
                        "Port RESOLVE":np.zeros(7),
                        "port ECO":np.zeros(7), 
                        "port Overall": np.zeros(7),
                        "nsa RESOLVE":np.zeros(7),
                        "nsa ECO":np.zeros(7), 
                        "nsa Overall": np.zeros(7)})
percent = percent.set_index('Category')

labels = {"defstarform": "Definite SF", "sftoagn": "SFing-AGN", 
          "composite": "Composite", "defseyf": "Seyfert", "defliner": "LINER", 
          "ambigagn": "Ambiguous AGN", "agntosf" : "AGN-to-SF"}

data = {"JHU RESOLVE" : resjhuflag, "JHU ECO": ecojhuflag, "JHU Overall": fulljhuflag,
        "Port RESOLVE" : resportflag, "port ECO": ecoportflag, "port Overall": fullportflag,
        "nsa RESOLVE" : resnsaflag, "nsa ECO": econsaflag, "nsa Overall": fullnsaflag,}



for label in labels.keys():
    for sample in data.keys():
        percent.loc[labels[label]][sample] = np.around(100.0*
                   np.sum(data[sample][label])/len(data[sample]),2)
        #percent.loc['Total'][sample] =        
print(percent[["JHU RESOLVE","JHU ECO","JHU Overall",
               "Port RESOLVE","port ECO","port Overall",
               "nsa RESOLVE","nsa ECO","nsa Overall"]].to_latex())
#print percent

#def proportion_confint(count, nobs, alpha=0.05, method='normal'):
#    pd_index = getattr(count, 'index', None)
#    if pd_index is not None and callable(pd_index):
#        # this rules out lists, lists have an index method
#        pd_index = None
#    count = np.asarray(count)
#    nobs = np.asarray(nobs)
#
#    q_ = count * 1. / nobs
#    alpha_2 = 0.5 * alpha
#
#    if method == 'normal':
#        std_ = np.sqrt(q_ * (1 - q_) / nobs)
#        dist = stats.norm.isf(alpha / 2.) * std_
#        ci_low = q_ - dist
#        ci_upp = q_ + dist
#    return ci_low, ci_upp        
#
#for x in percent.RESOLVE:
#    print proportion_confint(x*len(res)/100.0, len(res))        
#
#for x in percent.ECO:
#    print proportion_confint(x*len(eco)/100.0, len(eco))        
#
#for x in percent.Overall:
#    print proportion_confint(x*len(full)/100.0, len(full))        

resjhudwarf = resjhu[resjhu.logmstar < 9.5]
ecojhudwarf = ecojhu[ecojhu.logmstar < 9.5]
resjhuagn = resjhuflag.defagn| resjhuflag.composite | resjhuflag.sftoagn | \
            resjhuflag.agntosf
ecojhuagn = ecojhuflag.defagn| ecojhuflag.composite | ecojhuflag.sftoagn | \
            ecojhuflag.agntosf
ecojhudwarfagn = ecojhudwarf[ecojhuagn]
resjhudwarfagn = resjhudwarf[resjhuagn]

resportdwarf = resport[resport.logmstar < 9.5]
ecoportdwarf = ecoport[ecoport.logmstar < 9.5]
resportagn = resportflag.defagn| resportflag.composite | resportflag.sftoagn | \
            resportflag.agntosf
ecoportagn = ecoportflag.defagn| ecoportflag.composite | ecoportflag.sftoagn | \
            ecoportflag.agntosf
ecoportdwarfagn = ecoportdwarf[ecoportagn]
resportdwarfagn = resportdwarf[resportagn]

resnsadwarf = resnsa[resnsa.logmstar < 9.5]
econsadwarf = econsa[econsa.logmstar < 9.5]
resnsaagn = resnsaflag.defagn| resnsaflag.composite | resnsaflag.sftoagn | \
            resnsaflag.agntosf
econsaagn = econsaflag.defagn| econsaflag.composite | econsaflag.sftoagn | \
            econsaflag.agntosf
econsadwarfagn = econsadwarf[econsaagn]
resnsadwarfagn = resnsadwarf[resnsaagn]

indexes = ['JHU', 'Portsmouth', 'NSA', 'JHU OR Portsmouth OR NSA', 
           'JHU AND Portsmouth AND NSA']
colname = {'JHU': 'jhu', 'Portsmouth':'port', 'NSA':'nsa',
           'JHU OR Portsmouth OR NSA': 'jhu+port+nsa',
           'JHU AND Portsmouth AND NSA': 'jhu&port&nsa'}
print('\\begin{table*} \ \n \centering \
      \n \t \\begin{tabular}{C{3cm} | C{2cm} C{2cm} C{2cm} |C{2cm} C{2cm} C{2cm} } \
      \n \t \\toprule \ \n \t {} & \multicolumn{3}{c|}{RESOLVE} & \multicolumn{3}{c}{ECO} \\\\ \hline \
         \n \t Category &    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN &  \
    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN \\\\ \n \t \hline')
         
for i in range(len(indexes)):
    print ('\t          &  &   &  &  &   & \\\\')
    index = indexes[i] 
    
    if ('OR' in index): 
        rdwarf = len(list(set(resjhudwarf.index.values ) |set(resportdwarf.index.values ) |\
                          set(resnsadwarf.index.values)))
        rdwarfagn = len(list(set(resjhudwarfagn.index.values ) |set(resportdwarfagn.index.values ) \
                             |set(resnsadwarfagn.index.values)))
        edwarf = len(list(set(ecojhudwarf.index.values ) |set(ecoportdwarf.index.values ) |\
                          set(econsadwarf.index.values)))
        edwarfagn = len(list(set(ecojhudwarfagn.index.values ) |set(ecoportdwarfagn.index.values ) \
                             |set(econsadwarfagn.index.values)))

    elif('AND' in index):
        rdwarf = len(list(set(resjhudwarf.index.values ) & set(resportdwarf.index.values ) \
                          & set(resnsadwarf.index.values)))
        rdwarfagn = len(list(set(resjhudwarfagn.index.values ) & set(resportdwarfagn.index.values ) \
                             & set(resnsadwarfagn.index.values)))
        edwarf = len(list(set(ecojhudwarf.index.values ) & set(ecoportdwarf.index.values ) & \
                          set(econsadwarf.index.values)))
        edwarfagn = len(list(set(ecojhudwarfagn.index.values ) & set(ecoportdwarfagn.index.values ) \
                             & set(econsadwarfagn.index.values)))

    elif (index == 'JHU'):
        rdwarf = len(resjhudwarf)
        rdwarfagn = len(resjhudwarfagn)        
        edwarf = len(ecojhudwarf)
        edwarfagn = len(ecojhudwarfagn)        

    elif (index == 'Portsmouth'):
        rdwarf = len(resportdwarf)
        rdwarfagn = len(resportdwarfagn)        
        edwarf = len(ecoportdwarf)
        edwarfagn = len(ecoportdwarfagn)        

    elif (index == 'NSA'):
        rdwarf = len(resnsadwarf)
        rdwarfagn = len(resnsadwarfagn)        
        edwarf = len(econsadwarf)
        edwarfagn = len(econsadwarfagn)        

    rdwarfagnpc = round((100.0*rdwarfagn/rdwarf),2)
    r_edown, r_eup = 100.0*binom_conf_interval(rdwarfagn, rdwarf) - rdwarfagnpc
    r_edown = round(-r_edown,2)
    r_eup = round(r_eup,2)
    edwarfagnpc = round((100.0*edwarfagn/edwarf),2)
    e_edown, e_eup = 100.0*binom_conf_interval(edwarfagn, edwarf) - edwarfagnpc
    e_edown = round(-e_edown,2)
    e_eup = round(e_eup,2)

    print '\t'+index+' & '+ str(rdwarf)+' & ',str(rdwarfagn)+' & $'+str(rdwarfagnpc)+\
          '^{+'+str(r_eup)+'}'+'_{'+str(-r_edown)+'}$'+' & '\
          + str(edwarf)+' & ',str(edwarfagn)+' & $'+str(edwarfagnpc)+\
          '^{+'+str(e_eup)+'}'+'_{'+str(-e_edown)+'}$\\\\'
    print '\t \hline \n &  &   &  &  &   & \\\\'

print('\t \hline \n \t \end{tabular} \n \label{table:2} \n \end{table*}')

#print('\\begin{table*} \ \n \centering \
#      \n \t \\begin{tabular}{| C{2cm} | C{2cm} |C{2cm} |C{2cm} |C{2cm} |C{2cm} |C{2cm} |} \
#      \n \t \\toprule \ \n \t {} & \multicolumn{3}{|c|}{RESOLVE} & \multicolumn{3}{|c|}{ECO} \\\\ \hline \
#         \n \t Category &    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN &  \
#    Dwarf SEL Sample &  Dwarf SEL AGN &  \% of SEL AGN \\\\ \n \t \hline')
#
#indexes = ['Confidence $>=$ 0', 'Confidence $>=$ 1', 'Confidence $=$ 2']
#for i in [0,1,2]:
#    index = indexes[i]
#    if i == 2:
#        rdwarf = np.sum((resdwarfconf.jhu != 0) & (resdwarfconf.port != 0))
#        rdwarfagn = np.sum(resdwarfconf['confidence_level'] == 2)        
#        edwarf = np.sum((ecodwarfconf.jhu != 0) & (ecodwarfconf.port != 0))
#        edwarfagn = np.sum(ecodwarfconf['confidence_level'] == 2)        
#    else:        
#        rdwarf = len(resdwarfconf)
#        rdwarfagn = np.sum(resdwarfconf['confidence_level'] >= i)        
#        edwarf = len(ecodwarfconf)
#        edwarfagn = np.sum(ecodwarfconf['confidence_level'] >= i)        
#
#    rdwarfagnpc = round((100.0*rdwarfagn/rdwarf),2)
#    r_edown, r_eup = 100.0*binom_conf_interval(rdwarfagn, rdwarf) - rdwarfagnpc
#    r_edown = round(-r_edown,2)
#    r_eup = round(r_eup,2)
#
#    edwarfagnpc = round((100.0*edwarfagn/edwarf),2)
#    e_edown, e_eup = 100.0*binom_conf_interval(edwarfagn, edwarf) - edwarfagnpc
#    e_edown = round(-e_edown,2)
#    e_eup = round(e_eup,2)
#
#    print '\t'+index+' & '+ str(rdwarf)+' & ',str(rdwarfagn)+' & $'+str(rdwarfagnpc)+\
#          '^{+'+str(r_eup)+'}'+'_{'+str(-r_edown)+'}$'+' & '\
#          + str(edwarf)+' & ',str(edwarfagn)+' & $'+str(edwarfagnpc)+\
#          '^{+'+str(e_eup)+'}'+'_{'+str(-e_edown)+'}$\\\\'
#print('\t \hline \n \t \end{tabular} \n \label{table:2} \n \end{table*}')
