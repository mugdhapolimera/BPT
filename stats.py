# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 12:50:08 2019

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import os
from scipy import stats 

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
res = pd.read_pickle("RESOLVE_filter_new.pkl")
resflag = pd.read_csv("resolve_emlineclass_filter_new.csv")

eco = pd.read_pickle("ECO_filter_new.pkl")
ecoflag = pd.read_csv("eco_emlineclass_filter_new.csv")

full = pd.read_pickle("ECO+RESOLVE_filter_new.pkl")
fullflag = pd.read_csv("eco+resolve_emlineclass_filter.csv")

percent = pd.DataFrame({"Category": ["Definite SF", "SFing-AGN", "Composite",
                                     "Seyfert", "LINER", "Ambiguous AGN", 
                                     "AGN-to-SF"],
                        "RESOLVE":np.zeros(7), 
                        "ECO":np.zeros(7), 
                        "Overall": np.zeros(7)})
percent = percent.set_index('Category')

labels = {"defstarform": "Definite SF", "sftoagn": "SFing-AGN", 
          "composite": "Composite", "defseyf": "Seyfert", "defliner": "LINER", 
          "ambigagn": "Ambiguous AGN", "agntosf" : "AGN-to-SF"}

data = {"RESOLVE" : resflag, "ECO": ecoflag, "Overall": fullflag}

for label in labels.keys():
    for sample in data.keys():
        percent.loc[labels[label]][sample] = np.around(100.0*
                   np.sum(data[sample][label])/len(data[sample]),2)
        #percent.loc['Total'][sample] =        
print(percent[['RESOLVE', 'ECO', 'Overall']])#.to_latex())


def proportion_confint(count, nobs, alpha=0.05, method='normal'):
    pd_index = getattr(count, 'index', None)
    if pd_index is not None and callable(pd_index):
        # this rules out lists, lists have an index method
        pd_index = None
    count = np.asarray(count)
    nobs = np.asarray(nobs)

    q_ = count * 1. / nobs
    alpha_2 = 0.5 * alpha

    if method == 'normal':
        std_ = np.sqrt(q_ * (1 - q_) / nobs)
        dist = stats.norm.isf(alpha / 2.) * std_
        ci_low = q_ - dist
        ci_upp = q_ + dist
    return ci_low, ci_upp        

for x in percent.RESOLVE:
    print proportion_confint(x*len(res)/100.0, len(res))        

for x in percent.ECO:
    print proportion_confint(x*len(eco)/100.0, len(eco))        

for x in percent.Overall:
    print proportion_confint(x*len(full)/100.0, len(full))        
