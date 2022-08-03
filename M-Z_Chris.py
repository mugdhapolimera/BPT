# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 19:08:01 2021

@author: mugdhapolimera
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator, AutoMinorLocator
solarZ = 8.76
est = pd.read_csv("C:/Users/mugdhapolimera/github/nebulabayes/RESOLVE_jhu_snr5_coincident_open_M-7_0/RESOLVE_jhu_snr5_coincident_open_M-7_0_LOGZ.txt", 
                     sep = '\s+', names = ["name", "LOGZ", "err_up", "err_down"])
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5_dext_jhu.csv'
flags = "resolve_emlineclass_dext_snr5_jhu.csv"
res = pd.read_csv(inputfile)
emisline = pd.read_csv(flags)
emisline['name'] = emisline['galname']

est_res = pd.merge(est,res,on='name')
est_res = pd.merge(est_res,emisline,on='name')

res_sf = est_res[est_res['defstarform']==True]
res_agn = est_res[est_res['defstarform']==False]

mstar_sf = res_sf['logmstar']
mstar_agn = res_agn['logmstar']

logz_sf = res_sf['LOGZ']+solarZ
logz_agn = res_agn['LOGZ']+solarZ

m_stellar = np.linspace(8.5,11.5,1000)

tremonti = -1.492+1.847*m_stellar-0.08026*m_stellar**2.0
thomas_sf = 8.94 - np.log(1.0+(10**m_stellar/10**8.76)**-0.62)
pp04_N2 = 23.9049-5.62784*m_stellar+0.645142*m_stellar**2.0-0.0235065*m_stellar**3.0

fig = plt.figure(1)
plt.scatter(mstar_sf,logz_sf,s=0.5,color='blue',label='Def. SF')
plt.scatter(mstar_agn,logz_agn,s=0.5,color='red',label='AGN')
#plt.plot(m_stellar,tremonti,color='black',label='Tremonti')
#plt.plot(m_stellar,thomas_sf,color='black',label='Thomas SF',linestyle='--')
#plt.plot(m_stellar,pp04_N2,color='black',label='PP04 N2',linestyle=':')
plt.xlim(7.5,11.5)
plt.ylim(7.8,9.2)
plt.xlabel(r'log $M_*$')
plt.ylabel(r'12+log O/H')
plt.legend(frameon=False,loc='lower right')
#plt.show()
#plt.clf()