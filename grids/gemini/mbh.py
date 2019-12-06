# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 09:44:10 2019
@author: mugdhapolimera

Calculate the err(FWHM_Ha_Broad) needed to get estimates of M_BH within the 
error of the method
"""
import numpy as np

def prod_err(x,y, xerr, yerr):
    return x*y*np.sqrt((xerr/x)**2 + (yerr/y)**2)
    
def log_err(x,xerr):
    return xerr/(np.log(10)*x)

#gal = {'z' : 0.0325, 'flux_ha' : 582, 'flux_ha_err': 28, 
#       'FWHM': 598, 'FWHM_err':69, 'name': 'F'}
#gal = {'z' : 0.0327, 'flux_ha' : 118, 'flux_ha_err': 16, 
#       'FWHM': 636, 'FWHM_err':69, 'name': '11'}
#gal = {'z' : 0.02145, 'Lha' : 10**39.91, 'Lha_err': 0.018*np.log(10)*10**39.91, 
#       'FWHM': 1449.69, 'FWHM_err':55.50, 'name': 'rs1036'}
gal = {'z' : 0.0194, 'flux_ha' : 145.77, 'flux_ha_err': 10*2.64, 
       'FWHM': 500, 'FWHM_err':33*2.354, 'name': 'rs0010_theory'}

f = 3
f_err = 0.55
a = 0.45
a_err = 0.03
b = 2.06
b_err = 0.06
d = (3e5*gal['z']/70) * 3.086e+24
FWHM_ha = gal['FWHM'] / 1e3
#Lha = 10**40.09
if 'Lha' in gal.keys():
    Lha = gal['Lha'] /1e42
    Lha_err = gal['Lha_err'] /1e42
else:
    Lha = gal['flux_ha'] * 1e-17 *(4*3.14*d**2) / 1e42
    Lha_err = gal['flux_ha_err'] * 1e-17 *(4*3.14*d**2) / 1e42
FWHM_ha_err = gal['FWHM_err'] / 1e3


log_mbh = np.log10(f) + 6 + a*np.log10(Lha) + b*np.log10(FWHM_ha)

var_log_mbh = log_err(f,f_err)**2 + prod_err(a,np.log10(Lha/1e42),a_err,\
              log_err(Lha/1e42,Lha_err/1e42))**2 + \
              prod_err(b,np.log10(FWHM_ha/1000.0),b_err,\
                       log_err(FWHM_ha/1000.0,FWHM_ha_err))**2
var_log_mbh = log_err(f,f_err)**2 + prod_err(a,np.log10(Lha),a_err,\
              log_err(Lha,Lha_err))**2 + \
              prod_err(b,np.log10(FWHM_ha),b_err,\
                       log_err(FWHM_ha,FWHM_ha_err))**2

err_log_mbh = np.sqrt(var_log_mbh)
err_mbh = err_log_mbh*np.log(10)*(10**log_mbh)
print(log_mbh,err_log_mbh, err_mbh, (10**log_mbh)/err_mbh)

rblr = (10**(1.216 + 0.47*np.log10(Lha))) *8.4e-7 #kpc
print(rblr)

#Radius of sphere of influence
sig = 30*10**3 #m/s - assume from M*-sig_bulge relation
G = 6.67e-11
M_sun = 2e30
M = 10**6#log_mbh
r_soi = (G*M*M_sun/sig**2) * 3.2e-17 #pc
print(r_soi)