import numpy as np
from astropy.io import fits
import pandas as pd

#source code from Ashley Bittner and Margie Bruff
#last updated 03/27/2018 by Carlynn Ferguson

#ECO data merge
#'''
#ECOdata = np.genfromtxt('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_blend_dext_new.csv', delimiter=",",dtype=None,names=True)
#ECOdata = pd.DataFrame(data=ECOdata,index=ECOdata['name'])

ECOdata = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_inobssample.csv')
ECOdata.index = ECOdata.name
#ECOfits = fits.open('ECO_SDSS_dext.fits')
#data = np.array(ECOfits[1].data[np.argsort(ECOfits[1].data['NAME'])])
#data = data.byteswap().newbyteorder()
#ECOfits = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_SDSS_cont_flux.csv")#pd.DataFrame(data=data,index=data['NAME'])
ECOfits = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/mid_IR/ECO_WISE_good.csv")#pd.DataFrame(data=data,index=data['NAME'])
ECOfits.index = ECOfits.name

ECOdf = pd.merge(ECOdata,ECOfits,how="left",left_index=True,right_index=True)
ECOdf.to_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_WISE_good.csv')
#'''

#RESOLVE data merge

#RESOLVEdata = np.genfromtxt('C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_blend_dext_new.csv', delimiter=",",dtype=None,names=True)
#RESOLVEdata = pd.DataFrame(data=RESOLVEdata,index=RESOLVEdata['name'])
RESOLVEdata = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_inobssample.csv')
RESOLVEdata.index = RESOLVEdata.name

#RESOLVEfits = fits.open('RESOLVE_SDSS_dext.fits')
#data = np.array(RESOLVEfits[1].data[np.argsort(RESOLVEfits[1].data['NAME'])])
#data = data.byteswap().newbyteorder()
#RESOLVEfits = pd.DataFrame(data=data,index=data['NAME'])
#RESOLVEfits= pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_cont_flux.csv')
RESOLVEfits = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/mid_IR/RESOLVE_WISE_good.csv")
RESOLVEfits.index = RESOLVEfits.name
RESOLVEdf = pd.merge(RESOLVEdata,RESOLVEfits,how="left",left_index=True,right_index=True)
RESOLVEdf.to_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_WISE_good.csv')
