from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_blend_dext_new.csv'

res = pd.read_csv(inputfile)
catalog = readsav('../SDSS_spectra/resolvecatalog.dat')


cz = catalog['vlg'] #col[['cz']]
z = cz/3e5
radeg = catalog['ra']
rahr = radeg/15
rarad = radeg * (np.pi / 180.)

sfagnnames = res.index.values
sel = np.where(catalog['name'] == sfagnnames)
cz_targ = catalog['vlg'][sel]
z_targ = cz_targ/3e5
radeg_targ = catalog['ra'][sel]
rahr_targ = radeg_targ / 15
rarad_targ = radeg_targ * (np.pi / 180.)


thisr = 0.019222
thistheta = 0.424053

plt.figure(1)
plt.clf()
ax = plt.subplot(111, projection = 'polar')
ax.plot(rarad, z, 'b.', markersize = 3)
ax.plot(rarad_targ, z_targ, 'r.', markersize = 12)

#ax.annotate('rf0250', xy = (thisr, thistheta), xytext = (0.05, 0.05), textcoords = 'figure fraction', arrowprops=dict(facecolor='black', shrink=0.05), horizontalalignment='left', verticalalignment='bottom')

#ax.annotate('rf0250', xy=(thisr, thistheta),  xycoords='axes fraction',
                #horizontalalignment='center', verticalalignment='center')


xT=plt.xticks()[0]
xL=['0h',r'3h',r'6h',r'9h',\
    r'12h',r'15h',r'18h',r'21h']
plt.xticks(xT, xL, fontsize = 20)
plt.yticks(fontsize = 20)

#vlgmin=min(resvlg[where(inobssample)])
#vlgmax=max(resvlg[where(inobssample)]
#vlgboundf=[vlgmin,vlgmax,fltarr(1000)+vlgmax,vlgmin,fltarr(1000)+vlgmin]
#raradboundf=[330.,330.,findgen(1000)*(405.-330.)/999.+330.,405.,findgen(1000)*(330.-405.)/999.+405.]*!pi/180.
#vlgbounds=[vlgmin,vlgmax,fltarr(1000)+vlgmax,vlgmin,fltarr(1000)+vlgmin]
#raradbounds=[131.25,131.25,findgen(1000)*(236.25-131.25)/999.+131.25,236.25,findgen(1000)*(131.25-236.25)/999.+236.25]*!pi/180.
#!psym=0
#oplot,vlgboundf/70.,raradboundf,/polar
#oplot,vlgbounds/70.,raradbounds,/polar
#plt.title('RESOLVE Fan Plot', y = 1.08)

#plt.show()

fullres = pd.read_csv('RESOLVE_liveMay2020.csv')
fullresndx = [np.where(catalog.name == x)[0][0] for x in list(fullres.name)]
fullresvlg = catalog.vlg[fullresndx]

volres = pd.read_csv('RESOLVE_inobssample.csv')
volresndx = [np.where(catalog.name == x)[0][0] for x in list(volres.name)]
volresvlg = catalog.vlg[volresndx]

elres = pd.read_csv('RESOLVE_full_hasnr5_dext_jhu.csv')
elresndx = [np.where(catalog.name == x)[0][0] for x in list(elres.name)]
elresvlg = catalog.vlg[elresndx]

selres = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv')
selresndx = [np.where(catalog.name == x)[0][0] for x in list(selres.name)]
selresvlg = catalog.vlg[selresndx]

fullresdf = pd.DataFrame(data = {'name': fullres.name, 'vlg': fullresvlg, 
                                 'rarad': fullres.radeg/180*np.pi})
fullresdf.to_csv('fullres_vlg_rarad.csv')

volresdf = pd.DataFrame(data = {'name': volres.name, 'vlg': volresvlg, 
                                 'rarad': volres.radeg/180*np.pi})
volresdf.to_csv('volres_vlg_rarad.csv')

elresdf = pd.DataFrame(data = {'name': elres.name, 'vlg': elresvlg, 
                                 'rarad': elres.radeg/180*np.pi})
elresdf.to_csv('elres_vlg_rarad_jhu.csv')

selresdf = pd.DataFrame(data = {'name': selres.name, 'vlg': selresvlg, 
                                 'rarad': selres.radeg/180*np.pi})
selresdf.to_csv('selres_vlg_rarad_jhu.csv')