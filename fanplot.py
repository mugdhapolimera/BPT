from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

res = pd.read_csv('RESOLVE_live13Oct2016.csv')
res_ext_t = Table.read('RESOLVE_SDSS_dext.fits')

res_ext = res_ext_t.to_pandas()
res['name'] = res['name'].str.strip() 
res_ext.columns.values[0] = 'name'
res_full = pd.merge(res, res_ext, on='name', how='left')

col = res_full[['name', 'radeg', 'dedeg', 'cz', 'dist2d_nn1', 'logmh', 'grpn', 'r90', 'r50']]

#targ = col[(col['dedeg'] < 5.)]
targsamp = col[(col['name'] == 'rf0250') | (col['name'] == 'rf0266') | (col['name'] == 'rs0463') | (col['name'] == 'rs1268') | (col['name'] == 'rs1398')]

cz = col[['cz']]
z = cz/3e5
radeg = col[['radeg']]
rahr = radeg/15
rarad = radeg * (np.pi / 180.)

cz_targ = targsamp[['cz']]
z_targ = cz_targ/3e5
radeg_targ = targsamp[['radeg']]
rahr_targ = radeg_targ / 15
rarad_targ = radeg_targ * (np.pi / 180.)

labels = ['rf0250', 'rf0266', 'rs0463', 'rs1268', 'rs1398']

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

vlgmin=min(resvlg[where(inobssample)])
vlgmax=max(resvlg[where(inobssample)]
vlgboundf=[vlgmin,vlgmax,fltarr(1000)+vlgmax,vlgmin,fltarr(1000)+vlgmin]
raradboundf=[330.,330.,findgen(1000)*(405.-330.)/999.+330.,405.,findgen(1000)*(330.-405.)/999.+405.]*!pi/180.
vlgbounds=[vlgmin,vlgmax,fltarr(1000)+vlgmax,vlgmin,fltarr(1000)+vlgmin]
raradbounds=[131.25,131.25,findgen(1000)*(236.25-131.25)/999.+131.25,236.25,findgen(1000)*(131.25-236.25)/999.+236.25]*!pi/180.
!psym=0
oplot,vlgboundf/70.,raradboundf,/polar
oplot,vlgbounds/70.,raradbounds,/polar
#plt.title('RESOLVE Fan Plot', y = 1.08)

#plt.show()
