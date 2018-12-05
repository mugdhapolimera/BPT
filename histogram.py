#histogram comparing completeness of data from ECO and resolve
#created by Carlynn Ferguson 4/3/2018

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

df = pd.read_pickle('RESOLVE.pkl')
res = np.array(df['logmstar'])
df = pd.read_pickle('ECO.pkl')
eco = np.array(df['logmstar'])

resvol = 52100.0 #Mpc^3
ecovol = 560800.0 #Mpc^3

(yres, xres, patches) = plt.hist(res, bins=10, color='m', linestyle= ':',
         histtype='step',normed=False, label='RESOLVE')

delta = xres[1]-xres[0]
xnew = xres[:10]+(delta/2)

(yeco, xeco, patches) = plt.hist(eco, bins=xres, color='c',
         histtype='step', normed=False, label='ECO')

print len(xnew)
print len(yres)

plt.figure(2)

plt.plot(xnew,yres/resvol, lw=3, color='m', label='RESOLVE')
plt.plot(xnew,yeco/ecovol, lw=3, color='c', label='ECO')

#plot lines at 9.1 and 10.1 to show dwarf galaxy region
plt.axvline(x=8.9, color='y', ls='--', lw=3, label = 'RESOLVE COMPLETENESS')

plt.legend(prop=dict(size=12))
plt.xlabel('Log Stellar Mass')
plt.ylabel('Per Survey Volume')

plt.show()
