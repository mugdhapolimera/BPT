import numpy as np
import pandas as pd
import csv
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cycler import cycler


Zs = ['agn_0000-Z_0050','agn_0000-Z_0100','agn_0000-Z_0150','agn_0000-Z_0200','agn_0000-Z_0300','agn_0000-Z_0400',\
  'agn_0000-Z_0500','agn_0000-Z_0700','agn_0000-Z_1000','agn_0000-Z_1500','agn_0000-Z_2000']


################################
#
# Read in and filter obs
#


res_data = '../../../data/resolve_plusmetandemline-Hood.pkl'
res_den = pd.read_pickle(res_data)

obsR_nii = res_den['nii_6584_flux_ext']
obsR_oiii = res_den['oiii_5007_flux_ext']
obsR_sii_sum = res_den['sii_6717_flux_ext'] + res_den['sii_6731_flux_ext']
obsR_oi = res_den['oi_6300_flux_ext']
obsR_h_alpha = res_den['h_alpha_flux_ext']
obsR_h_beta = res_den['h_beta_flux_ext']
obsR_heii_4686 = res_den['Flux_HeII_4685_ext']

obsR_n2ha = np.log10(obsR_nii/obsR_h_alpha)
obsR_o3hb = np.log10(obsR_oiii/obsR_h_beta)
obsR_s2ha = np.log10(obsR_sii_sum/obsR_h_alpha)
obsR_o1ha = np.log10(obsR_oi/obsR_h_alpha)


################################
#
# Define demarcation functions
#

def o3hbcomposite(log_NII_HA):
    return (0.61 / (log_NII_HA - 0.05)) + 1.3

def o3hbmain(log_NII_HA):
    return (0.61 / (log_NII_HA - 0.47)) + 1.19

refn2ha = np.linspace(-3.0, 0.3)
main = o3hbmain(refn2ha)
refn2ha = np.linspace(-3.0, 0.0)
comp = o3hbcomposite(refn2ha)




##################################
#
# Set up color mapping for BPASS
#

metals = 11
ionps = 15

metal_colors = [plt.cm.Blues(i) for i in np.linspace(0.15, 1, metals)]
u_colors = [plt.cm.Reds(i) for i in np.linspace(0.15, 1, ionps)]

def truncate_colormap(cmap, minval=0.15, maxval=1.0, n=100):
  	new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap

metal_colors_map = truncate_colormap(cm.Blues)
u_colors_map = truncate_colormap(cm.Reds)



##################################
#
# Sims
#

simdata = np.genfromtxt(Zs[0]+'.lin', dtype=float, delimiter='\t', names=True)

oiii = simdata['O__3_500684A']
hbeta = simdata['H__1_486133A']
nii = simdata['N__2_658345A']
halpha = simdata['H__1_656281A']

for i in range(1,len(Zs)):
	tmpdata = np.genfromtxt(Zs[i]+'.lin', dtype=float, delimiter='\t', names=True)
	oiii = np.concatenate((oiii,tmpdata['O__3_500684A']))
	hbeta = np.concatenate((hbeta,tmpdata['H__1_486133A']))
	nii = np.concatenate((nii,tmpdata['N__2_658345A']))
	halpha = np.concatenate((halpha,tmpdata['H__1_656281A']))

oiii = np.reshape(oiii,(metals,ionps))
hbeta = np.reshape(hbeta,(metals,ionps))
nii = np.reshape(nii,(metals,ionps))
halpha = np.reshape(halpha,(metals,ionps))




##################################
#
# Plot
#

fig = plt.figure(1)

label_size = 8
tick_label = 8
cbar_tick_label = 8
text_size = 8
xmin, xmax = -2.50001, 0.00001
ymin, ymax = -1.0, 1.5

sp1 = plt.subplot(111)

plt.scatter(obsR_n2ha,obsR_o3hb,c='k',s=0.5)

plt.plot(refn2ha,main,'k')
plt.plot(refn2ha,comp,'k--')

x = np.log10(np.divide(nii,halpha))
y = np.log10(np.divide(oiii,hbeta))

sp1.set_prop_cycle(cycler('color',metal_colors))
plt.plot(np.transpose(x),np.transpose(y))
sp1.set_prop_cycle(cycler('color',u_colors))
plt.plot(x,y,ls='--')

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
sp1.axes.tick_params(axis='both',which='both',top='on',right='on',direction='in')
plt.xticks([-2.5,-1.5,-0.5])
plt.yticks([-1.0,0.0,1.0])
plt.xticks(fontsize=tick_label); plt.yticks(fontsize=tick_label)
minorLocator = MultipleLocator(0.1)
sp1.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(0.1)
sp1.yaxis.set_minor_locator(minorLocator)

plt.xlabel(r'$[$N II$]$ 6584 / H$\alpha$',fontsize=label_size)
plt.ylabel(r'$[$O III$]$ 5007 / H$\beta$',fontsize=label_size)

sm = plt.cm.ScalarMappable(cmap=metal_colors_map,norm=colors.Normalize(vmin=0.05, vmax=2.0))
sm._A = []
smaxes = inset_axes(sp1, width=0.06, height=0.4, loc=3, bbox_to_anchor=(0.18, 0.24), bbox_transform=sp1.figure.transFigure)
cbar = plt.colorbar(sm,cax=smaxes)
cbar.ax.set_title(r'Z / Z$_{\odot}$',fontsize=cbar_tick_label)
cbar.set_ticks([0.05,2.0])
cbar.set_ticklabels([0.05,2.0])
cbar.ax.tick_params(labelsize=cbar_tick_label) 

sm = plt.cm.ScalarMappable(cmap=u_colors_map,norm=colors.Normalize(vmin=-4.0, vmax=-0.5))
sm._A = []
smaxes = inset_axes(sp1, width=0.06, height=0.4, loc=3, bbox_to_anchor=(0.30, 0.24), bbox_transform=sp1.figure.transFigure)
cbar = plt.colorbar(sm,cax=smaxes)
cbar.ax.set_title(r'log $U$',fontsize=cbar_tick_label)
cbar.set_ticks([-4.0,-0.5])
cbar.set_ticklabels([-4.0,-0.5])
cbar.ax.tick_params(labelsize=cbar_tick_label) 


plt.savefig("Z_U-binary-agn0000.png",dpi=600)


















