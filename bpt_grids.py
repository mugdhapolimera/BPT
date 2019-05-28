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
from astropy.table import Table  # Used in converting to pandas DataFrame 
from marginalize_grid import marginalize

################################
#
# Read in and filter obs
#


#res_data = '../../../data/resolve_plusmetandemline-Hood.pkl'
res_data = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter_new.pkl'
res_den = pd.read_pickle(res_data)

obsR_nii = res_den['nii_6584_flux']
obsR_oiii = res_den['oiii_5007_flux']
obsR_sii_sum = res_den['sii_6717_flux'] + res_den['sii_6731_flux']
obsR_oi = res_den['oi_6300_flux']
obsR_h_alpha = res_den['h_alpha_flux']
obsR_h_beta = res_den['h_beta_flux']
obsR_heii_4686 = res_den['Flux_HeII_4685']

obsR_n2ha = np.log10(obsR_nii/obsR_h_alpha)
obsR_o3hb = np.log10(obsR_oiii/obsR_h_beta)
obsR_s2ha = np.log10(obsR_sii_sum/obsR_h_alpha)
obsR_o1ha = np.log10(obsR_oi/obsR_h_alpha)


################################
#
# Define demarcation functions
#

#def o3hbcomposite(log_NII_HA):
#    return (0.61 / (log_NII_HA - 0.05)) + 1.3
#
#def o3hbmain(log_NII_HA):
#    return (0.61 / (log_NII_HA - 0.47)) + 1.19
#
#def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
#    return 1.33 + (0.73 / (log_OI_HA + 0.59))
#def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
#    return -0.59
#def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
#    return 1.30 + (0.72 / (log_SII_HA - 0.32))

def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
    return 1.30 + (0.72 / (log_SII_HA - 0.32))
def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
    return 0.76 + 1.89*log_SII_HA
def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
    return 1.33 + (0.73 / (log_OI_HA + 0.59))
def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
    return 1.3 + 1.18*log_OI_HA
def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
    return -0.59
def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
    return 1.19 + (0.61 / (log_NII_HA - 0.47))

refn2ha = np.linspace(-3.0, 0.35)
#comp = o3hbcomposite(refn2ha)
#main = o3hbmain(refn2ha)

refsiiha = np.linspace(-2, 0.3)
main_sii = s2hamain(refsiiha)

refoiha = np.linspace(-3.0, -0.59)
main_oi = o1hamain(refoiha)

##################################
#
# Sims
#

#grid = 'L10'
grid = 'BPASS'
if grid == 'BPASS':
    sim = pd.read_csv('C:/Users/mugdhapolimera/github/izi/Richardson-0-0_1-0agn-BPASS-Binary-CSF-n=1e2-40.0Myr-NichollsCE-D_G-RR14_Fstar_0_3.csv')
    #sim = sim0
    #lowz, highz = marginalize(sim0)
    sim= sim[(sim["LOGQ"] > 6.9) & (sim["LOGQ"] < 8.9)]
    sim= sim[(sim["LOGZ"] > np.log10(0.1))]

if grid == 'L10':
    gridfile = r'C:/Users/mugdhapolimera/github/izi/l09_high_csf_n1e2_6.0Myr_new.fits'
    grid0 = Table.read(gridfile, format='fits')
    sim = grid0.to_pandas()

bpt = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/resolve_emlineclass_filter_new.csv')
bpt.index = bpt['galname']


##################################
#
# Set up color mapping for BPASS
#

metals = len(np.unique(sim['LOGZ']))
ionps = len(np.unique(sim['LOGQ']))
if grid =='BPASS' : 
    agn = np.unique(sim['AGNFRAC'])
else:
    agn = [0]

metal_colors = [plt.cm.BuPu(i) for i in np.linspace(0.15, 1, metals)]
u_colors = [plt.cm.YlOrRd(i) for i in np.linspace(0.15, 1, ionps)]
z50_colors = ['green']
def truncate_colormap(cmap, minval=0.15, maxval=1.0, n=256):
  	new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap
#metal_colors_map = truncate_colormap(cm.Blues, n = metals)
#u_colors_map = truncate_colormap(cm.Reds, n = ionps)

metal_colors_map = truncate_colormap(cm.BuPu, n = metals)
u_colors_map = truncate_colormap(cm.YlOrRd, n = ionps)

#gnuplot, rainbow, winter, spring, viridis, magma
for frac in [0,0.5,1]:#0.16,0.32,0.5,1]:
    if "oiii5007" in sim.keys(): #grid == 'L10':
        simdata = sim[sim['AGNFRAC']==frac]
        oiii = np.array(simdata["oiii5007"])
        hbeta = np.array(simdata["hbeta"])
        nii = np.array(simdata["nii6584"])
        halpha = np.array(simdata["halpha"])
        oi = np.array(simdata["oi6300"])
        sii = np.array(simdata["sii6717"] + simdata["sii6731"])
        ndx = np.where(simdata.LOGZ == np.unique(simdata.LOGZ)[6])[0]
        
    else:
        simdata = sim[sim['AGNFRAC']==frac]
        oiii = np.array(simdata[" 'O__3_500684A'"])
        hbeta = np.array(simdata[" 'H__1_486133A'"])
        nii = np.array(simdata[" 'N__2_658345A'"])
        halpha = np.array(simdata[" 'H__1_656281A'"])
        oi = np.array(simdata[" 'O__1_630030A'"])
        sii = np.array(simdata[" 'S__2_671644A'"] + simdata[" 'S__2_673082A'"])
       
    oiii = oiii.reshape((metals,ionps))
    hbeta = hbeta.reshape((metals,ionps))
    nii = nii.reshape((metals,ionps))
    halpha = halpha.reshape((metals,ionps))
    oi = oi.reshape((metals,ionps))
    sii = sii.reshape((metals,ionps))

    if grid == 'L10':
        ndx = sim.LOGQ[:ionps].argsort()
        oiii = oiii[:,ndx]
        hbeta = hbeta[:,ndx]
        nii = nii[:,ndx]
        halpha = halpha[:,ndx]
        oi = oi[:,ndx]
        sii = sii[:,ndx]
#    oiii = oiii.reshape((ionps,metals))
#    hbeta = hbeta.reshape((ionps,metals))
#    nii = nii.reshape((ionps,metals))
#    halpha = halpha.reshape((ionps,metals))
#    oi = oi.reshape((ionps,metals))
#    sii = sii.reshape((ionps,metals))
    
    ##################################
    #
    # Plot
    #
    
    #fig = plt.figure('AGN Fraction '+str(frac) + ' at 1.0 Z_solar')
    label_size = 8
    tick_label = 8
    cbar_tick_label = 8
    text_size = 8
    xmin, xmax = -2.50001, 0.30001
    ymin, ymax = -1.0, 1.2
    
    f, (sp1, sp2, sp3) = plt.subplots(1,3, sharey = True)
    
    sp1.scatter(obsR_n2ha,obsR_o3hb,c='k',s=0.5)
    #plt.scatter(obsR_n2ha[bpt['sftoagn']],obsR_o3hb[bpt['sftoagn']],color='r',marker='s')
    
    sp1.plot(refn2ha, n2hamain(refn2ha),'k')
    sp1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),'k--')
    
    x = np.log10(np.divide(nii,halpha))
    y = np.log10(np.divide(oiii,hbeta))
#    z50_x = x[8]
#    z50_y = y[8]
    
    sp1.set_prop_cycle(cycler('color',metal_colors))
    sp1.plot(np.transpose(x),np.transpose(y), lw = 2)
    #sp1.set_prop_cycle(cycler('color',z50_colors))
    #plt.plot(np.transpose(z50_x),np.transpose(z50_y), lw = 5)
    sp1.set_prop_cycle(cycler('color',u_colors))
    sp1.plot(x,y,ls='--')
    
    sp1.set_xlim(xmin,xmax)
    sp1.set_ylim(ymin,ymax)
    sp1.axes.tick_params(axis='both',which='both',top='on',right='on',
                         direction='in')
    #plt.xticks([-2.5,-1.5,-0.5])
    #plt.yticks([-1.0,0.0,1.0])
    #plt.xticks(fontsize=tick_label); plt.yticks(fontsize=tick_label)
    minorLocator = MultipleLocator(0.1)
    sp1.xaxis.set_minor_locator(minorLocator)
    minorLocator = MultipleLocator(0.1)
    sp1.yaxis.set_minor_locator(minorLocator)
    
    sp1.set_xlabel(r'$\rm log([$N II$]$ 6584 / H$\alpha)$',fontsize=15)
    sp1.set_ylabel(r'$\rm log([$O III$]$ 5007 / H$\beta$)',fontsize=15)
    
    Z = 10**np.unique(sim.LOGZ)
    sm = plt.cm.ScalarMappable(cmap=metal_colors_map,
                               norm=colors.Normalize(vmin=min(Z), vmax=max(Z)))
    sm._A = []
    smaxes = inset_axes(sp1, width=0.06, height=1.0, loc=4, 
            bbox_to_anchor=(0.1, 0.17), bbox_transform=sp1.figure.transFigure)
    cbar = plt.colorbar(sm,cax=smaxes)
    cbar.ax.set_title(r'Z / Z$_{\odot}$',fontsize=cbar_tick_label)
    #cbar.set_ticks(10**np.unique(simdata['LOGZ']))
    #cbar.set_ticklabels(np.around(10**np.unique(simdata['LOGZ']),1))
    cbar.ax.tick_params(labelsize=cbar_tick_label) 
    
    q = np.unique(sim.LOGQ)
    sm = plt.cm.ScalarMappable(cmap=u_colors_map,
                               norm=colors.Normalize(vmin=min(q), vmax=max(q)))
    sm._A = []
    smaxes = inset_axes(sp1, width=0.06, height=1.0, loc=4, 
            bbox_to_anchor=(0.18, 0.17), bbox_transform=sp1.figure.transFigure)
    cbar = plt.colorbar(sm,cax=smaxes)
    cbar.ax.set_title(r'log $q$',fontsize=cbar_tick_label)
    #cbar.set_ticks([-4.0,-0.5])
    #cbar.set_ticklabels([-4.0,-0.5])
    cbar.ax.tick_params(labelsize=cbar_tick_label) 
    
    
    #fig = plt.figure(2)
    
    label_size = 8
    tick_label = 8
    cbar_tick_label = 8
    text_size = 8
    xmin, xmax = -2.00001, 0.50001
    ymin, ymax = -1.0, 1.2
    
    #sp2 = plt.subplot(132)
    
    sp2.scatter(obsR_s2ha,obsR_o3hb,c='k',s=0.5)
    #plt.scatter(obsR_s2ha[bpt['sftoagn']],obsR_o3hb[bpt['sftoagn']],
    #            color='r',marker='s')
    sp2.plot(refsiiha, s2hamain(refsiiha),'k')
#    plt.plot(refsiiha-0.09,main_sii,'k')
    
    x = np.log10(np.divide(sii,halpha))
    y = np.log10(np.divide(oiii,hbeta))
    z50_x = x[8]
    z50_y = y[8]
    
    sp2.set_prop_cycle(cycler('color',metal_colors))
    sp2.plot(np.transpose(x),np.transpose(y), lw = 2)
    sp2.set_prop_cycle(cycler('color',z50_colors))
    #plt.plot(np.transpose(z50_x),np.transpose(z50_y), lw = 5)
    sp2.set_prop_cycle(cycler('color',u_colors))
    sp2.plot(x,y,ls='--')
    
    sp2.set_xlim(xmin,xmax)
    sp2.set_ylim(ymin,ymax)
    sp2.axes.tick_params(axis='both',which='both',top='on',right='on',
                         direction='in')
    #plt.xticks([-2.5,-1.5,-0.5, 0.5])
    #plt.yticks([-1.0,0.0,1.0])
    #plt.xticks(fontsize=tick_label); plt.yticks(fontsize=tick_label)
    minorLocator = MultipleLocator(0.1)
    sp2.xaxis.set_minor_locator(minorLocator)
    minorLocator = MultipleLocator(0.1)
    sp2.yaxis.set_minor_locator(minorLocator)
    
    sp2.set_xlabel(r'$\rm log([$S II$]$ 6717 + 6731/ H$\alpha)$',fontsize=15)
    #plt.ylabel(r'$[$O III$]$ 5007 / H$\beta$',fontsize=15)

    if frac == 0:
        sp2.set_title('(a)', y = -0.15)
    elif frac == 0.5:
        sp2.set_title('(b)', y = -0.15)
    elif frac == 1.0:
        sp2.set_title('(c)', y = -0.15)

    
    #fig = plt.figure(3)
    
    label_size = 8
    tick_label = 8
    cbar_tick_label = 8
    text_size = 8
    xmin, xmax = -3.00001, 0.00001
    ymin, ymax = -1.0, 1.2
    
    #sp3 = plt.subplot(133)
    
    sp3.scatter(obsR_o1ha,obsR_o3hb,c='k',s=0.5)
#    plt.scatter(obsR_o1ha[bpt['sftoagn']],obsR_o3hb[bpt['sftoagn']],
#                color='r',marker='s')
    sp3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),'k')
    
    x = np.log10(np.divide(oi,halpha))
    y = np.log10(np.divide(oiii,hbeta))
    z50_x = x[8]
    z50_y = y[8]
    
    sp3.set_prop_cycle(cycler('color',metal_colors))
    sp3.plot(np.transpose(x),np.transpose(y), lw = 2)
    sp3.set_prop_cycle(cycler('color',z50_colors))
    #plt.plot(np.transpose(z50_x),np.transpose(z50_y), lw = 5)
    sp3.set_prop_cycle(cycler('color',u_colors))
    sp3.plot(x,y,ls='--')
    
    sp3.set_xlim(xmin,xmax)
    sp3.set_ylim(ymin,ymax)
    #sp1.axes.tick_params(axis='both',which='both',top='on',right='on',direction='in')
    #plt.xticks([-3.0, -2.5,-1.5,-0.5])
    #plt.yticks([-1.0,0.0,1.0])
    #plt.xticks(fontsize=tick_label); plt.yticks(fontsize=tick_label)
    minorLocator = MultipleLocator(0.1)
    sp3.xaxis.set_minor_locator(minorLocator)
    minorLocator = MultipleLocator(0.1)
    sp3.yaxis.set_minor_locator(minorLocator)
    
    sp3.set_xlabel(r'$\rm log ([$O I$]$ 6300 / H$\alpha$)',fontsize=15)
    #plt.ylabel(r'$[$O III$]$ 5007 / H$\beta$',fontsize=15)
    sp3.legend(title = 'AGN Fraction '+str(int(frac*100))+'%', loc = 'upper right')

    
#plt.savefig("Z_U-binary-agn0000.png",dpi=600)


















