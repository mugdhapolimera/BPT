import numpy as np
import pandas as pd
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
res_data = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5.pkl'
res_den = pd.read_pickle(res_data)
full = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5.csv')
full.index = full.name
port = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5_port.csv')
port.index = port.name
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

def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
    return 1.19 + (0.61 / (log_NII_HA - 0.47))
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

#comp = o3hbcomposite(refn2ha)
#main = o3hbmain(refn2ha)

refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, -0.4)
refsiiha = np.linspace(-2, 0.3,100)
main_sii = s2hamain(refsiiha)
main_oi = o1hamain(refoiha)

##################################
#
# Sims
#

#grid = 'L10'
grid = 'BPASS'
qfix = 0
if grid == 'BPASS':
    #CSF
    sim = pd.read_csv('C:/Users/mugdhapolimera/github/izi/Richardson-0-0_1-0agn-BPASS-Binary-CSF-n=1e2-40.0Myr-NichollsCE-D_G-RR14_Fstar_0_3.csv')
    
    #SSP
    #sim = pd.read_csv('C:/Users/mugdhapolimera/github/izi/Richardson-0-0_1-0agn-BPASS-Binary-SSP-n=1e2-1.0Myr-NichollsCE-D_G-RR14_Fstar_0_3.csv')
                      
    #Richardson-0-0_1-0agn-M-5_0-BPASS-Binary-CSF-n=1e3-40.0Myr-NichollsCE-D_G-RR14_Fstar_0_3-unified.csv')
    
    #sim = sim0
    #lowz, highz = marginalize(sim0)
    if qfix : 
        sim= sim[(sim["LOGQ"] == 7.227)] #> 6.9) & (sim["LOGQ"] < 8.9)]
        sim= sim[(sim["LOGZ"] > np.log10(0.1))] #== np.log10(0.3))] 
    else:
        sim= sim[(sim["LOGQ"] > 6.9) & (sim["LOGQ"] < 8.9)]
        sim= sim[(sim["LOGZ"] == np.log10(0.3))] 

if grid == 'L10':
    gridfile = r'C:/Users/mugdhapolimera/github/izi/l09_high_csf_n1e2_6.0Myr_new.fits'
    grid0 = Table.read(gridfile, format='fits')
    sim = grid0.to_pandas()

bpt = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/resolve_emlineclass_full_snr5.csv')
bpt.index = bpt['galname']

#rms = np.sqrt(np.mean(np.log10(res_den.oi_6300_flux_err[bpt['defstarform']]/res_den.oi_6300_flux[bpt['defstarform']]))**2)

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

metal_list = [0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1., 1.5, 2. ]
ionps_list = [6.9 , 7.2, 7.5, 7.7, 8. , 8.2, 8.5, 8.7]
f, (sp1, sp2, sp3) = plt.subplots(1,3, sharey = True)
fid = 1  
frac_color = {0: 'b', 0.5: 'fuchsia', 1:'brown'}
gals = nov_target#['rf0006', 'rf0073','rf0503']#,'rs0463']#, 'rs0814']
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
#        ndx = np.where(simdata.LOGZ == np.unique(simdata.LOGZ)[6])[0]
        
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
    if qfix: 
        ndx = np.where(simdata.LOGZ == np.log10(0.4))[0][0]
    else:
        ndx = np.where(simdata.LOGQ == 7.227)[0][0]
    if grid == 'L10':
        ndx = sim.LOGQ[:ionps].argsort()
        oiii = oiii[:,ndx]
        hbeta = hbeta[:,ndx]
        nii = nii[:,ndx]
        halpha = halpha[:,ndx]
        oi = oi[:,ndx]
        sii = sii[:,ndx]
    ##################################
    #
    # Plot
    #
    
    #fig = plt.figure('AGN Fraction '+str(frac) + ' at 1.0 Z_solar')
    label_size = 8
    tick_label = 8
    cbar_tick_label = 8
    text_size = 8
    xmin, xmax = -2.00001, 0.30001
    ymin, ymax = -1.0, 1.2
    
    if frac == 0:
        
        sp1.scatter(obsR_n2ha,obsR_o3hb,c='k',s=0.5, alpha = 0.3)

        
    #plt.scatter(obsR_n2ha[bpt['sftoagn']],obsR_o3hb[bpt['sftoagn']],color='r',marker='s')
    
    sp1.plot(refn2ha, n2hamain(refn2ha),'k')
    sp1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),'k--')
    
    x = np.log10(np.divide(nii,halpha))
    y = np.log10(np.divide(oiii,hbeta))
    if fid:
        if x.shape[0] == 1:
            x_fid = x[0,ndx]
            y_fid = y[0,ndx]
        else:
            x_fid = x[ndx,0]
            y_fid = y[ndx,0]

#    z50_x = x[8]
#    z50_y = y[8]
    
    sp1.set_prop_cycle(cycler('color',metal_colors))
#    sp1.plot(np.transpose(x),np.transpose(y), lw = 2)
    sp1.set_prop_cycle(cycler('color',z50_colors))
#    plt.plot(np.transpose(z50_x),np.transpose(z50_y), lw = 5)
    sp1.set_prop_cycle(cycler('color',u_colors))
#    sp1.plot(x,y,ls='--')
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
    if fid:
        sp1.scatter(x_fid,y_fid,marker = 's', s = 150, c = frac_color[frac])
    
    if frac == 1:
        sp1.set_xlabel(r'$\rm log([$N II$]$ 6584 / H$\alpha)$',fontsize=15)
    sp1.set_ylabel(r'$\rm log([$O III$]$ 5007 / H$\beta$)',fontsize=15)
    
    if metals > 1:
        Z = 10**np.unique(sim.LOGZ)
        sm = plt.cm.ScalarMappable(cmap=metal_colors_map,
                                   norm=colors.Normalize(vmin=min(Z), vmax=max(Z)))
        sm._A = []
        smaxes = inset_axes(sp1, width=0.06, height=1.0, loc=4, 
                bbox_to_anchor=(0.1, 0.17), bbox_transform=sp1.figure.transFigure)
        cbar = plt.colorbar(sm,cax=smaxes)
        cbar.ax.set_title(r'Z / Z$_{\odot}$',fontsize=cbar_tick_label)
        cbar.set_ticks(np.linspace(0.15,2.0, metals))
        cbar.set_ticklabels(metal_list)
        cbar.ax.tick_params(labelsize=cbar_tick_label) 
        
        q = np.unique(sim.LOGQ)
        sm = plt.cm.ScalarMappable(cmap=u_colors_map,
                                   norm=colors.Normalize(vmin=min(q), vmax=max(q)))
        sm._A = []
        smaxes = inset_axes(sp1, width=0.06, height=1.0, loc=4, 
                bbox_to_anchor=(0.18, 0.17), bbox_transform=sp1.figure.transFigure)
        cbar = plt.colorbar(sm,cax=smaxes)
        cbar.ax.set_title(r'log $q$',fontsize=cbar_tick_label)
        cbar.set_ticks(np.linspace(6.977,8.727, ionps))
        cbar.set_ticklabels(ionps_list)
        cbar.ax.tick_params(labelsize=cbar_tick_label) 
    
    
    #fig = plt.figure(2)
    
    label_size = 8
    tick_label = 8
    cbar_tick_label = 8
    text_size = 8
    xmin, xmax = -2.00001, 0.50001
    ymin, ymax = -1.5, 1.2
    
    #sp2 = plt.subplot(132)
    if frac == 0:
        
        sp2.scatter(obsR_s2ha,obsR_o3hb,c='k',s=1, alpha = 0.3)
    #plt.scatter(obsR_s2ha[bpt['sftoagn']],obsR_o3hb[bpt['sftoagn']],
    #            color='r',marker='s')
    sp2.plot(refsiiha, s2hamain(refsiiha),'k')
    sp2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
                  'k-.')

#    plt.plot(refsiiha-0.09,main_sii,'k')
    
    x = np.log10(np.divide(sii,halpha)) #+ 0.02
    y = np.log10(np.divide(oiii,hbeta)) #+ 0.04
    if fid:
        if x.shape[0] == 1:
            x_fid = x[0,ndx]
            y_fid = y[0,ndx]
        else:
            x_fid = x[ndx,0]
            y_fid = y[ndx,0]

#    z50_x = x[8]
#    z50_y = y[8]
    
#    cmap = plt.get_cmap('rainbow',len(simdata))
    #sp2.set_prop_cycle(cycler('color',metal_colors))
#    sp2.plot(np.transpose(x),np.transpose(y), 'g', lw = 2)
    #sp2.plot(x,y, 'g', lw = 2)
    
#    if x.shape[0] == 1:
#        sp2.scatter(x[0],y[0],marker = 's', c = np.array(simdata.LOGQ), 
#                cmap = cmap, s = 100)
#    else:
#        sp2.scatter(x[:,0],y[:,0],marker = 's', c = np.array(simdata.LOGZ), 
#                cmap = cmap, s = 100)
#    color = ['red','orange','yellow','green','blue','purple','magenta','yellowgreen',
#         'orchid']
#    if x.shape[0] == 1:
#        for i in range(len(simdata)-1):
#            sp2.plot(x[0][i:i+2],y[0][i:i+2], c = color[i], lw = 2)
#    else:
#        for i in range(len(simdata)-1):
#            sp2.plot(x[i:i+2],y[i:i+2], c = color[i], lw = 2)
    #    sp2.set_prop_cycle(cycler('color',z50_colors))
#    plt.plot(np.transpose(z50_x),np.transpose(z50_y), lw = 5)
    #sp2.set_prop_cycle(cycler('color',u_colors))
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
    if fid:
        sp2.scatter(x_fid,y_fid,marker = 's', s = 150, c = frac_color[frac])
    
    if frac == 1:
        sp2.set_xlabel(r'$\rm log([$S II$]$ 6717 + 6731/ H$\alpha)$',fontsize=15)
    #plt.ylabel(r'$[$O III$]$ 5007 / H$\beta$',fontsize=15)

#    if frac == 0:
#        sp2.set_title('(a)', y = -0.15)
#    elif frac == 0.5:
#        sp2.set_title('(b)', y = -0.15)
#    elif frac == 1.0:
#        sp2.set_title('(c)', y = -0.15)

    
    #fig = plt.figure(3)
    
    label_size = 8
    tick_label = 8
    cbar_tick_label = 8
    text_size = 8
    xmin, xmax = -2.50001, -0.40001
    ymin, ymax = -1.0, 1.00001
    
    #sp3 = plt.subplot(133)
    if frac == 0:
        sp3.scatter(obsR_o1ha,obsR_o3hb,c='k',s=1, alpha = 0.3)
#    plt.scatter(obsR_o1ha[bpt['sftoagn']],obsR_o3hb[bpt['sftoagn']],
#                color='r',marker='s')
    sp3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),'k')
    sp3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
                               'k-.')

    x = np.log10(np.divide(oi,halpha)) + 0.07
    y = np.log10(np.divide(oiii,hbeta)) + 0.04
    if fid:
        if x.shape[0] == 1:
            x_fid = x[0,ndx]
            y_fid = y[0,ndx]
        else:
            x_fid = x[ndx,0]
            y_fid = y[ndx,0]

#    z50_x = x[8]
#    z50_y = y[8]
    
    sp3.set_prop_cycle(cycler('color',metal_colors))
#    sp3.plot(np.transpose(x),np.transpose(y), lw = 2)
#    sp3.set_prop_cycle(cycler('color',z50_colors))
    #plt.plot(np.transpose(z50_x),np.transpose(z50_y), lw = 5)
    sp3.set_prop_cycle(cycler('color',u_colors))
#    sp3.plot(x,y,ls='--')
   
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
    if fid:
        sp3.scatter(x_fid,y_fid,marker = 's', s = 150, c = frac_color[frac])
    
    if frac == 1:
        sp3.set_xlabel(r'$\rm log ([$O I$]$ 6300 / H$\alpha$)',fontsize=15)
    #plt.ylabel(r'$[$O III$]$ 5007 / H$\beta$',fontsize=15)
    #sp3.legend(title = 'AGN Fraction '+str(int(frac*100))+'%', loc = 'upper right')

for gal in gals:
    plotflag = 0
    print gal
    if gal in nsaflag.index.values:
        if nsaflag.loc[gal]['sftoagn']:
            print 'NSA'
            niiha = np.log10(nsa.loc[gal]['nii_6584_flux']/nsa.loc[gal]['h_alpha_flux'])
            SII = nsa.loc[gal]['sii_6731_flux']
            siiha = np.log10(SII/nsa.loc[gal]['h_alpha_flux'])
            oiha = np.log10(nsa.loc[gal]['oi_6300_flux']/nsa.loc[gal]['h_alpha_flux'])
            oiiihb = np.log10(nsa.loc[gal]['oiii_5007_flux']/nsa.loc[gal]['h_beta_flux'])
            plotflag = 1
    elif (gal in portflag.index.values) & (gal != 'rs0472'):
        if portflag.loc[gal]['sftoagn']:
            print 'Port'
            niiha = np.log10(port.loc[gal]['Flux_NII_6583']/port.loc[gal]['Flux_Ha_6562'])
            SII = port.loc[gal]['Flux_SII_6716'] + port.loc[gal]['Flux_SII_6730']
            siiha = SII/port.loc[gal]['Flux_Ha_6562']
            oiiihb = np.log10(port.loc[gal]['Flux_OIII_5006']/port.loc[gal]['Flux_Hb_4861'])
            oiha = np.log10(port.loc[gal]['Flux_OI_6300']/port.loc[gal]['Flux_Ha_6562'])
            plotflag = 1
    if ~plotflag & (gal in jhuflag.index.values):
        if jhuflag.loc[gal]['sftoagn']:
            print 'JHU'
            niiha = np.log10(full.loc[gal]['nii_6584_flux']/full.loc[gal]['h_alpha_flux'])
            SII = full.loc[gal]['sii_6717_flux'] + full.loc[gal]['sii_6731_flux']
            siiha = np.log10(SII/full.loc[gal]['h_alpha_flux'])
            oiha = np.log10(full.loc[gal]['oi_6300_flux']/full.loc[gal]['h_alpha_flux'])
            oiiihb = np.log10(full.loc[gal]['oiii_5007_flux']/full.loc[gal]['h_beta_flux'])
    print niiha, siiha, oiha, oiiihb
    sp1.scatter(niiha,oiiihb, c='r',marker = '*', s= 200, zorder = 1)#marker = (6,2,0), s=200)
    sp2.scatter(siiha,oiiihb, c='r',marker = '*', s= 200)#marker = (6,2,0), s=200)
    sp3.scatter(oiha,oiiihb, c='r',marker = '*', s= 200, zorder = 1)#marker = (6,2,0), s=200)
sp2.scatter(siiha,oiiihb, c='r',marker = '*', s= 200, label = 'Targets for this proposal')
#    for i, txt in enumerate([gal]):
#        sp3.annotate(txt, (oiha, oiiihb))#marker = (6,2,0), s=200)
#niiha = np.log10(nsa.loc[sami]['nii_6584_flux']/nsa.loc[sami]['h_alpha_flux'])
#SII = nsa.loc[sami]['sii_6731_flux'] #+ nsa.loc[sami]['sii_6717_flux']
#siiha = np.log10(SII/nsa.loc[sami]['h_alpha_flux'])
#oiha = np.log10(nsa.loc[sami]['oi_6300_flux']/nsa.loc[sami]['h_alpha_flux'])
#oiiihb = np.log10(nsa.loc[sami]['oiii_5007_flux']/nsa.loc[sami]['h_beta_flux'])
#sp1.scatter(niiha,oiiihb, c='c',marker = '*', s= 100)#marker = (6,2,0), s=200)
#sp2.scatter(siiha,oiiihb, c='c',marker = '*', s= 100)#marker = (6,2,0), s=200)
#sp3.scatter(oiha,oiiihb, c='c',marker = '*', s= 100)#marker = (6,2,0), s=200)
niiha = np.log10(full.loc[sami]['nii_6584_flux']/full.loc[sami]['h_alpha_flux'])
SII = full.loc[sami]['sii_6731_flux'] + full.loc[sami]['sii_6717_flux']
siiha = np.log10(SII/full.loc[sami]['h_alpha_flux'])
oiha = np.log10(full.loc[sami]['oi_6300_flux']/full.loc[sami]['h_alpha_flux'])
oiiihb = np.log10(full.loc[sami]['oiii_5007_flux']/full.loc[sami]['h_beta_flux'])
sp1.scatter(niiha,oiiihb, c='lime',marker = '*', s= 100, edgecolor = 'k')
sp2.scatter(siiha,oiiihb, c='lime',marker = '*', s= 100, 
            label = 'SFing-AGN in SAMI', edgecolor = 'k')
sp3.scatter(oiha,oiiihb, c='lime',marker = '*', s= 100, edgecolor = 'k')
sp2.legend(loc = 'lower left', fontsize = 15)
#for i, txt in enumerate(sami):
#    sp3.annotate(txt, (oiha[i]+0.05, oiiihb[i]))#marker = (6,2,0), s=200)

#gal = 'rs0909'
#sp1.scatter(np.log10(nsa.loc[gal]['nii_6584_flux']/nsa.loc[gal]['h_alpha_flux']),
#                    np.log10(nsa.loc[gal]['oiii_5007_flux']/nsa.loc[gal]['h_beta_flux']),
#                             c='c',marker = '*', s= 100)#marker = (6,2,0), s=200)
#SII = nsa.loc[gal]['sii_6731_flux']
#sp2.scatter(np.log10(SII/nsa.loc[gal]['h_alpha_flux']),
#            np.log10(nsa.loc[gal]['oiii_5007_flux']/nsa.loc[gal]['h_beta_flux']),
#                     c='c',marker = '*', s= 100)#marker = (6,2,0), s=200)
#sp3.scatter(np.log10(nsa.loc[gal]['oi_6300_flux']/nsa.loc[gal]['h_alpha_flux']),
#            np.log10(nsa.loc[gal]['oiii_5007_flux']/nsa.loc[gal]['h_beta_flux']),
#                     c='c',marker = '*', s= 100)#marker = (6,2,0), s=200)
#gal = 'rs0323'
#sp1.scatter(np.log10(nsa.loc[gal]['nii_6584_flux']/nsa.loc[gal]['h_alpha_flux']),
#                    np.log10(nsa.loc[gal]['oiii_5007_flux']/nsa.loc[gal]['h_beta_flux']),
#                             c='b',marker = '*', s= 100)#marker = (6,2,0), s=200)
#SII = nsa.loc[gal]['sii_6731_flux']
#sp2.scatter(np.log10(SII/nsa.loc[gal]['h_alpha_flux']),
#            np.log10(nsa.loc[gal]['oiii_5007_flux']/nsa.loc[gal]['h_beta_flux']),
#                     c='b',marker = '*', s= 100)#marker = (6,2,0), s=200)
#sp3.scatter(np.log10(nsa.loc[gal]['oi_6300_flux']/nsa.loc[gal]['h_alpha_flux']),
#            np.log10(nsa.loc[gal]['oiii_5007_flux']/nsa.loc[gal]['h_beta_flux']),
#                     c='b',marker = '*', s= 100)#marker = (6,2,0), s=200)

#plt.savefig("Z_U-binary-agn0000.png",dpi=600)

#o1ha = np.array(np.log10(full.loc[target]['oi_6300_flux']/full.loc[target]['h_alpha_flux']))
#o3hb = np.array(np.log10(full.loc[target]['oiii_5007_flux']/full.loc[target]['h_beta_flux']))
#for i, txt in enumerate(target):
#    #print i,txt
#    sp3.annotate(txt, (o1ha[i], o3hb[i]))#marker = (6,2,0), s=200)














