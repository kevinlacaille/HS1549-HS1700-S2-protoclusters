import numpy as np
import matplotlib.pyplot as pl
import random as r

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def kcorr(z): #Barger+2014
    return 19.88 + 3.2*np.log10(z) + 1.13*(np.log10(z))**2 + 2.79*(np.log10(z))**3 + 2.58*(np.log10(z))**4

def fill_between(x, y1, y2=0, ax=None, **kwargs):
    """
        Plot filled region between `y1` and `y2`.

        This function works exactly the same as matplotlib's fill_between, except
        that it also plots a proxy artist (specifically, a rectangle of 0 size)
        so that it can be added it appears on a legend.
        """
    ax = ax if ax is not None else pl.gca()
    ax.fill_between(x, y1, y2, **kwargs)
    p = pl.Rectangle((0, 0), 0, 0, **kwargs)
    ax.add_patch(p)
    return p


#Import template
template = np.loadtxt('../SED/ALMA_compositeSED.dat',unpack=True)
wavelength_1 = template[0]
S_1 = template[1]# * 2.998e8 / wavelength_1**2# / max(template[1])
template = np.loadtxt('SED/ALESS_compositeSED_JMS2015_LBOL.dat',unpack=True)
wavelength_2= template[0]*1e-4
S_2 = template[1]# * 2.998e8 / (wavelength_2)**2# / max(template[1])
template = np.loadtxt('SED/smm2135_sed.dat',unpack=True)
wavelength_3= template[0]/(1+2.3259)
S_3 = template[1] / max(template[1])
pl.close()
'''
pl.figure()
pl.xlabel(r'$\lambda$ ($\mu$m)')
pl.ylabel(r'$S_{\nu}$')
# pl.xlim(0.05,1e4)
# pl.ylim(1e-6,1.2)
pl.xlim(0.05,1e4)
pl.ylim(1e-6,200)
#pl.loglog(wavelength_1,S_1)#,label='ALMA composite')
pl.loglog(wavelength_2,S_2, label='z=0')
pl.loglog(wavelength_2*(1+2),S_2, label='z=2')
#pl.loglog(wavelength_1*(1+2.3),S_1,label='ALMA composite')
# pl.loglog(wavelength_2*(1+2.3),S_2, label='ALESS composite')
# pl.loglog(wavelength_3,S_3, label='Eyelash')
pl.axvline(x=3.6,ls='-',c='k', label='Spitzer-IRAC band')
pl.axvline(x=4.5,ls='-',c='k')
pl.axvline(x=5.6,ls='-',c='k')
pl.axvline(x=8.0,ls='-',c='k')
pl.legend()
pl.savefig('ALMA_ALESS_SED.pdf',bbox_inches='tight')
pl.close()
'''

#import ALESS SMGs Simpson 2014
cat_ALESS = open('ALESS_SMGs_Simpson_2014.txt','r')
ID_ALESS = []
ch4_ALESS = []
e_ch4_ALESS = []
ch3_ALESS = []
e_ch3_ALESS = []
ch2_ALESS = []
e_ch2_ALESS = []
ch1_ALESS = []
e_ch1_ALESS = []
for line in cat_ALESS.readlines():
    tmp = line.split()
    ID_ALESS.append(tmp[0]+'_'+tmp[1])
    ch4_ALESS.append(float(tmp[-2]))
    e_ch4_ALESS.append(float(tmp[-1]))
    ch3_ALESS.append(float(tmp[-4]))
    e_ch3_ALESS.append(float(tmp[-3]))
    ch2_ALESS.append(float(tmp[-6]))
    e_ch2_ALESS.append(float(tmp[-5]))
    ch1_ALESS.append(float(tmp[-8]))
    e_ch1_ALESS.append(float(tmp[-7]))
cat_ALESS.close()
ID_ALESS = np.array(ID_ALESS)
ch4_ALESS = np.array(ch4_ALESS)
e_ch4_ALESS = np.array(e_ch4_ALESS)
ch3_ALESS = np.array(ch3_ALESS)
e_ch3_ALESS = np.array(e_ch3_ALESS)
ch2_ALESS = np.array(ch2_ALESS)
e_ch2_ALESS = np.array(e_ch2_ALESS)
ch1_ALESS = np.array(ch1_ALESS)
e_ch1_ALESS = np.array(e_ch1_ALESS)

S_ch4_ALESS = 10**(ch4_ALESS/2.5)
S_ch3_ALESS = 10**(ch3_ALESS/2.5)
S_ch2_ALESS = 10**(ch2_ALESS/2.5)
S_ch1_ALESS = 10**(ch1_ALESS/2.5)
e_S_ch4_ALESS = 10**(e_ch4_ALESS/2.5)
e_S_ch3_ALESS = 10**(e_ch3_ALESS/2.5)
e_S_ch2_ALESS = 10**(e_ch2_ALESS/2.5)
e_S_ch1_ALESS = 10**(e_ch1_ALESS/2.5)

z_ALESS = np.zeros(len(ID_ALESS))
cat_ALESS_z = open('ALESS_SMGs_Simpson_2014_z.txt','r')
ID_z_ALESS = []
# e_z_up_ALESS = []
# e_z_low_ALESS = []
for line in cat_ALESS_z.readlines():
    tmp = line.split()
    ID_z_ALESS.append(tmp[0]+'_'+tmp[1])
    z_ALESS[np.where(tmp[0]+'_'+tmp[1] == ID_ALESS)[0]]+=float(tmp[2])
cat_ALESS_z.close()
ID_z_ALESS = np.array(ID_z_ALESS)
z_ALESS = np.array(z_ALESS)
z_index_ALESS = np.where(z_ALESS>0)[0]
z_index_ALESS_1700 = np.where((z_ALESS>2.2) & (z_ALESS<2.4))[0]
# e_z_up_ALESS = np.array(e_z_up_ALESS)
# e_z_low_ALESS = np.array(e_z_low_ALESS)


CAT_SPIRE = open('../../h1700/Herschel_SPIRE/1700_SPIRE_fluxes.gaia','r')
RMS_250 = np.sqrt(6.7**2+3**3)*1e-3
RMS_350 = np.sqrt(7**2+2**2)*1e-3
RMS_500 = np.sqrt(7.1**2+2.6**2)*1e-3
ID = []
S_250 = []
S_350 = []
S_500 = []
#confusion + instrumental noise from Kato+2016: https://arxiv.org/pdf/1605.07370.pdf
for line in CAT_SPIRE.readlines()[2:]:
    tmp = line.split()
    ID.append(tmp[0])
    if float(tmp[1]) <= RMS_250:
        S_250.append(-1)
    else:
        S_250.append(float(tmp[1]))

    if float(tmp[2]) <= RMS_350:
        S_350.append(-1)
    else:
        S_350.append(float(tmp[2]))

    if float(tmp[3]) <= RMS_500:
        S_500.append(-1)
    else:
        S_500.append(float(tmp[3]))
S_250 = np.array(S_250)
e_S_250 = RMS_250*np.ones(len(ID))
S_350 = np.array(S_350)
e_S_350 = RMS_350*np.ones(len(ID))
S_500 = np.array(S_500)
e_S_500 = RMS_500*np.ones(len(ID))
ID = np.array(ID)
CAT_SPIRE.close()


#open file
cat_IR_1700 = open('../../h1700/mutli_wavelength_fixedmaybe.gaia','r')
#import 850um data
ID_1700 = []
ch4_1700 = []
ch3_1700 = []
ch2_1700 = []
ch1_1700 = []
K_1700 = []
z_1700 = []
for line in cat_IR_1700.readlines()[2:]:
    tmp = line.split()
    ID_1700.append(tmp[0])
    z_1700.append(float(tmp[-1]))
    # if tmp[0] == 'h1700_1':
    #     z_1700.append(2.816)
    # if tmp[0] == 'h1700_3_1':
    #     z_1700.append(2.318)
    # if tmp[0] == 'h1700_5_1': #photometric
    #     z_1700.append(2.3)
    # if tmp[0] == 'h1700_5_2':
    #     z_1700.append(2.313)
    # if tmp[0] == 'h1700_5_3': #photometric
    #     z_1700.append(2.3)
    # if tmp[0] == 'h1700_8_1':
    #     z_1700.append(2.303)
    # if tmp[0] == 'h1700_11': #photometric
    #     z_1700.append(2.3)
    # if tmp[0] == 'h1700_12':
    #     z_1700.append(2.72)
    # if tmp[0] == 'h1700_14':
    #     z_1700.append(2.306)
    # if tmp[0] == 'h1700_15_2':
    #     z_1700.append(2.30)
    # if tmp[0] == 'h1700_16':
    #     z_1700.append(1.575)
    # else:
    #     z_1700.append(0)
# if tmp[4] == '0':
#     ch4_1700.append(0)
# if tmp[5] == '0':
#     ch3_1700.append(0)
# if tmp[6] == '0':
#     ch2_1700.append(0)
# if tmp[7] == '0':
#     ch1_1700.append(0)
# if tmp[8] == '0':
#     K_1700.append(0)
# else:
    #adjusting for correct zero points from Damen+2011 - SIMPLE survey
    if tmp[4] == '0':
        ch4_1700.append(0)
    else:
        ch4_1700.append(float(tmp[4])+4.2165) #first number calibrated Damen+2011, second number recalibrated to Shapley+2005 paper
        # ch4_1700.append(float(tmp[4])+0.022840000000000416)#+4.2165) #first number calibrated Damen+2011, second number recalibrated to Shapley+2005 paper
    if tmp[5] == '0':
        ch3_1700.append(0)
        #13: 4.1712999999999987
        #5_1: 3.7771000000000008
        #5_2: 4.7298000000000009
        #17: 4.1448999999999998
        #18: 4.1451999999999991
    else:
        ch3_1700.append(float(tmp[5])+2.7884)
        # ch3_1700.append(float(tmp[5])+2.7884-2.307528571/2)
    if tmp[6] == '0':
        ch2_1700.append(0)
    else:
        ch2_1700.append(float(tmp[6])+4.1097)
        # ch2_1700.append(float(tmp[6])+0.18697857099999915)#+4.1097-2.307528571*2)
        #2_1: 0.038828570999999812
        #2_2: 0.25262857100000247
        #6_2: 0.1293285709999985
        #3_1: 0.364728570999997
        #3_2: 0.010428570999998499
        #3_3: 0.32592857099999861
    if tmp[7] == '0':
        ch1_1700.append(0)
    else:
        ch1_1700.append(float(tmp[7])+4.4000)
        # ch1_1700.append(float(tmp[7])-0.029799999999999827)#4.4000-4.28*1.035)
        #2_1: -0.21061666700000004
        #2_2: 0.062683333000002506
        #3_1: 0.11028333299999815
        #3_2: -0.1695166669999999
        #3_3: 0.1429833330000001
    if tmp[8] == '0':
        K_1700.append(0)
    else:
        K_1700.append(float(tmp[8]) + 1.85 - 1) #turn from K_vega to K_AB
cat_IR_1700.close()
ID_1700 = np.array(ID_1700)
ch4_1700 = np.array(ch4_1700)
ch3_1700 = np.array(ch3_1700)
ch2_1700 = np.array(ch2_1700)
ch1_1700 = np.array(ch1_1700)
K_1700 = np.array(K_1700)
z_1700 = np.array(z_1700)

cc_index_1700 = [np.where((ch2_1700>0) & (ch4_1700>0) & (ch1_1700>0))[0]]
z_1700_cluster_all = [np.where((z_1700>2.2) & (z_1700<2.4))[0]]
z_1700_cluster_IRAC = np.where((z_1700>2.2) & (z_1700<2.4) & ((ch2_1700>0) & (ch4_1700>0) & (ch1_1700>0)))[0]
z_1700_notcluster_all = [np.concatenate([np.where(z_1700<2.2)[0], np.where(z_1700>2.4)[0]])]
z_1700_notcluster_IRAC = np.concatenate([np.where((ch2_1700>0) & (ch4_1700>0) & (ch1_1700>0) & (z_1700<2.2) & (z_1700>0))[0],np.where((ch2_1700>0) & (ch4_1700>0) & (ch1_1700>0) & (z_1700>2.4))[0]])
z_1700_noz_IRAC = [np.where((ch2_1700>0) & (ch4_1700>0) & (ch1_1700>0) & (z_1700<2.2) & (z_1700<0))[0]]


index_1700_notcluster = np.array([np.where(ID_1700 == 'h1700_7_2')[0][0],np.where(ID_1700 == 'h1700_7_3')[0][0], np.where(ID_1700 == 'h1700_7_4')[0][0], np.where(ID_1700 == 'h1700_8_1')[0][0], np.where(ID_1700 == 'h1700_9')[0][0],np.where(ID_1700 == 'h1700_7_1')[0][0],np.where(ID_1700 == 'h1700_14')[0][0],np.where(ID_1700 == 'h1700_15_1')[0][0],np.where(ID_1700 == 'h1700_15_2')[0][0]])

# all:
# 1700.4, 1700.5\_2, 1700.7\_1,\_2,\_3,\_4, 1700.8, 1700.9, 1700.14, 1700.15\_1,\_2, 1700.16
# detections:
# 1700.4 (2.318), 1700.5\_2 (2.303), 1700.7\_1 (2.313), 1700.16 (1.575), 1700.17 (2.306),
# non detections
# 1700.7\_2,\_3,\_4, 1700.8, 1700.9, 1700.14, 1700.15\_1,\_2

z_15 = [np.where(z_ALESS[z_index_ALESS]<1.5)[0]]
z_1525 = [np.where((z_ALESS[z_index_ALESS]>1.5) & (z_ALESS[z_index_ALESS]<2.5))[0]]
z_25 = [np.where(z_ALESS[z_index_ALESS]>2.5)[0]]

colour_index = np.where((S_500>0)&(S_350>0)&(S_250>0))[0]
S_500_withdata = S_500[colour_index]
S_350_withdata = S_350[colour_index]
S_250_withdata = S_250[colour_index]
e_S_250_withdata = e_S_250[colour_index]
e_S_350_withdata = e_S_350[colour_index]
e_S_500_withdata = e_S_500[colour_index]
ID_withdata = ID[colour_index]


#CC figure
fig = pl.figure()#figsize=(10,4))
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font',size=15)
pl.rc('mathtext', default='regular')
ax = fig.add_subplot(1,1,1)
ax.set_xscale('linear')
ax.set_yscale('linear')


zspace = np.linspace(0.5,4,500)
S_all = [S_1,S_2,S_3]
wavelength_all = [wavelength_1,wavelength_2,wavelength_3]
label = ['ALMA composite', 'ALESS composite', 'Eyelash']
col = ['b','k','k']
for i in range(1):
    i=2
    y = []
    x = []
    # y_protoz = []
    # x_protoz = []
    # y_z3 = []
    # x_z3 = []
    # y_z15 = []
    # x_z15 = []
    #S = S_all[i]
    S = movingaverage(S_all[i],50)
    wavelength = wavelength_all[i]
    # if i == 2:
    #     wavelength/=(1+2.35)
    for z in zspace:
        S500 = S[np.where(wavelength*(1+z)>=500.)[0][0]]
        S350 = S[np.where(wavelength*(1+z)>=350.)[0][0]]
        S250 = S[np.where(wavelength*(1+z)>=250.)[0][0]]

        y.append(S500/S350)
        x.append(S350/S250)

    y = np.array(y)
    x = np.array(x)

    #y = movingaverage(y, 5)

#    pl.scatter(x,y,s=10,marker='.',c=zspace,cmap=pl.cm.get_cmap(cmap),vmin=vmin,vmax=vmax,zorder=10, label=label[i])
#    ax.plot(x,y,c=col[i],ls='-',label=label[i],zorder=0)
    ax.plot(x[np.where((zspace<1.5))[0]],y[np.where((zspace<1.5))[0]],c='b',ls='-',zorder=0)
    ax.plot(x[np.where((zspace>1.5) & (zspace<2.5))[0]],y[np.where((zspace>1.5) & (zspace<2.5))[0]],c='g',ls='-',zorder=0)
    ax.plot(x[np.where((zspace>2.5))[0]],y[np.where((zspace>2.5))[0]],c='r',ls='-',zorder=0)
    ax.plot(x[np.where((zspace>2.2) & (zspace<2.4))[0]],y[np.where((zspace>2.2) & (zspace<2.4))[0]],c='m',ls='-',zorder=0)


# draw box
# ax.hlines(y=0.85,linestyle='-',color='k',xmin=1.15,xmax=1.66,zorder=0)
# ax.hlines(y=1.53,linestyle='-',color='k',xmin=1.15,xmax=1.66,zorder=0)
# #pl.axhline(y=3.57,ls='--',c='k')
# ax.vlines(x=1.15,linestyle='-',color='k',ymin=0.85,ymax=1.53,zorder=0)
# ax.vlines(x=1.66,linestyle='-',color='k',ymin=0.85,ymax=1.53,zorder=0)
# #pl.axvline(x=1.82,ls='--',c='k')
# ax.text(0.68,1.35,r'$2.2<z<2.4$',color='k')
    if i==2:
        for j in np.linspace(0.5,4,8):
            pl.text(x[np.where(zspace<=j)[0][-1]],y[np.where(zspace<=j)[0][-1]],'z='+str(j),color='orange')
            #pl.axvline(x[np.where(zspace<=j)[0][-1]],ls='--',c='k')
# uncertanties propogated from confusion and intrument noise
# these are used in the X and Y values for CC plot
e_350_250 = 0.6*S_350_withdata/S_250_withdata*np.sqrt((e_S_350_withdata/S_350_withdata)**2 + (e_S_250_withdata/S_250_withdata)**2)
e_500_350 = 0.6*S_500_withdata/S_350_withdata*np.sqrt((e_S_350_withdata/S_350_withdata)**2 + (e_S_500_withdata/S_500_withdata)**2)
X_up = S_350_withdata / S_250_withdata + e_350_250
X_down = S_350_withdata / S_250_withdata - e_350_250
X = S_350_withdata/S_250_withdata
Y = S_500_withdata/S_350_withdata
print 'Photometric redshifts and comparison to known data'
det_lim=0.003
photo_z_1700 = []
e_photo_z_1700 = []
delta_z_1700 = []
print 'SCUBA-2 ID', '\t','photo-z','\t', 'e_photo-z', '\t', 'IRAC ID','\t', "'known-z'",'\t','delta-z (photo-known)'
for i in range(len(ID_withdata)):
    samez = 'N'
    for j in range(len(ID_1700)):
        if ID_1700[j].split('_')[1] == ID_withdata[i]:
            photo_z = round(zspace[np.where(abs(x-S_350_withdata[i]/S_250_withdata[i])<det_lim)[0][-1]],3)
            e_photo_z = zspace[np.where(abs(x-X_up[i])<det_lim)[0][-1]] - photo_z
            e_photo_z_1700.append(e_photo_z)
            photo_z_1700.append(photo_z)
            if z_1700[j] < 0:
                delta_z = z_1700[j]
            else:
                delta_z = round(z_1700[j],3) - photo_z
            delta_z_1700.append(delta_z)
            if round(e_photo_z,1) >= abs(round(delta_z,1)):
                samez = 'Y'
            print ID_withdata[i],'\t',round(photo_z,1),'\t', round(e_photo_z,1), '\t', ID_1700[j], '\t',round(z_1700[j],3), '\t',round(delta_z,1), '\t', samez
photo_z_1700 = np.array(photo_z_1700)
delta_z_1700 = np.array(delta_z_1700)
e_photo_z_1700 = np.array(e_photo_z_1700)

ax.scatter(X,Y, s=250, marker='*',edgecolors='k',facecolors='k',zorder=10)
ax.errorbar(X,Y,xerr=e_350_250,yerr=e_500_350,ls='none',c='k',alpha=0.5)
for i in range(len(ID_withdata)):
    pl.text(X[i]-0.05,Y[i]-0.09,ID_withdata[i])

# CC_1700_uplims_index = np.where((ch2_1700[index_1700_notcluster]>0) & (ch4_1700[index_1700_notcluster]>0))[0]
# CC_1700_uplims_x = 10**(-0.4*(ch2_1700[index_1700_notcluster]-ch1_1700[index_1700_notcluster]))[CC_1700_uplims_index]
# CC_1700_uplims_y = 10**(-0.4*(ch4_1700[index_1700_notcluster]-ch2_1700[index_1700_notcluster]))[CC_1700_uplims_index]
# ch1_1700_lowlim_mag = np.min(ch1_1700[np.where(ch1_1700>20.1)[0]]) #np.mean(ch1_1700[np.where(ch1_1700>0)[0]])
# ch1_1700_lowlim_array = np.array(len(CC_1700_uplims_index)*[ch1_1700_lowlim_mag])
# ch1_1700_min = 10**(-0.4*(ch2_1700[index_1700_notcluster][CC_1700_uplims_index] - ch1_1700_lowlim_array))
# ax.scatter(len(CC_1700_uplims_y)*[0.6],CC_1700_uplims_y,s=250, marker='*',edgecolors='k',facecolors='k',zorder=10)
# pl.errorbar(len(CC_1700_uplims_y)*[0.6],CC_1700_uplims_y, xerr=[len(CC_1700_uplims_y)*[0],len(CC_1700_uplims_y)*[0.08]], xlolims=True, ls='None', c='k')



ax.set_ylabel(r'$S_{500}/S_{350}$')
ax.set_xlabel(r'$S_{350}/S_{250}$')

ax.set_xlim(0.3,2)
ax.set_ylim(0.3,1.5)
# ax.legend(fontsize=12,loc=2,ncol=1,frameon=True,numpoints=1,scatterpoints=1)

# texty=0.32
# ax.text(1.6,texty,'HS1700',color='k',size=16)


#pl.xticks([0.6,1,2])
# ax.set_yticks([0.5,1,2,3,4])
# ax.set_yticklabels(['0.5','1','2','3','4'])
#ax.set_xticklabels(['0.6','1','2'])
#pl.rcParams['axes.formatter.min_exponent'] = 100

#pl.savefig('../../Figures/Colour/CM-CC_1700.pdf', bbox_inches='tight')
#pl.show()
#pl.close()
#pl.rcParams['axes.formatter.min_exponent'] = 100


pl.savefig('../../Figures/Colour/1700_SPIRE_colour.png')#, bbox_inches='tight')
#pl.show()
pl.close()



# print 'IDs with no '
# for i in ID[np.where((S_500<0)&(S_350<0)&(S_250<0))[0]]:
