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


# z_ALESS = np.zeros(len(ID_ALESS))
# cat_ALESS_z = open('ALESS_SMGs_Simpson_2014_z.txt','r')
# ID_z_ALESS = []
# # e_z_up_ALESS = []
# # e_z_low_ALESS = []
# for line in cat_ALESS_z.readlines():
#     tmp = line.split()
#     ID_z_ALESS.append(tmp[0]+'_'+tmp[1])
#     z_ALESS[np.where(tmp[0]+'_'+tmp[1] == ID_ALESS)[0]]+=float(tmp[8])
# cat_ALESS_z.close()
# ID_z_ALESS = np.array(ID_z_ALESS)
# z_ALESS = np.array(z_ALESS)
# z_index_ALESS = np.where(z_ALESS>0)[0]
# # e_z_up_ALESS = np.array(e_z_up_ALESS)
# # e_z_low_ALESS = np.array(e_z_low_ALESS)


#open file
cat_IR_1549 = open('../../h1549/mutli_wavelength.gaia','r')
#import 850um data
ID_1549 = []
ch4_1549 = []
ch3_1549 = []
ch2_1549 = []
ch1_1549 = []
K_1549 = []
z_1549 = []
for line in cat_IR_1549.readlines()[2:]:
    tmp = line.split()
    ID_1549.append(tmp[0])
    z_1549.append(float(tmp[-1]))
    if float(tmp[4]) == 0:
        ch4_1549.append(0)
    else:
        ch4_1549.append(float(tmp[4])+4.2165)
    if float(tmp[5]) == 0:
        ch3_1549.append(0)
    else:
        ch3_1549.append(float(tmp[5])+2.7884)
    if float(tmp[6]) == 0:
        ch2_1549.append(0)
    else:
        ch2_1549.append(float(tmp[6])+4.1097)
    if float(tmp[7]) == 0:
        ch1_1549.append(0)
    else:
        ch1_1549.append(float(tmp[7])+4.4000)
    if float(tmp[8]) == 0:
        K_1549.append(0)
    else:
        K_1549.append(float(tmp[8]) + 1.9) #turn from K_vega to K_AB
cat_IR_1549.close()
ID_1549 = np.array(ID_1549)
ch4_1549 = np.array(ch4_1549)
ch3_1549 = np.array(ch3_1549)
ch2_1549 = np.array(ch2_1549)
ch1_1549 = np.array(ch1_1549)
K_1549 = np.array(K_1549)
z_1549 = np.array(z_1549)
cc_index_1549 = [np.where((ch2_1549>0) & (ch4_1549>0))[0]]

z_1549_cluster_all = [np.where((z_1549>2.75) & (z_1549<2.95))[0]]
z_1549_cluster_IRAC = np.where((z_1549>2.75) & (z_1549<2.95) & ((ch2_1549>0) & (ch4_1549>0)))[0]
z_1549_notcluster_all = [np.concatenate([np.where(z_1549<2.75)[0], np.where(z_1549>2.95)[0]])]
z_1549_notcluster_IRAC = np.concatenate([np.where((ch2_1549>0) & (ch4_1549>0) & (z_1549<2.75) & (z_1549>0))[0],np.where((ch2_1549>0) & (ch4_1549>0) & (z_1549>2.95))[0]])
z_1549_noz_IRAC = [np.where((ch2_1549>0) & (ch4_1549>0) & (z_1549<2.75) & (z_1549<0))[0]]

index_1549_notcluster = np.array([np.where(ID_1549 == 'h1549_10')[0][0]])

S_ch4_1549 = 10**(ch4_1549/2.5)
S_ch3_1549 = 10**(ch3_1549/2.5)
S_ch2_1549 = 10**(ch2_1549/2.5)
S_ch1_1549 = 10**(ch1_1549/2.5)
S_K_1549 = 10**(K_1549/2.5)

#open file
cat_IR_1700 = open('../../h1700/mutli_wavelength.gaia','r')
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


#CC figure
fig = pl.figure(figsize=(10,4))
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font',size=15)
pl.rc('mathtext', default='regular')
ax = fig.add_subplot(1,2,1)
ax.set_xscale('log')
ax.set_yscale('log')




z_15 = [np.where(z_ALESS[z_index_ALESS]<1.5)[0]]
z_1525 = [np.where((z_ALESS[z_index_ALESS]>1.5) & (z_ALESS[z_index_ALESS]<2.5))[0]]
z_25 = [np.where(z_ALESS[z_index_ALESS]>2.5)[0]]

### Derivation on errorbars
# d(10^(ch2-ch1)) = (d(10^-0.4(ch2-ch1))/dch2 * e_ch2)^2 + (d(10^-0.4(ch2-ch1))/dch1 * e_ch1)^2
#  = (-0.4*ln(10)*10^-0.4(ch2-ch1)*e_ch2)^2 + (-0.4*ln(10)*10^-0.4(ch2-ch1)*e_ch1)^2
#  = (0.4^2*ln(10)^2)^2*10^-0.4(ch2-ch1) * (e_ch2^2 + e_ch1^2)^2
#  d(10^(ch2-ch1)) = (0.4*ln(10))^2*10^-0.4(ch2-ch1) * sqrt(e_ch2^2 + e_ch1^2)


# pl.scatter(10**(-0.4*(ch2_ALESS-ch1_ALESS)),10**(-0.4*(ch4_ALESS-ch2_ALESS)),s=20,facecolors='m',edgecolors='m', label='ALESS SMGs')
ax.scatter(10**(-0.4*(ch2_ALESS[z_index_ALESS][z_15]-ch1_ALESS[z_index_ALESS][z_15])),10**(-0.4*(ch4_ALESS[z_index_ALESS][z_15]-ch2_ALESS[z_index_ALESS][z_15])),s=20,facecolors='b',edgecolors='b',alpha=0.5)# label=r'$z<1.5$',alpha=0.5)
pl.errorbar(10**(-0.4*(ch2_ALESS-ch1_ALESS))[z_index_ALESS][z_15], 10**(-0.4*(ch4_ALESS-ch2_ALESS))[z_index_ALESS][z_15], xerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch2_ALESS-ch1_ALESS)) * np.sqrt(e_ch2_ALESS**2 + e_ch1_ALESS**2))[z_index_ALESS][z_15], yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_ALESS-ch2_ALESS)) * np.sqrt(e_ch4_ALESS**2 + e_ch2_ALESS**2))[z_index_ALESS][z_15], ls='None',c='b',alpha=0.5)
ax.scatter(10**(-0.4*(ch2_ALESS[z_index_ALESS][z_1525]-ch1_ALESS[z_index_ALESS][z_1525])),10**(-0.4*(ch4_ALESS[z_index_ALESS][z_1525]-ch2_ALESS[z_index_ALESS][z_1525])),s=20,facecolors='g',edgecolors='g',alpha=0.5)#label=r'$1.5<z<2.5$',alpha=0.5)
pl.errorbar(10**(-0.4*(ch2_ALESS-ch1_ALESS))[z_index_ALESS][z_1525], 10**(-0.4*(ch4_ALESS-ch2_ALESS))[z_index_ALESS][z_1525], xerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch2_ALESS-ch1_ALESS)) * np.sqrt(e_ch2_ALESS**2 + e_ch1_ALESS**2))[z_index_ALESS][z_1525], yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_ALESS-ch2_ALESS)) * np.sqrt(e_ch4_ALESS**2 + e_ch2_ALESS**2))[z_index_ALESS][z_1525], ls='None',c='g',alpha=0.5)
ax.scatter(10**(-0.4*(ch2_ALESS[z_index_ALESS][z_25]-ch1_ALESS[z_index_ALESS][z_25])),10**(-0.4*(ch4_ALESS[z_index_ALESS][z_25]-ch2_ALESS[z_index_ALESS][z_25])),s=20,facecolors='r',edgecolors='r',alpha=0.5)#label=r'$z>2.5$',alpha=0.5)
pl.errorbar(10**(-0.4*(ch2_ALESS-ch1_ALESS))[z_index_ALESS][z_25], 10**(-0.4*(ch4_ALESS-ch2_ALESS))[z_index_ALESS][z_25], xerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch2_ALESS-ch1_ALESS)) * np.sqrt(e_ch2_ALESS**2 + e_ch1_ALESS**2))[z_index_ALESS][z_25], yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_ALESS-ch2_ALESS)) * np.sqrt(e_ch4_ALESS**2 + e_ch2_ALESS**2))[z_index_ALESS][z_25], ls='None',c='r',alpha=0.5)#pl.scatter(10**(-0.4*(ch2_1700[cc_index_1700]-ch1_1700[cc_index_1700])),10**(-0.4*(ch4_1700[cc_index_1700]-ch2_1700[cc_index_1700])),s=80, marker='*', facecolors='k', edgecolors='k', label='HS1700')
vmin=2.2
vmax=2.4
cmap='Spectral_r'#'RdYlBu_r'
# cax=pl.scatter(np.concatenate([10**(-0.4*(ch2_ALESS[z_index_ALESS]-ch1_ALESS[z_index_ALESS])),10**(-0.4*(ch2_1700[cc_index_1700]-ch1_1700[cc_index_1700]))]),np.concatenate([10**(-0.4*(ch4_ALESS[z_index_ALESS]-ch2_ALESS[z_index_ALESS])),10**(-0.4*(ch4_1700[cc_index_1700]-ch2_1700[cc_index_1700]))]),c=np.concatenate([z_ALESS[z_index_ALESS],z_1700[cc_index_1700]]), cmap=pl.cm.get_cmap(cmap),marker='o',s=1,vmin=vmin,vmax=vmax)
# pl.scatter(10**(-0.4*(ch2_ALESS[z_index_ALESS]-ch1_ALESS[z_index_ALESS])),10**(-0.4*(ch4_ALESS[z_index_ALESS]-ch2_ALESS[z_index_ALESS])),s=10,c=z_ALESS[z_index_ALESS],cmap=pl.cm.get_cmap(cmap),label='ALESS SMGs',vmin=vmin,vmax=vmax,zorder=10)

# test = []
# for i in range(len(e_ch2_ALESS[np.where(e_ch2_ALESS>0)[0]])):
#     test.append(abs(r.gauss(-0.05,0.05)))
# pl.hist(e_ch2_ALESS[np.where(e_ch2_ALESS>0)[0]])
# pl.hist(np.array(test))
# print test
# pl.show()

# test = []
# for i in range(len(e_ch2_ALESS[np.where(e_ch2_ALESS>0)[0]])):
#     test.append(abs(r.gauss(-0.05,0.05))+0.01)
# print test
# test = np.array(test)
# pl.hist(e_ch2_ALESS[np.where(e_ch2_ALESS>0)[0]])
# pl.hist(test)
# pl.show()

rand_err_1700_x = np.array([0.06052617381624862, 0.047466240358435045, 0.14565001024596663, 0.14187395541007097, 0.05256994025141485, 0.03679948193947722, 0.14690369885541193, 0.01061445536201775, 0.08112709942700791, 0.04742829628652853, 0.0777275083122609, 0.08047671998654328, 0.019200400193467727, 0.01612633885564229, 0.017309180014689667, 0.09220814285506398, 0.019249078079014527, 0.07106154031616317, 0.13252483084364577, 0.026752395270954916, 0.06244887886097028, 0.05960571552210332, 0.09559089653749837, 0.06339532151753972, 0.025408801016298242, 0.0718774319000873, 0.034778115243207085, 0.011812009914087615, 0.09034582973841555, 0.06021406360623804, 0.08604826304533274, 0.027800113501988792, 0.07593829431856411, 0.10743176813186464, 0.10689848027090026, 0.028914858697167943, 0.024008783393633534, 0.12433487380451562, 0.038623316183355356, 0.02646393025010084, 0.018817318614396533, 0.12355603564145128, 0.16097813146338213, 0.13911407643517582, 0.046830191931906226, 0.08542069271778381, 0.0698406429871322, 0.0969122873410515, 0.05051634754228539, 0.10800286132775379, 0.06232522662343778, 0.012655809487499855, 0.024696016266573008, 0.037628330560278755, 0.04422312513747239, 0.08760509997925409, 0.023937364855624678, 0.0209451752898894, 0.1303999367646611, 0.04228857711476672, 0.035213616101014276, 0.010072767680294364, 0.07636902097954104, 0.11292916970260665, 0.03770135672593327, 0.07395490726643161, 0.12738361531017506, 0.019166752531648838, 0.070447368598301, 0.010611225493807197, 0.05065194270693315, 0.10863780116447708, 0.04914173538407655])
rand_err_1700_y = np.array([0.12329141019672503, 0.06400492577085044, 0.038439878064725534, 0.10593660736700261, 0.016113403914556658, 0.030969704702063137, 0.03399615758417491, 0.015889877023365653, 0.07459957400048446, 0.08282953458686156, 0.09136202010900234, 0.06910419836266408, 0.0726769666057652, 0.12739290505849446, 0.06477204679783445, 0.06804701004423683, 0.10542434468143379, 0.0910563665412577, 0.06282714687674561, 0.0810642323935348, 0.192171198226271, 0.02520419105111356, 0.08313035538337653, 0.051808261845834765, 0.05135045451862135, 0.07809389951458857, 0.17991509491275065, 0.03253287544415223, 0.05105180017681247, 0.10253784602938103, 0.0151635481963715, 0.14156061026253441, 0.0876261460486051, 0.08207208971820437, 0.02274872878951454, 0.1896751458181895, 0.08899791668452532, 0.12620996457830297, 0.05349413524789995, 0.09131156756381098, 0.025403810435441897, 0.07640333061555225, 0.017131159142505716, 0.06088569359148386, 0.023604752470009337, 0.08947592565523471, 0.03355436540857107, 0.08282747325822372, 0.06643843390356936, 0.12415973022735148, 0.03248980415233448, 0.05641815765230703, 0.18178596052189827, 0.09594513236667145, 0.04272677742694453, 0.08689701654745781, 0.12111391162204785, 0.07611241443723975, 0.024208055521125288, 0.19044136292597336, 0.08798846003590889, 0.11304062778041786, 0.049393624035697026, 0.05929985663652515, 0.05768486012153358, 0.0647972031074591, 0.09870152828627661, 0.09611940785046114, 0.1463328485599144, 0.023617355984723372, 0.14399451406243755, 0.11832411914365382, 0.12307310222179801])
rand_err_1549_x = np.array([0.06855296722716439, 0.025872505979617287, 0.1282122153849773, 0.14698830580201377, 0.03194817937945236, 0.03491102185699824, 0.09951523694507731, 0.10386080626906553, 0.03925682457209818, 0.16647259365677464, 0.1348778443372043, 0.13123107007531423, 0.12086403338446516, 0.04313863263379337, 0.07545260267954589, 0.06401104542235564, 0.08435603141931385, 0.04173808868237791, 0.04051935977918145, 0.061325011773086795, 0.03852636280947324, 0.0864641319432526, 0.010786096078072722, 0.07484853819384796, 0.09062788836491369, 0.14576018882014796, 0.023471325335310912, 0.12478750868068729, 0.07954862068473532, 0.06456199085140946, 0.03640977375934027, 0.027521941732073528, 0.08027404828087845, 0.1099229027590334, 0.012848807010481652, 0.09050178307053502, 0.03194171915404266, 0.0724767097151418, 0.03319797722020731, 0.055337156619837055, 0.1670406174339644, 0.06575080756846165, 0.0624972040035378, 0.07894322573401953, 0.12095394070587132, 0.03072695721849896, 0.025923526229455686, 0.10414335029309289, 0.020776301381936327, 0.06488538008893044, 0.030812549749462834, 0.06169009858174754, 0.04627249060026108, 0.06490596034080516, 0.15206969192670977, 0.01673134305361656, 0.03948237739099979, 0.07991449408351667, 0.03822544187803662, 0.012251791792937316, 0.015346385337709911, 0.06077323859871323, 0.07698415983717231, 0.1683670636577168, 0.07338103883982665, 0.03825462681070111, 0.04536316174104864, 0.1119331774301585, 0.047748384928661236, 0.03088792656145651, 0.07488814864488978, 0.042237003353944404, 0.10868237405361812])
rand_err_1549_y = np.array([0.08843759926008775, 0.15162062391452957, 0.02097155776737692, 0.1028656681953122, 0.06504549869294356, 0.12333175236909792, 0.02240125468775083, 0.07653090486917727, 0.06819705490666825, 0.026402752017484953, 0.024995876966813463, 0.1032812078682189, 0.014045466667157244, 0.08824264393041863, 0.0700954947916613, 0.026008705800483423, 0.10647457917193312, 0.06635654432734184, 0.0625181791977513, 0.10847737113844645, 0.02528761554937735, 0.07197958555784068, 0.12480326354023707, 0.13081193798136964, 0.12808519474265195, 0.04545013364688007, 0.07256582384871572, 0.06492547220998153, 0.06727023066365902, 0.19815741651447843, 0.08476053727967416, 0.07306058309982814, 0.04099702059431125, 0.09372178727995917, 0.04011401610858313, 0.10964977193903973, 0.1545813447207956, 0.0238198229513285, 0.061877123135018, 0.07120840884812094, 0.02746312301856841, 0.06972393257261784, 0.015242586153529946, 0.017024187965951176, 0.0853185015339937, 0.15506857820521786, 0.09647677297547884, 0.025329528737255823, 0.07257176582041168, 0.014660420088378643, 0.011885015298354817, 0.062035896702984865, 0.10193929577987296, 0.11100036586290678, 0.12574437393973073, 0.121624863844178, 0.016292784283045973, 0.03979049364541579, 0.037834507098497984, 0.13628020753859, 0.11488806244148754, 0.05398110504624414, 0.022663288702858657, 0.036742459140693524, 0.0900361414012542, 0.05751192998076496, 0.14221751659612625, 0.12384554122739404, 0.011956178042068643, 0.14000937497424462, 0.07547897559713927, 0.13949056489353567, 0.022618338641255835])


ax.scatter(10**(-0.4*(ch2_1700[z_1700_noz_IRAC]-ch1_1700[z_1700_noz_IRAC])),10**(-0.4*(ch4_1700[z_1700_noz_IRAC]-ch2_1700[z_1700_noz_IRAC])),s=250, marker='*',edgecolors='k',facecolors='w',zorder=10)
pl.errorbar(10**(-0.4*(ch2_1700[z_1700_noz_IRAC]-ch1_1700[z_1700_noz_IRAC])),10**(-0.4*(ch4_1700[z_1700_noz_IRAC]-ch2_1700[z_1700_noz_IRAC])), xerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch2_1700[z_1700_noz_IRAC]-ch1_1700[z_1700_noz_IRAC]))) * np.sqrt(rand_err_1700_x[z_1700_noz_IRAC]**2 + rand_err_1700_x[z_1700_noz_IRAC]**2), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1700[z_1700_noz_IRAC]-ch2_1700[z_1700_noz_IRAC]))) * np.sqrt(rand_err_1700_y[z_1700_noz_IRAC]**2 + rand_err_1700_y[z_1700_noz_IRAC]**2), ls='None',c='k',alpha=0.5)

ax.scatter(10**(-0.4*(ch2_1700[z_1700_cluster_IRAC]-ch1_1700[z_1700_cluster_IRAC])),10**(-0.4*(ch4_1700[z_1700_cluster_IRAC]-ch2_1700[z_1700_cluster_IRAC])),s=250, marker='*',facecolors='g',edgecolors='k',zorder=10)
pl.errorbar(10**(-0.4*(ch2_1700[z_1700_cluster_IRAC]-ch1_1700[z_1700_cluster_IRAC])),10**(-0.4*(ch4_1700[z_1700_cluster_IRAC]-ch2_1700[z_1700_cluster_IRAC])), xerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch2_1700[z_1700_cluster_IRAC]-ch1_1700[z_1700_cluster_IRAC]))) * np.sqrt(rand_err_1700_x[z_1700_cluster_IRAC]**2 + rand_err_1700_x[z_1700_cluster_IRAC]**2), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1700[z_1700_cluster_IRAC]-ch2_1700[z_1700_cluster_IRAC]))) * np.sqrt(rand_err_1700_y[z_1700_cluster_IRAC]**2 + rand_err_1700_y[z_1700_cluster_IRAC]**2), ls='None',c='g',alpha=1)

''' This is 1700.16 with z=1.575 and doesn't fall within the 1.5<z<2.5 redshift bin, so it should be green'''
ax.scatter(10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC]))[0],10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC]))[0],s=250, marker='*', facecolors='g',edgecolors='k',zorder=10)
pl.errorbar(10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC])),10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC])), xerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC]))) * np.sqrt(rand_err_1700_x[z_1700_notcluster_IRAC]**2 + rand_err_1700_x[z_1700_notcluster_IRAC]**2), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC]))) * np.sqrt(rand_err_1700_y[z_1700_notcluster_IRAC]**2 + rand_err_1700_y[z_1700_notcluster_IRAC]**2), ls='None',c='g',alpha=1)

ax.scatter(10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC]))[1],10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC]))[1],s=250, marker='*', facecolors='red',edgecolors='k',zorder=10)
pl.errorbar(10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC]))[1],10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC]))[1], xerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC][1]-ch1_1700[z_1700_notcluster_IRAC][1]))) * np.sqrt(rand_err_1700_x[z_1700_notcluster_IRAC][1]**2 + rand_err_1700_x[z_1700_notcluster_IRAC][1]**2), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC][1]-ch2_1700[z_1700_notcluster_IRAC][1]))) * np.sqrt(rand_err_1700_y[z_1700_notcluster_IRAC][1]**2 + rand_err_1700_y[z_1700_notcluster_IRAC][1]**2), ls='None',c='r',alpha=1)

ax.scatter(10**(-0.4*(ch2_1700[index_1700_notcluster]-ch1_1700[index_1700_notcluster])),10**(-0.4*(ch4_1700[index_1700_notcluster]-ch2_1700[index_1700_notcluster])),s=250, marker='*',edgecolors='k',facecolors='k',zorder=10)
pl.errorbar(10**(-0.4*(ch2_1700[index_1700_notcluster]-ch1_1700[index_1700_notcluster])),10**(-0.4*(ch4_1700[index_1700_notcluster]-ch2_1700[index_1700_notcluster])), xerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch2_1700[index_1700_notcluster]-ch1_1700[index_1700_notcluster]))) * np.sqrt(rand_err_1700_x[index_1700_notcluster]**2 + rand_err_1700_x[index_1700_notcluster]**2), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1700[index_1700_notcluster]-ch2_1700[index_1700_notcluster]))) * np.sqrt(rand_err_1700_y[index_1700_notcluster]**2 + rand_err_1700_y[index_1700_notcluster]**2), ls='None',c='k',alpha=1)


CC_1700_uplims_index = np.where((ch2_1700[index_1700_notcluster]>0) & (ch4_1700[index_1700_notcluster]>0))[0]
CC_1700_uplims_x = 10**(-0.4*(ch2_1700[index_1700_notcluster]-ch1_1700[index_1700_notcluster]))[CC_1700_uplims_index]
CC_1700_uplims_y = 10**(-0.4*(ch4_1700[index_1700_notcluster]-ch2_1700[index_1700_notcluster]))[CC_1700_uplims_index]
ch1_1700_lowlim_mag = np.min(ch1_1700[np.where(ch1_1700>20.1)[0]]) #np.mean(ch1_1700[np.where(ch1_1700>0)[0]])
ch1_1700_lowlim_array = np.array(len(CC_1700_uplims_index)*[ch1_1700_lowlim_mag])
ch1_1700_min = 10**(-0.4*(ch2_1700[index_1700_notcluster][CC_1700_uplims_index] - ch1_1700_lowlim_array))
ax.scatter(len(CC_1700_uplims_y)*[0.6],CC_1700_uplims_y,s=250, marker='*',edgecolors='k',facecolors='k',zorder=10)
pl.errorbar(len(CC_1700_uplims_y)*[0.6],CC_1700_uplims_y, xerr=[len(CC_1700_uplims_y)*[0],len(CC_1700_uplims_y)*[0.08]], xlolims=True, ls='None', c='k')

# pl.scatter(10**(-0.4*(ch2_1700[z_1700_cluster_IRAC]-ch1_1700[z_1700_cluster_IRAC])),10**(-0.4*(ch4_1700[z_1700_cluster_IRAC]-ch2_1700[z_1700_cluster_IRAC])),s=200, marker='*', c=z_1700[z_1700_cluster_IRAC], cmap=pl.cm.get_cmap(cmap),vmin=vmin,vmax=vmax,edgecolors='k',zorder=10)
# pl.scatter(10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC])),10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC])),s=200, marker='*', c=z_1700[z_1700_notcluster_IRAC], cmap=pl.cm.get_cmap(cmap), label='HS1700',vmin=vmin,vmax=vmax,edgecolors='k',zorder=10)


#pl.scatter(10**(-0.4*(ch2_1700[z_1700_notcluster]-ch1_1700[z_1700_notcluster])),10**(-0.4*(ch4_1700[z_1700_notcluster]-ch2_1700[z_1700_notcluster])),s=200, marker='*',edgecolors='k',facecolors='k')
# cax=pl.scatter(np.concatenate([10**(-0.4*(ch2_ALESS[z_index_ALESS_1700]-ch1_ALESS[z_index_ALESS_1700])),10**(-0.4*(ch2_1700[z_1700_cluster]-ch1_1700[z_1700_cluster]))]),np.concatenate([10**(-0.4*(ch4_ALESS[z_index_ALESS_1700]-ch2_ALESS[z_index_ALESS_1700])),10**(-0.4*(ch4_1700[z_1700_cluster]-ch2_1700[z_1700_cluster]))]),c=np.concatenate([z_ALESS[z_index_ALESS_1700],z_1700[z_1700_cluster]]), cmap=pl.cm.get_cmap(cmap),marker='o',s=1,vmin=vmin,vmax=vmax)
# pl.scatter(10**(-0.4*(ch2_ALESS[z_index_ALESS_1700]-ch1_ALESS[z_index_ALESS_1700])),10**(-0.4*(ch4_ALESS[z_index_ALESS_1700]-ch2_ALESS[z_index_ALESS_1700])),s=10,c=z_ALESS[z_index_ALESS_1700],cmap=pl.cm.get_cmap(cmap),label='ALESS SMGs',vmin=vmin,vmax=vmax)
# pl.scatter(10**(-0.4*(ch2_1700[z_1700_cluster]-ch1_1700[z_1700_cluster])),10**(-0.4*(ch4_1700[z_1700_cluster]-ch2_1700[z_1700_cluster])),s=200, marker='*', c=z_1700[z_1700_cluster], cmap=pl.cm.get_cmap(cmap), label='HS1700',vmin=vmin,vmax=vmax,edgecolors='k')

# template = np.loadtxt('SED/ALESS_compositeSED_JMS2015_LBOL.dat',unpack=True)
# wavelength = template[0]*1e-4
# S = template[1]

zspace = np.linspace(0.6,6,500)
S_all = [S_1,S_2,S_3]
wavelength_all = [wavelength_1,wavelength_2,wavelength_3]
#label = ['ALMA composite', 'ALESS composite', 'Eyelash']
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
        Sch4 = S[np.where(wavelength*(1+z)>=8.0)[0][0]]
        Sch3 = S[np.where(wavelength*(1+z)>=5.8)[0][0]]
        Sch2 = S[np.where(wavelength*(1+z)>=4.5)[0][0]]
        Sch1 = S[np.where(wavelength*(1+z)>=3.6)[0][0]]

        y.append(Sch4/Sch2)
        x.append(Sch2/Sch1)

    y = np.array(y)
    x = np.array(x)

    #y = movingaverage(y, 5)

#    pl.scatter(x,y,s=10,marker='.',c=zspace,cmap=pl.cm.get_cmap(cmap),vmin=vmin,vmax=vmax,zorder=10, label=label[i])
#    ax.plot(x,y,c=col[i],ls='-',label=label[i],zorder=0)
    ax.plot(x[np.where((zspace<1.5))[0]],y[np.where((zspace<1.5))[0]],c='b',ls='-',zorder=0)
    ax.plot(x[np.where((zspace>1.5) & (zspace<2.5))[0]],y[np.where((zspace>1.5) & (zspace<2.5))[0]],c='g',ls='-',zorder=0)
    ax.plot(x[np.where((zspace>2.5))[0]],y[np.where((zspace>2.5))[0]],c='r',ls='-',zorder=0)

ax.hlines(y=0.85,linestyle='-',color='k',xmin=1.15,xmax=1.66,zorder=0)
ax.hlines(y=1.53,linestyle='-',color='k',xmin=1.15,xmax=1.66,zorder=0)
#pl.axhline(y=3.57,ls='--',c='k')
ax.vlines(x=1.15,linestyle='-',color='k',ymin=0.85,ymax=1.53,zorder=0)
ax.vlines(x=1.66,linestyle='-',color='k',ymin=0.85,ymax=1.53,zorder=0)
#pl.axvline(x=1.82,ls='--',c='k')
ax.text(0.68,1.35,r'$2.2<z<2.4$',color='k')

'''
#cbar = fig.colorbar(cax,ticks=range(vmin,vmax+1),label='redshift')#,cmap='Spectral_r',vmin=0,vmax=6, ticks=[0,1,2,3,4,5,6])
cbar = fig.colorbar(cax,ticks=[2.2,2.3,2.4],label='redshift')#,cmap='Spectral_r',vmin=0,vmax=6, ticks=[0,1,2,3,4,5,6])
# cbar_index = []
# for i in range(vmin,vmax+1):
#     cbar_index.append(str(i))
# cbar.ax.set_yticklabels(cbar_index)  # vertically oriented colorbar
cbar.ax.set_yticklabels(['2.2','2.3','2.4'])  # vertically oriented colorbar
# pl.xlim(0.5,2.1)
# pl.ylim(0.3,4.5)
'''

ax.set_ylabel(r'$S_{8.0}/S_{4.5}$')
ax.set_xlabel(r'$S_{4.5}/S_{3.6}$')

ax.set_xlim(0.5,2.3)
ax.set_ylim(0.3,4.5)
# ax.legend(fontsize=12,loc=2,ncol=1,frameon=True,numpoints=1,scatterpoints=1)

texty=0.32
ax.text(0.515,texty,'HS1700',color='k',size=16)


#pl.xticks([0.6,1,2])
ax.set_yticks([0.5,1,2,3,4])
ax.set_yticklabels(['0.5','1','2','3','4'])
#ax.set_xticklabels(['0.6','1','2'])
#pl.rcParams['axes.formatter.min_exponent'] = 100

#pl.savefig('../../Figures/Colour/CM-CC_1700.pdf', bbox_inches='tight')
#pl.show()
#pl.close()




'''CM PLOT'''
#CM figure
ax2 = fig.add_subplot(1,2,2)
ax2.set_xscale('log')
ax2.set_yscale('log')

fudgefactor = 30.0

# pl.scatter(10**((23.9-ch2_1549[cc_index_1549])/2.5)*fudgefactor,10**((23.9-ch4_1549[cc_index_1549])/2.5) / (10**((23.9-ch2_1549[cc_index_1549])/2.5)),s=250,marker='*',facecolors='k',edgecolors='k',alpha=0.5)
# pl.scatter(10**((23.9-ch2_ALESS[z_index_ALESS])/2.5),10**((23.9-ch4_ALESS[z_index_ALESS])/2.5) / (10**((23.9-ch2_ALESS[z_index_ALESS])/2.5)),s=10,marker='o',facecolors='k',edgecolors='k',alpha=0.5)

# pl.scatter(10**(-0.4*(ch2_ALESS-ch1_ALESS)),10**(-0.4*(ch4_ALESS-ch2_ALESS)),s=20,facecolors='m',edgecolors='m', label='ALESS SMGs')
ax2.scatter(10**((23.9-ch2_ALESS[z_index_ALESS][z_15])/2.5),10**((23.9-ch4_ALESS[z_index_ALESS][z_15])/2.5) / (10**((23.9-ch2_ALESS[z_index_ALESS][z_15])/2.5)),s=20,facecolors='b',edgecolors='b', alpha=0.5)#label=r'$z<1.5$',alpha=0.5)
pl.errorbar(10**((23.9-ch2_ALESS)/2.5)[z_index_ALESS][z_15], 10**(-0.4*(ch4_ALESS-ch2_ALESS))[z_index_ALESS][z_15], xerr=((0.4*np.log(10))**2 * 10**((23.9-ch2_ALESS)/2.5) * np.sqrt(e_ch2_ALESS**2))[z_index_ALESS][z_15], yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_ALESS-ch2_ALESS)) * np.sqrt(e_ch4_ALESS**2 + e_ch2_ALESS**2))[z_index_ALESS][z_15], ls='None',c='b',alpha=0.5)#pl.scatter(10**(-0.4*(ch2_1700[cc_index_1700]-ch1_1700[cc_index_1700])),10**(-0.4*(ch4_1700[cc_index_1700]-ch2_1700[cc_index_1700])),s=80, marker='*', facecolors='k', edgecolors='k', label='HS1700')
ax2.scatter(10**((23.9-ch2_ALESS[z_index_ALESS][z_1525])/2.5),10**((23.9-ch4_ALESS[z_index_ALESS][z_1525])/2.5) / (10**((23.9-ch2_ALESS[z_index_ALESS][z_1525])/2.5)),s=20,facecolors='g',edgecolors='g', alpha=0.5)#label=r'$1.5<z<2.5$',alpha=0.5)
pl.errorbar(10**((23.9-ch2_ALESS)/2.5)[z_index_ALESS][z_1525], 10**(-0.4*(ch4_ALESS-ch2_ALESS))[z_index_ALESS][z_1525], xerr=((0.4*np.log(10))**2 * 10**((23.9-ch2_ALESS)/2.5) * np.sqrt(e_ch2_ALESS**2))[z_index_ALESS][z_1525], yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_ALESS-ch2_ALESS)) * np.sqrt(e_ch4_ALESS**2 + e_ch2_ALESS**2))[z_index_ALESS][z_1525], ls='None',c='g',alpha=0.5)#pl.scatter(10**(-0.4*(ch2_1700[cc_index_1700]-ch1_1700[cc_index_1700])),10**(-0.4*(ch4_1700[cc_index_1700]-ch2_1700[cc_index_1700])),s=80, marker='*', facecolors='k', edgecolors='k', label='HS1700')
ax2.scatter(10**((23.9-ch2_ALESS[z_index_ALESS][z_25])/2.5),10**((23.9-ch4_ALESS[z_index_ALESS][z_25])/2.5) / (10**((23.9-ch2_ALESS[z_index_ALESS][z_25])/2.5)),s=20,facecolors='r',edgecolors='r', alpha=0.5)#label=r'$z>2.5$',alpha=0.5)
pl.errorbar(10**((23.9-ch2_ALESS)/2.5)[z_index_ALESS][z_25], 10**(-0.4*(ch4_ALESS-ch2_ALESS))[z_index_ALESS][z_25], xerr=((0.4*np.log(10))**2 * 10**((23.9-ch2_ALESS)/2.5) * np.sqrt(e_ch2_ALESS**2))[z_index_ALESS][z_25], yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_ALESS-ch2_ALESS)) * np.sqrt(e_ch4_ALESS**2 + e_ch2_ALESS**2))[z_index_ALESS][z_25], ls='None',c='r',alpha=0.5)#pl.scatter(10**(-0.4*(ch2_1700[cc_index_1700]-ch1_1700[cc_index_1700])),10**(-0.4*(ch4_1700[cc_index_1700]-ch2_1700[cc_index_1700])),s=80, marker='*', facecolors='k', edgecolors='k', label='HS1700')



#pl.scatter(10**((23.9-ch2_1700[cc_index_1700])/2.5) *fudgefactor,10**((23.9-ch4_1700[cc_index_1700])/2.5) / (10**((23.9-ch2_1700[cc_index_1700])/2.5)),s=250,marker='*',facecolors='k',edgecolors='k',label='HS1700')
vmin=2.2
vmax=2.4
cmap='Spectral_r'#'RdYlBu_r'
# cax=pl.scatter(np.concatenate([10**(-0.4*(ch2_ALESS[z_index_ALESS]-ch1_ALESS[z_index_ALESS])),10**(-0.4*(ch2_1700[cc_index_1700]-ch1_1700[cc_index_1700]))]),np.concatenate([10**(-0.4*(ch4_ALESS[z_index_ALESS]-ch2_ALESS[z_index_ALESS])),10**(-0.4*(ch4_1700[cc_index_1700]-ch2_1700[cc_index_1700]))]),c=np.concatenate([z_ALESS[z_index_ALESS],z_1700[cc_index_1700]]), cmap=pl.cm.get_cmap(cmap),marker='o',s=1,vmin=vmin,vmax=vmax)
# pl.scatter(10**(-0.4*(ch2_ALESS[z_index_ALESS]-ch1_ALESS[z_index_ALESS])),10**(-0.4*(ch4_ALESS[z_index_ALESS]-ch2_ALESS[z_index_ALESS])),s=10,c=z_ALESS[z_index_ALESS],cmap=pl.cm.get_cmap(cmap),label='ALESS SMGs',vmin=vmin,vmax=vmax,zorder=10)
# pl.scatter(10**(-0.4*(ch2_1700[z_1700_noz_IRAC]-ch1_1700[z_1700_noz_IRAC])),10**(-0.4*(ch4_1700[z_1700_noz_IRAC]-ch2_1700[z_1700_noz_IRAC])),s=250, marker='*',edgecolors='k',facecolors='w',zorder=10)
# pl.scatter(10**(-0.4*(ch2_1700[z_1700_cluster_IRAC]-ch1_1700[z_1700_cluster_IRAC])),10**(-0.4*(ch4_1700[z_1700_cluster_IRAC]-ch2_1700[z_1700_cluster_IRAC])),s=250, marker='*',facecolors='g',edgecolors='k',zorder=10)
# pl.scatter(10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC]))[0],10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC]))[0],s=250, marker='*', facecolors='blue',edgecolors='k',zorder=10)
# pl.scatter(10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC]))[1],10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC]))[1],s=250, marker='*', facecolors='red',edgecolors='k',zorder=10)


ax2.scatter(10**((23.9-ch2_1549[z_1549_noz_IRAC])/2.5)*fudgefactor,10**((23.9-ch4_1549[z_1549_noz_IRAC])/2.5) / (10**((23.9-ch2_1549[z_1549_noz_IRAC])/2.5)),s=250, marker='*',edgecolors='k',facecolors='w',zorder=10)
pl.errorbar(10**((23.9-ch2_1549[z_1549_noz_IRAC])/2.5)*fudgefactor,10**((23.9-ch4_1549[z_1549_noz_IRAC])/2.5) / (10**((23.9-ch2_1549[z_1549_noz_IRAC])/2.5)), xerr=((0.4*np.log(10))**2 * 10**((23.9-ch2_1549[z_1549_noz_IRAC])/2.5) * np.sqrt(rand_err_1549_x[z_1549_noz_IRAC]**2)), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1549[z_1549_noz_IRAC]-ch2_1549[z_1549_noz_IRAC])) * np.sqrt(rand_err_1549_y[z_1549_noz_IRAC]**2 + rand_err_1549_y[z_1549_noz_IRAC]**2)), ls='None',c='k',alpha=0.5)

ax2.scatter(10**((23.9-ch2_1549[z_1549_cluster_IRAC])/2.5)*fudgefactor,10**((23.9-ch4_1549[z_1549_cluster_IRAC])/2.5) / (10**((23.9-ch2_1549[z_1549_cluster_IRAC])/2.5)),s=250, marker='*',facecolors='r',edgecolors='k',zorder=10)
pl.errorbar(10**((23.9-ch2_1549[z_1549_cluster_IRAC])/2.5)*fudgefactor,10**((23.9-ch4_1549[z_1549_cluster_IRAC])/2.5) / (10**((23.9-ch2_1549[z_1549_cluster_IRAC])/2.5)), xerr=((0.4*np.log(10))**2 * 10**((23.9-ch2_1549[z_1549_cluster_IRAC])/2.5) * np.sqrt(rand_err_1549_x[z_1549_cluster_IRAC]**2)), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1549[z_1549_cluster_IRAC]-ch2_1549[z_1549_cluster_IRAC])) * np.sqrt(rand_err_1549_y[z_1549_cluster_IRAC]**2 + rand_err_1549_y[z_1549_cluster_IRAC]**2)), ls='None',c='r',alpha=1)

ax2.scatter(10**((23.9-ch2_1549[z_1549_notcluster_IRAC])/2.5)*fudgefactor,10**((23.9-ch4_1549[z_1549_notcluster_IRAC])/2.5) / (10**((23.9-ch2_1549[z_1549_notcluster_IRAC])/2.5)),s=250, marker='*', facecolors='g',edgecolors='k',zorder=10)
pl.errorbar(10**((23.9-ch2_1549[z_1549_notcluster_IRAC])/2.5)*fudgefactor,10**((23.9-ch4_1549[z_1549_notcluster_IRAC])/2.5) / (10**((23.9-ch2_1549[z_1549_notcluster_IRAC])/2.5)), xerr=((0.4*np.log(10))**2 * 10**((23.9-ch2_1549[z_1549_notcluster_IRAC])/2.5) * np.sqrt(rand_err_1549_x[z_1549_notcluster_IRAC]**2)), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1549[z_1549_notcluster_IRAC]-ch2_1549[z_1549_notcluster_IRAC])) * np.sqrt(rand_err_1549_y[z_1549_notcluster_IRAC]**2 + rand_err_1549_y[z_1549_notcluster_IRAC]**2)), ls='None',c='g',alpha=1)

#pl.text(19.3,1.2,ID_1549[z_1549_notcluster_IRAC][0],color='g')
ax2.scatter(10**((23.9-ch2_1549[index_1549_notcluster])/2.5)*fudgefactor,10**((23.9-ch4_1549[index_1549_notcluster])/2.5) / (10**((23.9-ch2_1549[index_1549_notcluster])/2.5)),s=250, marker='*',edgecolors='k',facecolors='k',zorder=10)
pl.errorbar(10**((23.9-ch2_1549[index_1549_notcluster])/2.5)*fudgefactor,10**((23.9-ch4_1549[index_1549_notcluster])/2.5) / (10**((23.9-ch2_1549[index_1549_notcluster])/2.5)), xerr=((0.4*np.log(10))**2 * 10**((23.9-ch2_1549[index_1549_notcluster])/2.5) * np.sqrt(rand_err_1549_x[index_1549_notcluster]**2)), yerr=(0.4*np.log(10))**2 * (10**(-0.4*(ch4_1549[index_1549_notcluster]-ch2_1549[index_1549_notcluster])) * np.sqrt(rand_err_1549_y[index_1549_notcluster]**2 + rand_err_1549_y[index_1549_notcluster]**2)), ls='None',c='k',alpha=1)

# pl.scatter(10**(-0.4*(ch2_1700[z_1700_cluster_IRAC]-ch1_1700[z_1700_cluster_IRAC])),10**(-0.4*(ch4_1700[z_1700_cluster_IRAC]-ch2_1700[z_1700_cluster_IRAC])),s=200, marker='*', c=z_1700[z_1700_cluster_IRAC], cmap=pl.cm.get_cmap(cmap),vmin=vmin,vmax=vmax,edgecolors='k',zorder=10)
# pl.scatter(10**(-0.4*(ch2_1700[z_1700_notcluster_IRAC]-ch1_1700[z_1700_notcluster_IRAC])),10**(-0.4*(ch4_1700[z_1700_notcluster_IRAC]-ch2_1700[z_1700_notcluster_IRAC])),s=200, marker='*', c=z_1700[z_1700_notcluster_IRAC], cmap=pl.cm.get_cmap(cmap), label='HS1700',vmin=vmin,vmax=vmax,edgecolors='k',zorder=10)


#pl.scatter(10**(-0.4*(ch2_1700[z_1700_notcluster]-ch1_1700[z_1700_notcluster])),10**(-0.4*(ch4_1700[z_1700_notcluster]-ch2_1700[z_1700_notcluster])),s=200, marker='*',edgecolors='k',facecolors='k')
# cax=pl.scatter(np.concatenate([10**(-0.4*(ch2_ALESS[z_index_ALESS_1700]-ch1_ALESS[z_index_ALESS_1700])),10**(-0.4*(ch2_1700[z_1700_cluster]-ch1_1700[z_1700_cluster]))]),np.concatenate([10**(-0.4*(ch4_ALESS[z_index_ALESS_1700]-ch2_ALESS[z_index_ALESS_1700])),10**(-0.4*(ch4_1700[z_1700_cluster]-ch2_1700[z_1700_cluster]))]),c=np.concatenate([z_ALESS[z_index_ALESS_1700],z_1700[z_1700_cluster]]), cmap=pl.cm.get_cmap(cmap),marker='o',s=1,vmin=vmin,vmax=vmax)
# pl.scatter(10**(-0.4*(ch2_ALESS[z_index_ALESS_1700]-ch1_ALESS[z_index_ALESS_1700])),10**(-0.4*(ch4_ALESS[z_index_ALESS_1700]-ch2_ALESS[z_index_ALESS_1700])),s=10,c=z_ALESS[z_index_ALESS_1700],cmap=pl.cm.get_cmap(cmap),label='ALESS SMGs',vmin=vmin,vmax=vmax)
# pl.scatter(10**(-0.4*(ch2_1700[z_1700_cluster]-ch1_1700[z_1700_cluster])),10**(-0.4*(ch4_1700[z_1700_cluster]-ch2_1700[z_1700_cluster])),s=200, marker='*', c=z_1700[z_1700_cluster], cmap=pl.cm.get_cmap(cmap), label='HS1700',vmin=vmin,vmax=vmax,edgecolors='k')

# template = np.loadtxt('SED/ALESS_compositeSED_JMS2015_LBOL.dat',unpack=True)
# wavelength = template[0]*1e-4
# S = template[1]



zspace = np.linspace(0.6,6,500)
S_all = [S_1,S_2,S_3]
wavelength_all = [wavelength_1,wavelength_2,wavelength_3]
#label = ['ALMA composite', 'ALESS composite', 'Eyelash']
col = ['b','k','k']
SED_fudgefactor = 30
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
    S = movingaverage(S_all[i],25)
    wavelength = wavelength_all[i]
    # if i == 2:
    #     wavelength/=(1+2.35)
    for z in zspace:
        Sch4 = S[np.where(wavelength*(1+z)>=8.0)[0][0]]
        Sch3 = S[np.where(wavelength*(1+z)>=5.8)[0][0]]
        Sch2 = S[np.where(wavelength*(1+z)>=4.5)[0][0]]
        Sch1 = S[np.where(wavelength*(1+z)>=3.6)[0][0]]

        y.append(Sch4/Sch2)
        x.append(Sch2*1e3 * SED_fudgefactor)

    y = np.array(y)
    x = np.array(x)

    #y = movingaverage(y, 5)

#    pl.scatter(x,y,s=10,marker='.',c=zspace,cmap=pl.cm.get_cmap(cmap),vmin=vmin,vmax=vmax,zorder=10, label=label[i])
    ax2.plot(x,y,c=col[i],ls='-',zorder=0)#label=label[i],zorder=0)
    ax2.plot(x[np.where((zspace<1.5))[0]],y[np.where((zspace<1.5))[0]],c='b',ls='-',zorder=0)
    ax2.plot(x[np.where((zspace>1.5) & (zspace<2.5))[0]],y[np.where((zspace>1.5) & (zspace<2.5))[0]],c='g',ls='-',zorder=0)
    ax2.plot(x[np.where((zspace>2.5))[0]],y[np.where((zspace>2.5))[0]],c='r',ls='-',zorder=0)

# pl.hlines(y=0.85,linestyle='-',color='orange',xmin=1.15,xmax=1.66,zorder=0)
# pl.hlines(y=1.53,linestyle='-',color='orange',xmin=1.15,xmax=1.66,zorder=0)
# #pl.axhline(y=3.57,ls='--',c='k')
# pl.vlines(x=1.15,linestyle='-',color='orange',ymin=0.85,ymax=1.53,zorder=0)
# pl.vlines(x=1.66,linestyle='-',color='orange',ymin=0.85,ymax=1.53,zorder=0)
# #pl.axvline(x=1.82,ls='--',c='k')
# pl.text(1.5,0.7,'2.2<z<2.4',color='orange')

'''
#cbar = fig.colorbar(cax,ticks=range(vmin,vmax+1),label='redshift')#,cmap='Spectral_r',vmin=0,vmax=6, ticks=[0,1,2,3,4,5,6])
cbar = fig.colorbar(cax,ticks=[2.2,2.3,2.4],label='redshift')#,cmap='Spectral_r',vmin=0,vmax=6, ticks=[0,1,2,3,4,5,6])
# cbar_index = []
# for i in range(vmin,vmax+1):
#     cbar_index.append(str(i))
# cbar.ax.set_yticklabels(cbar_index)  # vertically oriented colorbar
cbar.ax.set_yticklabels(['2.2','2.3','2.4'])  # vertically oriented colorbar
# pl.xlim(0.5,2.1)
# pl.ylim(0.3,4.5)
'''
ax2.axhline(y=1.60,linestyle='-',color='k')
ax2.axhline(y=1.96,linestyle='-',color='k')
ax2.text(80,2,r'$2.75<z<2.95$')
# ax2.scatter(10**((23.9-ch2_ALESS)/2.5)[np.where((z_ALESS>2.75) & (z_ALESS<2.95))[0]],10**(-0.4*(ch4_ALESS-ch2_ALESS))[np.where((z_ALESS>2.75) & (z_ALESS<2.95))[0]],s=200,marker='s',facecolors='m',edgecolors='m',zorder=0)
# ax2.scatter(10**((23.9-ch2_ALESS)/2.5)[np.where((z_ALESS>2.95))[0]],10**(-0.4*(ch4_ALESS-ch2_ALESS))[np.where((z_ALESS>2.95))[0]],s=100,marker='D',facecolors='c',edgecolors='c',zorder=1)



#ax2.ylabel(r'$S_{8.0}/S_{4.5}$')
ax2.set_yticklabels([])
ax2.set_xlabel(r'$S_{4.5}$ ($\mu$Jy)')

# pl.xlim(0.5,2.1)
ax2.set_ylim(0.3,4.5)
#pl.legend(fontsize=12,loc=1,ncol=1,frameon=True,numpoints=1,scatterpoints=1)


#pl.xticks([0.6,1,2])
ax2.set_yticks([0.5,1,2,3,4])
# ax.set_yticklabels(['0.5','1','2','3','4'])
#ax.set_xticklabels(['0.6','1','2'])
pl.rcParams['axes.formatter.min_exponent'] = 100


ax2.text(1,texty,'HS1549',color='k',size=16)


ax2.scatter(0,0,s=250,marker='*',edgecolors='k',facecolors='b',label=r'$z<1.5$')
ax2.scatter(0,0,s=250,marker='*',edgecolors='k',facecolors='g',label=r'$1.5<z<2.5$')
ax2.scatter(0,0,s=250,marker='*',edgecolors='k',facecolors='r',label=r'$z>2.5$')
ax2.scatter(0,0,s=250,marker='*',edgecolors='k',facecolors='k',label=r'$z\neq z_{proto}$')
ax2.scatter(0,0,s=250,marker='*',edgecolors='k',facecolors='w',label=r'no $z$')
ax2.legend(fontsize=12,loc=4,ncol=1,frameon=True,numpoints=1,scatterpoints=1)

pl.savefig('../../Figures/Colour/CCCM_17001549.pdf', bbox_inches='tight')
#pl.show()
pl.close()
