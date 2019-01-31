import numpy as np
import pylab as pl
import coords as co
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FixedFormatter
from scipy.optimize import curve_fit as cf
from numpy.random import uniform
from scipy import integrate

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



def sch(S):
    return (3300.0/3.7)*(S/3.7)**(-1.4)*np.exp(-S/3.7)

def S2CLS(S):
    return (7180/2.5)*(S/2.5)**(-1.5)*np.exp(-S/2.5)

from scipy.special import erf
def test_S2CLS(S):
    return -22705.15360*np.exp(-.4000000000*S)/np.sqrt(S)-25452.43730*erf(.6324555320*np.sqrt(S)) + 25870


#############
# HS1549 DATA
#############
RA_map_1549 = co.convHMS('15:51:53.582')
DEC_map_1549 = co.convDMS('19:10:50.30')
data_1549 = open('../h1549/flux_4sig.gaia','r')

ID_1549 = []
RA_1549 = []
DEC_1549 = []
flux_1549 = []
for line in data_1549.readlines()[2:]:
    tmp = line.split()
    if tmp[0] == '1': #SMA resolved source
        ID_1549.append('1.1')
        ID_1549.append('1.2')
        ID_1549.append('1.3')
        RA_1549.append(co.convHMS('15:51:53.8'))
        RA_1549.append(co.convHMS('15:51:53.2'))
        RA_1549.append(co.convHMS('15:51:52.5'))
        DEC_1549.append(co.convDMS('19:11:09.9'))
        DEC_1549.append(co.convDMS('19:10:59.1'))
        DEC_1549.append(co.convDMS('19:11:03.9'))
        flux_1549.append(9.4)
        flux_1549.append(5.6)
        flux_1549.append(8.8)
    else:
        ID_1549.append(tmp[0])
        RA_1549.append(co.convHMS(tmp[1]))
        DEC_1549.append(co.convDMS(tmp[2]))
        flux_1549.append(float(tmp[4])*1.1)
ID_1549 = np.array(ID_1549)
RA_1549 = np.array(RA_1549)
DEC_1549 = np.array(DEC_1549)
flux_1549 = np.array(flux_1549)
data_1549.close()

r_1549 = 60*np.sqrt(((RA_map_1549-RA_1549)*np.cos(19*np.pi/180.0))**2 + (DEC_map_1549 - DEC_1549)**2)
Srange_1549 = np.logspace(np.log10(0.99*min(flux_1549)), np.log10(0.99*max(flux_1549)),7)

N_1549 = []
A_1549 = []
for i in Srange_1549:
    if i == Srange_1549[0]:
        A_1549.append(np.pi*(1.65/60.0)**2)
        N_1549.append(len(np.where(r_1549[np.where(flux_1549>i)[0]] <= 1.65)[0]))
    elif i == Srange_1549[1]:
        A_1549.append(np.pi*(2.375/60.0)**2)
        N_1549.append(len(np.where(r_1549[np.where(flux_1549>i)[0]] <= 2.375)[0]))
    elif i == Srange_1549[2]:
        A_1549.append(np.pi*(3.7/60.0)**2)
        N_1549.append(len(np.where(r_1549[np.where(flux_1549>i)[0]] <= 3.7)[0]))
    else:
        N_1549.append(len(np.where(flux_1549>i)[0]))
        A_1549.append(np.pi*(6.0/60.0)**2) #max area at which a source could be detected
N_1549 = np.array(N_1549)
A_1549 = np.array(A_1549)
counts_1549 = N_1549 / A_1549#(np.pi*(R_1549/60.0)**2)
e_counts_1549 = counts_1549/np.sqrt(N_1549)

print "N_1549 = " + str(N_1549)
print "r_1549 = " + str(np.sqrt(A_1549/np.pi)*60)
print "flux_1549" + str(Srange_1549)

'''NEED TO CORRECT AREA'''
'''NEED CORRECT H1549 INNER DATA'''

#############
# HS1700 DATA
#############
RA_map_1700 = 255.2675900821597
DEC_map_1700 = 64.202893133802775
data_1700 = open('../h1700/flux_4sig.gaia','r')

ID_1700 = []
RA_1700 = []
DEC_1700 = []
flux_1700 = []
for line in data_1700.readlines()[2:]:
    tmp = line.split()
    ID_1700.append(tmp[0])
    RA_1700.append(co.convHMS(tmp[1]))
    DEC_1700.append(co.convDMS(tmp[2]))
    flux_1700.append(float(tmp[4])*1.1)
ID_1700 = np.array(ID_1700)
RA_1700 = np.array(RA_1700)
DEC_1700 = np.array(DEC_1700)
flux_1700 = np.array(flux_1700)
data_1700.close()

r_1700 = 60*np.sqrt(((RA_map_1700-RA_1700)*np.cos(64*np.pi/180.0))**2 + (DEC_map_1700 - DEC_1700)**2)
Srange_1700 = np.logspace(np.log10(0.99*min(flux_1700)), .99*np.log10(max(flux_1700)),11)#np.arange(0.99*min(flux_1700), max(flux_1700),1)

N_1700 = []
A_1700 = []
for i in Srange_1700:
    if i == Srange_1700[0]:
        A_1700.append(np.pi*(3.0/60.0)**2)
        N_1700.append(len(np.where(r_1700[np.where(flux_1700>i)[0]] <= 3.0)[0]))
    elif i == Srange_1700[1]:
        A_1700.append(np.pi*(3.83/60.0)**2)
        N_1700.append(len(np.where(r_1700[np.where(flux_1700>i)[0]] <= 3.83)[0]))
    else:
        N_1700.append(len(np.where(flux_1700>i)[0]))
        A_1700.append(np.pi*(np.max(r_1700)/60.0)**2)
N_1700 = np.array(N_1700)
A_1700 = np.array(A_1700)
counts_1700 = N_1700 / A_1700#(np.pi*(R_1700/60.0)**2)
e_counts_1700 = counts_1700/np.sqrt(N_1700)

Srange_1700 = np.concatenate([Srange_1700[:-2],[Srange_1700[-1]]])
counts_1700 = np.concatenate([counts_1700[:-2],[counts_1700[-1]]])
e_counts_1700 = np.concatenate([e_counts_1700[:-2],[e_counts_1700[-1]]])

#get S2CLS function
S2CLS_function = []
S2CLS_space = np.logspace(np.log10(2),np.log10(20),100)
for i in S2CLS_space:
    S2CLS_function.append(integrate.quad(S2CLS,i,np.inf)[0])
S2CLS_function = np.array(S2CLS_function)

#print 'HS1549 overdensity (flux, overdensity, max, min):'
#for i in range(len(Srange_1549)):
#    print round(Srange_1549[i],2), round(counts_1549[i] / integrate.quad(S2CLS,Srange_1549[i],np.inf)[0],1), round((e_counts_1549[i]) / integrate.quad(S2CLS,Srange_1549[i],np.inf)[0],1)
#print 'HS1700 overdensity (flux, overdensity, max, min):'
#for i in range(len(Srange_1700)):
#    print round(Srange_1700[i],2), round(counts_1700[i] / integrate.quad(S2CLS,Srange_1700[i],np.inf)[0],1), round((e_counts_1700[i]) / integrate.quad(S2CLS,Srange_1700[i],np.inf)[0],1)


######
# PLOT
######
fig = pl.figure(figsize=(5,8))
pl.subplots_adjust(hspace=0,wspace=0)
pl.rc('font',size=15)
pl.rc('mathtext', default='regular')
ax = fig.add_subplot(2,1,1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([2,3,4,5,6,7,8,9,10,20])
#ax.set_xticklabels(['2','3','4','5','6','7','8','9','10','20'])
pl.ylabel(r'N(>S$_{850}$) (deg$^{-2}$)')
#pl.xlabel('S$_{850}$ (mJy)')

#        ax.scatter(R_850, p_850_inner,s=60,marker='o',facecolors='w',edgecolors='magenta')
ax.scatter(Srange_1549, counts_1549, s=120, marker='*', facecolors='r', edgecolors='r', label='HS1549')
ax.scatter(Srange_1700, counts_1700, s=120, marker='*', facecolors='b', edgecolors='b', label='HS1700')
ax.plot(Srange_1549, counts_1549, 'r-')
ax.plot(Srange_1700, counts_1700, 'b-')
yerr_1549 = [[],[]]
yerr_1700 = [[],[]]
for i in range(len(e_counts_1549)):
    if e_counts_1549[i] != counts_1549[i]:
        yerr_1549[0].append(e_counts_1549[i])
        yerr_1549[1].append(e_counts_1549[i])
    else:
        yerr_1549[0].append(0)
        yerr_1549[1].append(e_counts_1549[i])
for i in range(len(e_counts_1700)):
    if e_counts_1700[i] != counts_1700[i]:
        yerr_1700[0].append(e_counts_1700[i])
        yerr_1700[1].append(e_counts_1700[i])
    else:
        yerr_1700[0].append(0)
        yerr_1700[1].append(e_counts_1700[i])
pl.errorbar(Srange_1549, counts_1549, yerr = yerr_1549, ls='.', c='r')
pl.errorbar(Srange_1700, counts_1700, yerr = yerr_1700, ls='.', c='b')
ax.errorbar(Srange_1549[-1], counts_1549[-1], yerr=[[counts_1549[-1]*.6],[0]], lolims=True, ls='.', c='r')
ax.errorbar(Srange_1700[-2:], counts_1700[-2:], yerr=[counts_1700[-2:]*.6,[0,0]], lolims=True, ls='.', c='b')

#pl.plot(np.logspace(0,1.3,100), sch(np.logspace(0,1.3,100)), 'k--') #CASEY 2013
#pl.plot(np.logspace(0,1.3,100), S2CLS(np.logspace(0,1.3,100)), 'k-') #S2CLS
#
#pl.plot(np.logspace(0,1.3,100), 2.1*S2CLS(np.logspace(0,1.3,100)), 'g-') #S2CLS fit
#pl.plot(np.logspace(0,1.3,100), 1.8*sch(np.logspace(0,1.3,100)), 'g--') #CASEY 2013 fit

#pl.plot(np.linspace(1,6.045149145,10),6827./np.linspace(1,6.045149145,10)**(7./5.), 'k--',label='SHADES') #SHADES
#pl.plot(np.linspace(6.045149145,20,10),8.791217434E5/np.linspace(6.045149145,20,10)**(4.1), 'k--') #SHADES

ax.scatter(np.arange(3,16), [1012,508,272,152,85,47,26,15,9,5.5,3.2,2.4,1.8], s=80, marker='o', facecolors='w', edgecolors='k', label='S2CLS')
pl.errorbar(np.arange(3,16), [1012,508,272,152,85,47,26,15,9,5.5,3.2,2.4,1.8], [19.6,12.3,8.5,6.2,4.7,3.6,2.8,2.2,1.8,1.5,1.2,1.1,1], ls='.',c='k')
#pl.plot(np.arange(3,16), [1012-19.2,508-12,272-8.2,152-6,85-4.4,47-3.3,26-2.5,15-1.9,9-1.5,5.5-1.2,3.2-0.9,2.4-0.8,1.8-0.7], 'k-')

pl.plot(S2CLS_space, 1.2*S2CLS_function, 'k-')

pl.xlim(2,20)
pl.ylim(1,2E4)
pl.text(2.2,2,'Entire field',size=16)
ax.legend(fontsize=15,loc=1,ncol=1,frameon=False,numpoints=1,scatterpoints=1)





'''CORE'''




#############
# HS1549 DATA
#############
RA_map_1549 = co.convHMS('15:51:53.582')
DEC_map_1549 = co.convDMS('19:10:50.30')
data_1549 = open('../h1549/flux_4sig.gaia','r')

ID_1549 = []
RA_1549 = []
DEC_1549 = []
flux_1549 = []
for line in data_1549.readlines()[2:]:
    tmp = line.split()
    if tmp[0] == '1': #SMA resolved source
        ID_1549.append('1.1')
        ID_1549.append('1.2')
        ID_1549.append('1.3')
        RA_1549.append(co.convHMS('15:51:53.8'))
        RA_1549.append(co.convHMS('15:51:53.2'))
        RA_1549.append(co.convHMS('15:51:52.5'))
        DEC_1549.append(co.convDMS('19:11:09.9'))
        DEC_1549.append(co.convDMS('19:10:59.1'))
        DEC_1549.append(co.convDMS('19:11:03.9'))
        flux_1549.append(9.4)
        flux_1549.append(5.6)
        flux_1549.append(8.8)
    else:
        ID_1549.append(tmp[0])
        RA_1549.append(co.convHMS(tmp[1]))
        DEC_1549.append(co.convDMS(tmp[2]))
        flux_1549.append(float(tmp[4])*1.1)
ID_1549 = np.array(ID_1549)
RA_1549 = np.array(RA_1549)
DEC_1549 = np.array(DEC_1549)
flux_1549 = np.array(flux_1549)
data_1549.close()

r_1549 = 60*np.sqrt(((RA_map_1549-RA_1549)*np.cos(19*np.pi/180.0))**2 + (DEC_map_1549 - DEC_1549)**2)

ID_1549 = np.array(ID_1549)[np.where(r_1549<=1.5)[0]]
RA_1549 = np.array(RA_1549)[np.where(r_1549<=1.5)[0]]
DEC_1549 = np.array(DEC_1549)[np.where(r_1549<=1.5)[0]]
flux_1549 = np.array(flux_1549)[np.where(r_1549<=1.5)[0]]

r_1549 = r_1549[np.where(r_1549<=1.5)[0]]


Srange_1549 = np.logspace(np.log10(0.99*min(flux_1549)),np.log10(0.99*max(flux_1549)),4)

N_1549 = []
for i in Srange_1549:
    N_1549.append(len(np.where(flux_1549>i)[0]))
N_1549 = np.array(N_1549)
A_1549 = np.pi*(1.5/60.0)**2
counts_1549 = N_1549 / A_1549
e_counts_1549 = counts_1549/np.sqrt(N_1549)

'''NEED TO CORRECT AREA'''
'''NEED CORRECT H1549 INNER DATA'''

#############
# HS1700 DATA
#############
RA_map_1700 = 255.2675900821597
DEC_map_1700 = 64.202893133802775
data_1700 = open('../h1700/flux_4sig.gaia','r')

ID_1700 = []
RA_1700 = []
DEC_1700 = []
flux_1700 = []
for line in data_1700.readlines()[2:]:
    tmp = line.split()
    ID_1700.append(tmp[0])
    RA_1700.append(co.convHMS(tmp[1]))
    DEC_1700.append(co.convDMS(tmp[2]))
    flux_1700.append(float(tmp[4])*1.1)
ID_1700 = np.array(ID_1700)
RA_1700 = np.array(RA_1700)
DEC_1700 = np.array(DEC_1700)
flux_1700 = np.array(flux_1700)
data_1700.close()

r_1700 = 60*np.sqrt(((RA_map_1700-RA_1700)*np.cos(64*np.pi/180.0))**2 + (DEC_map_1700 - DEC_1700)**2)

ID_1700 = np.array(ID_1700)[np.where(r_1700<=1.5)[0]]
RA_1700 = np.array(RA_1700)[np.where(r_1700<=1.5)[0]]
DEC_1700 = np.array(DEC_1700)[np.where(r_1700<=1.5)[0]]
flux_1700 = np.array(flux_1700)[np.where(r_1700<=1.5)[0]]

r_1700 = r_1700[np.where(r_1700<=1.5)[0]]

Srange_1700 =  Srange_1700 = [2.571657,3.8,5.2,6.73476804] #np.logspace(np.log10(0.99*min(flux_1700)),np.log10(0.99*max(flux_1700)),3)

N_1700 = []
for i in Srange_1700:
    N_1700.append(len(np.where(flux_1700>i)[0]))
N_1700 = np.array(N_1700)
A_1700 = np.pi*(1.5/60.0)**2
counts_1700 = N_1700 / A_1700
e_counts_1700 = counts_1700/np.sqrt(N_1700)

Srange_1700 = np.concatenate([Srange_1700[:-2],[Srange_1700[-1]]])
counts_1700 = np.concatenate([counts_1700[:-2],[counts_1700[-1]]])
e_counts_1700 = np.concatenate([e_counts_1700[:-2],[e_counts_1700[-1]]])

#print 'HS1549 overdensity (flux, overdensity, max, min):'
#for i in range(len(Srange_1549)):
#    print round(Srange_1549[i],2), round(counts_1549[i] / integrate.quad(S2CLS,Srange_1549[i],np.inf)[0],1), round((e_counts_1549[i]) / integrate.quad(S2CLS,Srange_1549[i],np.inf)[0],1)
#print 'HS1700 overdensity (flux, overdensity, max, min):'
#for i in range(len(Srange_1700)):
#    print round(Srange_1700[i],2), round(counts_1700[i] / integrate.quad(S2CLS,Srange_1700[i],np.inf)[0],1), round((e_counts_1700[i]) / integrate.quad(S2CLS,Srange_1700[i],np.inf)[0],1)



ax2 = fig.add_subplot(2,1,2)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xticks([2,3,4,5,6,7,8,9,10,20])
ax2.set_xticklabels(['2','3','4','5','6','7','8','9','10','20'])
pl.ylabel(r'N(>S$_{850}$) (deg$^{-2}$)')
pl.xlabel('S$_{850}$ (mJy)')

#        ax.scatter(R_850, p_850_inner,s=60,marker='o',facecolors='w',edgecolors='magenta')
ax2.scatter(Srange_1549, counts_1549, s=120, marker='*', facecolors='r', edgecolors='r', label='HS1549')
ax2.scatter(Srange_1700, counts_1700, s=120, marker='*', facecolors='b', edgecolors='b', label='HS1700')
ax2.plot(Srange_1549, counts_1549, 'r-')
ax2.plot(Srange_1700, counts_1700, 'b-')
yerr_1549 = [[],[]]
yerr_1700 = [[],[]]
for i in range(len(e_counts_1549)):
    if e_counts_1549[i] != counts_1549[i]:
        yerr_1549[0].append(e_counts_1549[i])
        yerr_1549[1].append(e_counts_1549[i])
    else:
        yerr_1549[0].append(0)
        yerr_1549[1].append(e_counts_1549[i])
for i in range(len(e_counts_1700)):
    if e_counts_1700[i] != counts_1700[i]:
        yerr_1700[0].append(e_counts_1700[i])
        yerr_1700[1].append(e_counts_1700[i])
    else:
        yerr_1700[0].append(0)
        yerr_1700[1].append(e_counts_1700[i])
pl.errorbar(Srange_1549, counts_1549, yerr = yerr_1549, ls='.', c='r')
pl.errorbar(Srange_1700, counts_1700, yerr = yerr_1700, ls='.', c='b')
ax2.errorbar(Srange_1549[-1], counts_1549[-1], yerr=[[counts_1549[-1]*.6],[0]], lolims=True, ls='.', c='r')
ax2.errorbar(Srange_1700[-1], counts_1700[-1], yerr=[[counts_1700[-1]*.6],[0]], lolims=True, ls='.', c='b')
ax2.scatter(np.arange(3,16), [1012,508,272,152,85,47,26,15,9,5.5,3.2,2.4,1.8], s=80, marker='o', facecolors='w', edgecolors='k', label='S2CLS')
pl.errorbar(np.arange(3,16), [1012,508,272,152,85,47,26,15,9,5.5,3.2,2.4,1.8], [19.6,12.3,8.5,6.2,4.7,3.6,2.8,2.2,1.8,1.5,1.2,1.1,1], ls='.',c='k')



pl.plot(S2CLS_space, 1.2*S2CLS_function, 'k-')

pl.text(2.2,2,'Core region',size=16)
#ax2.legend(fontsize=15,loc=1,ncol=1,frameon=False,numpoints=1,scatterpoints=1)
pl.xlim(2,20)
pl.ylim(1,2E4)

pl.savefig('../Figures/flux_counts_both.pdf',bbox_inches='tight')
#pl.show()
pl.close()



