import numpy as np
import matplotlib.pyplot as pl
import coords as co
from scipy.optimize import curve_fit as cf

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

def SFR_low(x):
    SFR = 7.8E-10 * x
    return SFR

def SFR_high(x):
    SFR = 7.8E-10*x * (7.76E-11*x)**0.048
    return SFR

#h1700 24um data
data_1700 = np.loadtxt('../../Protocluster_paper/H1700/Data/h1700_24um_cat.txt',unpack=True)
bins_1700 = 10**data_1700[0] #log(F_24um)
sub1_1700 = data_1700[1] #these are the red and green subpopulations that are shown in the figure
sub2_1700 = data_1700[2] #but don't use these or you'll double count
Rband_1700 = data_1700[3] #undetected in 24um have estimates of their 24um flux from the R-band mag
all_1700 = data_1700[4] #everything detected at 24um
width_1700 = np.logspace(0.05+0.15,3.5+0.15,len(data_1700[0])) - bins_1700

#850um data
ID_850 = []
RA_850 = []
DEC_850 = []
cat_850 = open('../h1549/4sig.gaia','r')
for line in cat_850.readlines()[2:]:
    tmp = line.split()
    ID_850.append(tmp[0])
    RA_850.append(co.convHMS(tmp[1]))
    DEC_850.append(co.convDMS(tmp[2]))
ID_850 = np.array(ID_850)
RA_850 = np.array(RA_850)
DEC_850 = np.array(DEC_850)
cat_850.close()

#proto-z LBG data for HS1549
ID = []
RA = []
DEC = []
RA_1549 = co.convHMS('15:51:52.604')
DEC_1549 = co.convDMS('19:11:22.90')
flux = []
z = []
cat = open('../h1549/MIPS/h1549_24um_proto_cat.gaia','r')
for line in cat.readlines()[2:]:
    tmp = line.split()
    ID.append(tmp[0])
    RA.append(float(tmp[1]))
    DEC.append(float(tmp[2]))
    flux.append(float(tmp[3]))
    z.append(float(tmp[4]))
ID = np.array(ID)
RA = np.array(RA)
DEC = np.array(DEC)
'''theres a difference adding here versus adding in data before...'''
flux = np.array(flux)*1000 +1000*0.009*1.6 #convert to uJy and add 9*1.6uJy because map isnt flat
z = np.array(z)
cat.close()

#append all 850 ID'd LBG's and obstructed LBGs
ID_smg = []
smg = []
ID_obstruct = []
obstruct = []
false_AGN = []

ID_obstruct = ['GNB1617','GNB2267','GNB2044','GNB1521','GNB4231','GNB1521','GNB','GNB4724','GNB2614','GNB2249','GNB3532','GNB2753','GNB4','GNB4231','GNB1300a','GNB126','GNB2044']

false_AGN.append(np.where(ID == 'GNB4950a')[0][0])
false_AGN.append(np.where(ID == 'GNB4950b')[0][0])

obstruct.append(np.where(ID == 'GNB1617')[0][0])
obstruct.append(np.where(ID == 'GNB2267')[0][0])
obstruct.append(np.where(ID == 'GNB2044')[0][0])
obstruct.append(np.where(ID == 'GNB1521')[0][0])
obstruct.append(np.where(ID == 'GNB4231')[0][0])
obstruct.append(np.where(ID == 'GNB1521')[0][0])
obstruct.append(np.where(ID == 'GNB520')[0][0])
obstruct.append(np.where(ID == 'GNB4724')[0][0])
obstruct.append(np.where(ID == 'GNB2614')[0][0])
obstruct.append(np.where(ID == 'GNB2249')[0][0])
obstruct.append(np.where(ID == 'GNB3532')[0][0])
obstruct.append(np.where(ID == 'GNB2753')[0][0])
obstruct.append(np.where(ID == 'GNB4')[0][0])
obstruct.append(np.where(ID == 'GNB4231')[0][0])
obstruct.append(np.where(ID == 'GNB1300a')[0][0])
obstruct.append(np.where(ID == 'GNB126')[0][0])
obstruct.append(np.where(ID == 'GNB2044')[0][0])
obstruct.append(np.where(ID == 'GNB4022b')[0][0])
obstruct.append(np.where(ID == 'GNB4049')[0][0])

'''
for i in range(len(RA_850)):
    for j in range(len(RA)):
        if i == 0:
            if 3600*np.sqrt(((RA_850[i]-RA[j])*np.cos(19*np.pi/180.0))**2 + (DEC_850[i] - DEC[j])**2) <= 12:
                print ID_850[i],ID[j]
                ID_smg.append(ID_850[i])
                smg.append(j)
        else:
            if 3600*np.sqrt(((RA_850[i]-RA[j])*np.cos(19*np.pi/180.0))**2 + (DEC_850[i] - DEC[j])**2) <= 7.5:
                print ID_850[i],ID[j]
                ID_smg.append(ID_850[i])
                smg.append(j)
'''

ID_smg = ['MD17','D14','Q1549','GNB3999']

smg.append(np.where(ID == 'MD17')[0][0])
smg.append(np.where(ID == 'D14')[0][0])
smg.append(np.where(ID == 'Q1549')[0][0])
smg.append(np.where(ID == 'GNB3999')[0][0])


#all LBGs not substructed nor 850 IDs
good = [i for i in range(len(RA)) if i not in smg+obstruct+false_AGN]

flux_good = flux[good][np.where(flux[good]>0)[0]]
flux_smg = flux[smg]
flux_obstruct = flux[obstruct]
flux_zero = flux[np.where(flux<0)[0]]
flux_false = flux[false_AGN]

z_good = z[good][np.where(flux[good]>0)[0]]
r_1549 = 60*np.sqrt(((RA_1549 - RA)*np.cos(19.*np.pi/180))**2 + (DEC_1549 - DEC)**2)
r_1549_good = r_1549[good][np.where(flux[good]>0)[0]]

'''1700 time'''
cat_op = open('../../Protocluster_paper/H1700/Data/Catalogues_new/q1700_all_radec_python.gaia','r')
ID_op = []
RA_op = []
DEC_op = []
z_1700 = []
for line in cat_op.readlines()[2:]:
    tmp = line.split()
    ID_op.append(tmp[0])
    RA_op.append(co.convHMS(tmp[1]))
    DEC_op.append(co.convDMS(tmp[2]))
    z_1700.append(float(tmp[4]))
ID_op = np.array(ID_op)
RA_op = np.array(RA_op)
DEC_op = np.array(DEC_op)
z_1700 = np.array(z_1700)

RA_c_1700,DEC_c_1700 = np.mean(RA_op[np.where(abs(z_1700-2.3)<=.05)[0]]),np.mean(DEC_op[np.where(abs(z_1700-2.3)<=.05)[0]])

r_1700 = 60*np.sqrt(((RA_c_1700 - RA_op)*np.cos(64.*np.pi/180))**2 + (DEC_c_1700 - DEC_op)**2)
z_1700_proto = z_1700[np.where(abs(z_1700-2.3) <=0.05)[0]]
r_1700_proto = r_1700[np.where(abs(z_1700-2.3) <=0.05)[0]]
ID_op_proto = ID_op[np.where(abs(z_1700-2.3) <=0.05)[0]]

r_1700_proto_core = r_1700_proto[np.where(r_1700_proto<=1.5)[0]]
ID_op_proto_core = ID_op_proto[np.where(r_1700_proto<=1.5)[0]]

#flux_1700 = np.array([0,0,0.0112353,0.0476978,0.0490396,0.0490396,0,0,0.0112353,0.0281812,0.0112353])*1000 +1000*0.009*1.6
flux_1700 = np.array([ 0,0,25.6353,62.0978,63.4396,63.4396,0,0,25.6353,42.5812,25.6353])
DL_1700 = 18773.7*3.08567758E22 #meters
Hz = 2.998E8/24E-6
SFR_1700_core = sum(SFR_high(4*np.pi*(DL_1700**2)*(flux_1700*1E-6*1E-26)/3.846E26 * Hz))


#24um histogram
fig = pl.figure(figsize=(6,6))
fig.subplots_adjust(hspace=0,wspace=0)
pl.rc('font',size=15)
pl.rc('mathtext', default='regular')
ax1 = fig.add_subplot(2,1,1)

pl.ylabel('N')
pl.setp(ax1.get_xticklabels(), visible=False)

pl.plot(100,10,'w.',label='HS1549')
ax1.hist(flux_good, bins = np.logspace(np.log10(1),np.log10(max(flux)+10),25),color='blue')
ax1.hist(flux_good[np.where(flux_good<=10)[0]], bins = np.logspace(np.log10(1),np.log10(max(flux)+10),25),color='purple')
ax1.hist(flux_obstruct, bins = np.logspace(np.log10(1),np.log10(max(flux)+10),25),color='magenta',alpha=0)
ax1.hist(flux_obstruct, bins = np.logspace(np.log10(1),np.log10(max(flux)+10),25),facecolor='None',hatch='/')
ax1.hist(flux_obstruct, bins = np.logspace(np.log10(1),np.log10(max(flux)+10),25),facecolor='None',hatch='\\')
ax1.hist(flux_false, bins = np.logspace(np.log10(1),np.log10(max(flux)+10),25),color='red')
#ax1.hist(flux_smg, bins = np.logspace(np.log10(1),np.log10(max(flux)+10),25),color='green')
#ax1.hist(flux[index], bins=np.logspace(np.log10(1.0),np.log10(max(flux)+10),25),color='purple')
#data_24 = ax1.hist(flux[index][high], bins=np.logspace(np.log10(1.0),np.log10(max(flux)+10),25))
#data_24 = ax1.hist(flux[good][high], bins=np.logspace(np.log10(1.0),np.log10(max(flux)+10),25))
pl.legend(fontsize=15,loc=0,frameon=False,numpoints=1)
pl.gca().set_xscale("log")
#ax1.set_yticks([0,1,2,3,4])
#ax1.set_yticklabels([0,1,2,3,4])
pl.xlim(1,3E3)

ax2 = fig.add_subplot(2,1,2)
pl.ylabel('N')
pl.xlabel(r'$S_{24 \mu m}$ ($\mu$Jy)')
pl.plot(1,1,'w.',label='HS1700')
ax2.bar(bins_1700[4:-4],all_1700[4:-4],color='blue',width=width_1700[4:-4])
ax2.bar(bins_1700,Rband_1700,color='purple',width=width_1700)
ax2.bar(bins_1700,sub1_1700,color='red',width=width_1700)
#ax2.bar(bins_1700,sub2_1700,color='green',width=width_1700)

pl.gca().set_xscale("log")
ax2.set_yticks([0,5,10,15])
ax2.set_yticklabels([0,5,10,15])
ax2.set_xticks([1,10,100,1000,10000])
ax2.set_xticklabels([1,10,100,1000,10000])
pl.xlim(1,2E3)

pl.legend(fontsize=15,loc=0,frameon=False,numpoints=1)

pl.savefig('../Figures/24um.pdf')#,bbox_inches='tight')
pl.savefig('../Figures/24um.png',dpi=200)#,bbox_inches='tight')
#pl.show()
pl.close()


#SFR calculations

DL_1700 = 18773.7*3.08567758E22 #meters
DL_1549 = 24366.1*3.08567758E22 #meters
Hz = 2.998E8/24E-6

S_1700 = []
for i in range(len(all_1700[4:-4])):
    tmp = np.ones(int(all_1700[4:-4][i]))*bins_1700[4:-4][i]
    for j in tmp:
        S_1700.append(j)
S_1700 = np.array(S_1700)

# L =4*pi*D_L^2*S
L_1549_good = 4*np.pi*(DL_1549**2)*(flux_good[np.where(flux_good>=10)[0]]*1E-6*1E-26)/3.846E26 * Hz #only above RMS
dL_1549_good = np.array(len(L_1549_good)*[4*np.pi*(DL_1549**2)*(3*0.008*1E-3*1E-26)/3.846E26 * Hz]) #2RMS

L_1549_good_low = 4*np.pi*(DL_1549**2)*(flux_good[np.where(flux_good<10)[0]]*1E-6*1E-26)/3.846E26 * Hz #only above RMS
dL_1549_good_low = np.array(len(L_1549_good_low)*[4*np.pi*(DL_1549**2)*(3*0.008*1E-3*1E-26)/3.846E26 * Hz]) #2RMS

L_1700 = 4*np.pi*(DL_1700**2)*(S_1700*1E-6*1E-26)/3.846E26 * Hz #in 24um
dL_1700 = np.array(len(L_1700)*[4*np.pi*(DL_1700**2)*(3*0.0064*1E-3*1E-26)/3.846E26 * Hz]) #2RMS


#SFR from Reike+2009
SFR_1549_good = []
dSFR_1549_good = []
for i in range(len(L_1549_good)):
    if L_1549_good[i]>=6E8 and L_1549_good[i]<1.3E10:
        SFR_1549_good.append(SFR_low(L_1549_good[i]))
        dSFR_1549_good.append(SFR_low(dL_1549_good[i]))
    else:
        SFR_1549_good.append(SFR_high(L_1549_good[i]))
        dSFR_1549_good.append(SFR_high(dL_1549_good[i]))
SFR_1549_good = np.array(SFR_1549_good)
dSFR_1549_good = np.array(dSFR_1549_good)

SFR_1549_good_low = []
dSFR_1549_good_low = []
for i in range(len(L_1549_good_low)):
    if L_1549_good_low[i]>=6E8 and L_1549_good_low[i]<1.3E10:
        SFR_1549_good_low.append(SFR_low(L_1549_good_low[i]))
        dSFR_1549_good_low.append(SFR_low(dL_1549_good_low[i]))
    else:
        SFR_1549_good_low.append(SFR_high(L_1549_good_low[i]))
        dSFR_1549_good_low.append(SFR_high(dL_1549_good_low[i]))
SFR_1549_good_low = np.array(SFR_1549_good_low)
dSFR_1549_good_low = np.array(dSFR_1549_good_low)

SFR_1700_all = []
dSFR_1700_all = []
for i in range(len(L_1700)):
    if L_1700[i]>=6E8 and L_1700[i]<1.3E10:
        SFR_1700_all.append(SFR_low(L_1700[i]))
        dSFR_1700_all.append(SFR_low(dL_1700[i]))
    else:
        SFR_1700_all.append(SFR_high(L_1700[i]))
        dSFR_1700_all.append(SFR_high(dL_1700[i]))
SFR_1700_all = np.array(SFR_1700_all)
dSFR_1700_all = np.array(dSFR_1700_all)


SFR_1549_all = np.concatenate([SFR_1549_good,SFR_1549_good_low])
dSFR_1549_all = np.concatenate([dSFR_1549_good,dSFR_1549_good_low])


SFR_1549 = 2*(sum(SFR_1549_good) + sum(SFR_1549_good_low) )#remove IR-obscured and AGN #+ 10*len(flux_zero) + 10*len(flux_obstruct) + 150*len(false_AGN)) #2x to adjust for position of PAH emission line at z=2.85
SFR_1700 = sum(SFR_1700_all) #remove IR-obscured and AGN + 10.0*(sum(Rband_1700)) + 150*sum(sub2_1700)

boost_1549 = 75*(8.28704 + 8.9232 + 5.36911) + 75*6.37961*1.1 #Barger+2014 SFR for 3 inside HS1549_1 (from SMA fluxes) + HS1549_2*1.1 to compensate for matched filter
boost_1700 = 75*(10.9+6.3+6.6+3.7) #all SMGs with spec-z 4,5,7,17, already matched filter corrected
SFR_1549_boost = SFR_1549 + boost_1549
SFR_1549_boost_conservative = SFR_1549/2.0 + boost_1549
SFR_1700_boost = SFR_1700 + boost_1700

dboost_1549 = 17*(8.28704 + 8.9232 + 5.36911) + 17*6.37961*1.1
dboost_1700 = 17*(10.9+6.3+6.6+3.7)

#uncertainties April 6,2016 notebook
dSFR_1549_m = (SFR_1549/(1.5*1))*np.sqrt((np.sqrt(sum(dSFR_1549_all**2))/SFR_1549)**2 + 0.5**2)
dSFR_1549_m_boost_conservative = np.sqrt((SFR_1549_boost_conservative/(1.5*1.0)*np.sqrt((75*4*0.6/SFR_1549_boost_conservative)**2 + 0.5**2))**2 + dSFR_1549_m**2)
dSFR_1549_m_boost = np.sqrt((SFR_1549_boost/(1.5*1.0)*np.sqrt((75*4*0.6/SFR_1549_boost)**2 + 0.5**2))**2 + dSFR_1549_m**2)




dSFR_1700_m = (SFR_1700/1.0)*np.sqrt((np.sqrt(sum(dSFR_1700_all**2))/SFR_1700)**2 + 0.5**2)
#np.sqrt((np.sqrt(sum(dSFR_1700_all**2))/SFR_1700)**2 + (0.5)**2)*SFR_1700
dSFR_1700_m_boost = np.sqrt((SFR_1700_boost/1.0*np.sqrt((200*4*0.5/SFR_1700_boost)**2 + 0.5**2))**2 + dSFR_1700_m**2)

N_1700 = sum(all_1700[4:-4])+sum(Rband_1700)
N_1549 = len(SFR_1549_all)+len(flux_zero)

print "\n" + "24um/R-band SFRs:"
print "HS1549: " + str(int(SFR_1549)) + " " + str(int(SFR_1549_boost))
print "HS1700: " + str(int(SFR_1700)) + " " + str(int(SFR_1700_boost)) + "\n"

#SFR values from stack and boost values from stack plus 200*S_850 (seen above)
SFR_1549_S2 = 11100.0
dSFR_1549_S2 = 1100.0
SFR_1549_boost_S2 = 11100.0 + boost_1549
dSFR_1549_boost_S2 = np.sqrt(690**2 + 1100**2)

SFR_1700_S2 = 1400.0
dSFR_1700_S2 = 400.0
SFR_1700_boost_S2 = 1400.0 + boost_1700
dSFR_1700_boost_S2 = np.sqrt(430**2 + 400**2)

dSFR_1549_m_S2 = (SFR_1549_S2/1.5*1)*np.sqrt((dSFR_1549_S2/SFR_1549_S2)**2 + 0.5**2)
dSFR_1549_m_boost_S2 = np.sqrt((SFR_1549_boost_S2/(1.5*1.0)*np.sqrt((200*4*0.6/SFR_1549_boost_S2)**2 + 0.5**2))**2 + dSFR_1549_m_S2**2)

dSFR_1700_m_S2 = (SFR_1700_boost_S2/1.0)*np.sqrt((dSFR_1700_S2/SFR_1700_S2)**2 + 0.5**2)
dSFR_1700_m_boost_S2 = np.sqrt((SFR_1700_boost_S2/(1.0)*np.sqrt((200*4*0.5/SFR_1700_boost_S2)**2 + 0.5**2))**2 + dSFR_1700_m_S2**2)

dSFR_1549 = np.sqrt(sum(dSFR_1549_all**2))
dSFR_1549_boost = np.sqrt((dSFR_1549/SFR_1549)**2 + (dboost_1549/boost_1549)**2) * SFR_1549_boost #np.sqrt((4*200*0.6)**2 + dSFR_1549**2)

dSFR_1700 = np.sqrt(sum(dSFR_1700_all**2))
dSFR_1700_boost = np.sqrt((dSFR_1700/SFR_1700)**2 + (dboost_1549/boost_1549)**2) * SFR_1700_boost #np.sqrt((4*200*0.5)**2 + dSFR_1700**2)

print "24um/R-band SFRDs:"
print "HS1549: " + str(int(SFR_1549_boost / 4.2)) + " +/- " + str(int(SFR_1549_boost/4.2 * np.sqrt((dSFR_1549_boost / SFR_1549_boost)**2 + 0.1**1)))
print "HS1700: " + str(int(SFR_1700_boost / 4.2)) + " +/- " + str(int(SFR_1700_boost/4.2 * np.sqrt((dSFR_1700_boost / SFR_1700_boost)**2 + 0.1**1))) + "\n"


print "SCUBA-2 SFRs:"
print "HS1549: " + str(int(SFR_1549_S2)) + " " + str(int(SFR_1549_boost_S2))
print "HS1700: " + str(int(SFR_1700_S2)) + " " + str(int(SFR_1700_boost_S2)) + "\n"

#volume for HS1700 = 4.2 Mpc^3 (from Kato+2016)
#use same volume for HS1549 size dz doesnt change much between z~2-4 and dz_1549 ~ dz_1700 = 0.005
#assumed 10% error on volume
print "SCUBA-2 SFRDs:"
print "HS1549: " + str(int(SFR_1549_boost_S2 / 4.2)) + " +/- " + str(int(SFR_1549_boost_S2/4.2 * np.sqrt((dSFR_1549_boost_S2 / SFR_1549_boost_S2)**2 + 0.1**1)))
print "HS1700: " + str(int(SFR_1700_boost_S2 / 4.2)) + " +/- " + str(int(SFR_1700_boost_S2/4.2 * np.sqrt((dSFR_1700_boost_S2 / SFR_1700_boost_S2)**2 + 0.1**1))) + "\n"


SFR_XCSJ2215 = [130,140,140,160,360,130,150,180,180,160,160]

zspace = np.linspace(0,5.5,200)

fig = pl.figure()
ax = pl.gca()
pl.rc('font',size=16)
pl.rc('mathtext', default='regular')
#ax.set_xlim(0,max(zspace))
#ax.set_ylim(0.1,1E6)
xlim = ax.set_xlim(0,4)
ax.set_ylim(1,1E5)
ax.set_yscale('log')
ax.set_xscale('linear')

#Popesso groups/poor clusters
#ax.plot(zspace,213*zspace**(1.33),'m-')
pop_group_top = np.concatenate([(213+44)*zspace[:37]**(1.33-0.34),(213+44)*zspace[37:]**(1.33+0.34)])
pop_group_bottom = np.concatenate([(213-44)*zspace[:37]**(1.33+0.34),(213-44)*zspace[37:]**(1.33-0.34)])
#fill_between(zspace,pop_group_top,pop_group_bottom,facecolor='magenta',alpha=0.2,label='Popesso groups')

#Popesso clusters
#ax.plot(zspace,66*zspace**(1.77),'k-')
last_popesso_index = np.where(abs(zspace-1.6)<0.02)[0][0] #z ~ 1.6
pop_cluster_top = np.concatenate([(66+23)*zspace[:37]**(1.77-0.36),(66+23)*zspace[37:]**(1.77+0.36)])
pop_cluster_bottom = np.concatenate([(66-23)*zspace[:37]**(1.77+0.36),(66-23)*zspace[37:]**(1.77-0.36)])
#fill_between(zspace,pop_cluster_top,pop_cluster_bottom,facecolor='black',edgecolor='black',alpha=0.2,label='Popesso clusters')
fill_between(zspace[:last_popesso_index],pop_cluster_top[:last_popesso_index],pop_cluster_bottom[:last_popesso_index],facecolor='black',edgecolor='black',alpha=0.3,label='Popesso clusters')

fill_between(zspace[last_popesso_index:],pop_cluster_top[last_popesso_index:],pop_cluster_bottom[last_popesso_index:],facecolor='black',edgecolor='black',alpha=0.1)
fill_between(zspace[last_popesso_index:],pop_cluster_top[last_popesso_index:],pop_cluster_bottom[last_popesso_index:],facecolor='None',alpha=1, hatch='x')



#Red-Sequence Cluster Survey-1 (Webb+2013)
ax.scatter([0.37,0.56,0.725,0.915],[13,25,60,80],s=60,marker='o',facecolors='none',edgecolors='k',label='RCS') #clusters
ax.errorbar([0.37,0.56,0.725,0.915],[13,25,60,80],yerr = [3,7,35,65],c='k',ls='')

#GEMINI CLUSTER ASTROPHYSICS SPECTROSCOPIC SURVEY clusters (Muzzin+2012)
ax.plot([1.025,1.025],[88,215],'k.-')
ax.scatter([1.025,1.025],[88,215],s=60,marker='s',facecolors='none',edgecolors='k',label='GCLASS') #clusters

#XCSJ2215 (Ma+2015)
#the majority of the SMGs discovered in XCS J2215 lie within the core region of the cluster interspersed with the passive cluster population
ax.scatter(1.46,460,s=60,marker='X',facecolors='none',edgecolors='k')
ax.errorbar(1.46,460,yerr = [[210],[150]],c='k',ls='')
#XCSJ2215 (Ma+2015) boosted with 24um sources from Hilton+2010, and [OII] emitters from Hayashi+2010
ax.scatter(1.46,950,s=60,marker='X',facecolors='none',edgecolors='k',label='XCS2215')
ax.errorbar(1.46,950,yerr = 320,c='k',ls='')


#ClJ0218 (Smail+2014)
#cluster core dominated by passive galaxies, with the SMGs lying on the outskirts.
ax.scatter(1.6,800,s=60,marker='^',facecolors='none',edgecolors='k',label='ClJ0218')
ax.errorbar(1.6,800,yerr = 500,c='k',ls='')

#HDF structure (Chapman+2009)
ax.scatter(1.99,3400,s=60,marker='D',facecolors='none',edgecolors='k',label='HDF1.99') #structure
ax.errorbar(1.99,3400,yerr = 2100,c='k',ls='')

#Cl J1449+0856 (Strazzullo+2018)
ax.scatter(2.0,1300,s=60,marker='>',facecolors='none',edgecolors='k',label='ClJ1449',zorder=10)
ax.errorbar(2.0,1300,yerr = 400,c='k',ls='',zorder=10)

#MRC1138 (Dannerbauer+2014)
ax.scatter(2.16,1030,s=100,marker='p',facecolors='none',edgecolors='k',label='MRC1138') #protocluster
ax.errorbar(2.16,1030,yerr = 390,c='k',ls='')



#H1700
'''24um/R-band'''
#####ax.scatter(2.3,SFR_1700/1.0,s=250,marker=(5,1),facecolors='none',edgecolors='blue')
#####ax.scatter(2.3,SFR_1700_boost/1.0,s=250,marker=(5,1),facecolors='blue',edgecolors='blue',label='HS1700 (this work)')
#####ax.errorbar(2.3,SFR_1700_boost/1.0,yerr = dSFR_1700_m_boost,c='b',ls='')

##ax.errorbar(2.3,SFR_1700/1.0,yerr = dSFR_1700_m,c='b',ls='')

#SMGs only
ax.scatter(2.3,boost_1700/1.0,s=250,marker=(5,1),facecolors='blue',edgecolors='blue',label='HS1700 (this work)')
ax.errorbar(2.3,boost_1700/1.0,yerr = dboost_1700,c='b',ls='')

#SMGs + LBGs+NB-emitters - full conversion
ax.scatter(2.3,SFR_1700_boost/1.0,s=250,marker=(5,1),facecolors='none',edgecolors='blue')
ax.errorbar(2.3,SFR_1700_boost/1.0,yerr = dSFR_1700_m_boost,c='b',ls='')

'''S2'''
#ax.scatter(2.3,SFR_1700_S2/1.0,s=200,marker=(5,1),facecolors='none',edgecolors='blue')
#ax.scatter(2.3,SFR_1700_boost_S2/1.0,s=200,marker=(5,1),facecolors='blue',edgecolors='blue',label='HS1700')
#ax.errorbar(2.3,SFR_1700/1.0,yerr = dSFR_1700_m,c='b',ls='')
#ax.errorbar(2.3,SFR_1700_boost_S2/1.0,yerr = dSFR_1700_m_boost_S2,c='b',ls='')

#ax.scatter(2.3,SFR_1700_core/1.0,s=250,marker=(5,1),facecolors='none',edgecolors='black') #estimated from mean from 1549's core
#ax.scatter(2.3,(SFR_1700_core + (6.2+6.3+4.6)*200)/1.0,s=250,marker=(5,1),facecolors='black',edgecolors='black') #estimated from mean from 1549's core



#PCL1002+0222 (Casey+2015)
ax.scatter(2.47,6400,s=60,marker='v',facecolors='none',edgecolors='k',label='PCL1002+0222')
ax.errorbar(2.47,6400,yerr = 2500,c='k',ls='')



#H1549
'''24um/R-band'''
#####ax.scatter(2.85,SFR_1549/(1),s=250,marker=(5,1),facecolors='none',edgecolors='red')
#####ax.scatter(2.85,SFR_1549_boost/(1),s=250,marker=(5,1),facecolors='red',edgecolors='red',label='HS1549 (this work)')
#####ax.errorbar(2.85,SFR_1549_boost/(1),yerr = dSFR_1549_m_boost*1.5,c='r',ls='')
##ax.errorbar(2.85,SFR_1549/(1.5*1),yerr = dSFR_1549,c='r',ls='')
##ax.errorbar(2.85,SFR_1549/(1.5*1),yerr = dSFR_1549_m,c='r',ls='')

#SMGs only
ax.scatter(2.85,boost_1549/(1),s=250,marker=(5,1),facecolors='red',edgecolors='red',label='HS1549 (this work)')
ax.errorbar(2.85,boost_1549/(1),yerr = dboost_1549,c='r',ls='')

#SMGs + LBGs+NB-emitters/2 - conservative: accounting for AGN+K-correction
ax.scatter(2.85,SFR_1549_boost_conservative/1.0,s=250,marker=(5,1),facecolors='none',edgecolors='red',hatch='////')
#ax.errorbar(2.85,SFR_1549_boost_conservative/1.0,yerr = dSFR_1549_m_boost_conservative,c='r',ls='')

#SMGs + LBGs+NB-emitters - full conversion
ax.scatter(2.85,SFR_1549_boost/(1),s=250,marker=(5,1),facecolors='none',edgecolors='red')
ax.errorbar(2.85,SFR_1549_boost/(1),yerr = dSFR_1549_m_boost*1.5,c='r',ls='')


'''S2'''
#ax.scatter(2.85,SFR_1549_S2/(1.5*1),s=200,marker=(5,1),facecolors='none',edgecolors='red')
#ax.scatter(2.85,SFR_1549_boost_S2/(1.5*1),s=200,marker=(5,1),facecolors='red',edgecolors='red',label='HS1549')
#ax.errorbar(2.85,SFR_1549/(1.5*1),yerr = dSFR_1549,c='r',ls='')
#ax.errorbar(2.85,SFR_1549/(1.5*1),yerr = dSFR_1549_m,c='r',ls='')
#ax.errorbar(2.85,SFR_1549_boost_S2/(1.5*1),yerr = dSFR_1549_m_boost_S2,c='r',ls='')

# ax.scatter(2.85,2670./(1.5*1),s=250,marker=(5,1),facecolors='none',edgecolors='black') #core
# ax.scatter(2.85,(4760.+2670.)/(1.5*1),s=250,marker=(5,1),facecolors='black',edgecolors='black') #core




#SA22 protocluster (Steidel+1998)
# ax.scatter(3.09,1600,s=60,marker='d',facecolors='none',edgecolors='k',label='SA22') #protocluster
# ax.errorbar(3.09,1600,yerr = 1300,c='k',ls='')
ax.scatter(3.09,2600,s=60,marker='d',facecolors='none',edgecolors='k',label='SA22')
ax.errorbar(3.09,2600,yerr = 1300*1.5,c='k',ls='')

'''
#Oteo+2018 https://arxiv.org/pdf/1709.02809.pdf
#SFR_tot = 14,400 Msun/yr
#M_cl = 4.4 x 10^13 Msun
#z=4.002
ax.scatter(4.002,14400./0.15,s=60,marker='<',facecolors='none',edgecolors='k',label='DRC')
ax.errorbar(4.002,14400./0.15, yerr=0.6*14400./0.15, c='k',ls='')

#SPT2349
ax.scatter(4.3040,16500./0.15,s=60,marker='>',facecolors='none',edgecolors='k',label='SPT2349')
ax.errorbar(4.3040,16500./0.15, yerr=np.sqrt((0.7/1.16)**2 + 0.28**2)*16500./0.15, c='k',ls='')



#AzTEC-3 overdensity (Capak+2011)
ax.scatter(5.03,4E5,s=60,marker='H',facecolors='none',edgecolors='k',label='AzTEC-3 overdensity')
ax.errorbar(5.03,4E5,yerr = 3E5,c='k',ls='')
'''



'''
#Popesso cluster sample (Popesso+2012)
pl.plot([0.175,0.226,0.228,0.329,0.351,0.539,0.780,0.830,0.838],[],'ko',label='Popesso clusters (Popesso+2012)')
pl.errorbar([0.175,0.226,0.228,0.329,0.351,0.539,0.780,0.830,0.838],[],yerr = [[],[]],c='k',ls='o')

#Popesso groups/poor cluster sample (Popesso+2012)
pl.plot([0.1,0.22,0.35,0.395,0.70,0.735,0.85,1.02,1.61],[],'ko',label='Popesso groups/poor clusters (Popesso+2012)')
pl.errorbar([0.1,0.22,0.35,0.395,0.70,0.735,0.85,1.02,1.61],[],yerr = [[],[]],c='m',ls='o')
'''


'''
#lookback times generated from http://www.astro.ucla.edu/%7Ewright/DlttCalc.html
ax2 = ax.twiny()
ax2.set_xticks([0,0.075,0.159,0.253,0.362,0.487,0.636,0.818,1.046,1.347,1.772,2.434,3.678])
ax2.set_xticklabels(['0','1','2','3','4','5','6','7','8','9','10','11','12'])
ax2.set_xlabel('Lookback Time (Gyr)')
xlim2  = ax2.set_xlim(0,3.5)
if xlim2 != xlim:
    print "inconsitant x-limits"
'''



# ax.scatter(2.27,11632./1,s=60,marker='d',facecolors='none',edgecolors='k',label='Bootes')
# ax.scatter(0.67,620./10,s=60,marker='d',facecolors='none',edgecolors='k',label='EGS')
# ax.scatter(2.05,4924./1,s=60,marker='d',facecolors='none',edgecolors='k',label='Lockman')
# ax.scatter(1.04,1631./10,s=60,marker='d',facecolors='none',edgecolors='k',label='CDF-S')

# # #ax.errorbar(3.09,2600,yerr = 1300*1.5,c='k',ls='')
# # low-z 6e14 - 1e15
# # panck clumps 1e14

zspace = np.linspace(0,6,100)
zall = [2.3,2.85,0.37,0.56,0.725,0.915,1.025,1.6,1.46,1.99,2.47,2.0]#, 4.002, 4.3040]
sfrmall = [SFR_1700_boost/1,SFR_1549_boost/(1),13,25,60,80,215,800,950,3400,6400,1300]#, 14400./0.15, 16500./0.15]
#ax.plot(zall,sfrmall,'ko')

#pl.text(2.9,4e4,r'$(1+z)^7$')
#pl.text(3.15,1.2e4,r'$\propto(1+z)^7$')
pl.text(2.65,4e4,r'$\propto(1+z)^7$')

def power(x,b):
    return (1+x)**b
params,cov = cf(power,zall,sfrmall)
b = params[0]
ax.plot(zspace,power(zspace,b),'k')

#ax.text(2.98,10000,r'$\propto(1+z)^7$',size=14)


#ax.plot(zspace,213*zspace**(1.33),'m-')
ax.set_ylabel(r'$ SFR_{tot}/M_{cl}$ $(M_{\odot} yr^{-1} / 10^{14}M_{\odot})$')
ax.set_xlabel('z')

### fill_between and pl.fill douplicate legend label, so do this to only get one label ###
legend = ax.legend(fontsize=9,loc=4,ncol=2,frameon=True,numpoints=1,scatterpoints=1)
handles, labels = ax.get_legend_handles_labels()
from collections import OrderedDict
by_label = OrderedDict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), fontsize=9,loc=4,ncol=2,frameon=True,numpoints=1,scatterpoints=1)


pl.savefig('../Figures/SFR.pdf',bbox_inches='tight')
'''pl.savefig('test.pdf',bbox_inches='tight')'''
#pl.show()
pl.close()




'''
fig = pl.figure()
ax = pl.gca()
fig.subplots_adjust(hspace=0,wspace=0)
pl.rc('font',size=16)
pl.rc('mathtext', default='regular')
ax.set_xlim(0,10)
ax.set_ylim(1E2,1E5)
ax.set_yscale('log')
ax.set_xscale('linear')

#H1700
ax.scatter(1.0,SFR_1700,s=200,marker=(5,1),facecolors='none',edgecolors='blue')
ax.scatter(1.0,SFR_1700_boost,s=200,marker=(5,1),facecolors='blue',edgecolors='blue',label='HS1700')
ax.errorbar(1.0,SFR_1700_boost,xerr = 0.5*6,yerr = dSFR_1700_boost,c='b',ls='')

#H1549
ax.scatter(1.5*1,SFR_1549,s=200,marker=(5,1),facecolors='none',edgecolors='red')
ax.scatter(1.5*1,SFR_1549_boost,s=200,marker=(5,1),facecolors='red',edgecolors='red',label='HS1549')
ax.errorbar(1.5*1,SFR_1549_boost,xerr=0.5*1.5*6,yerr = dSFR_1549_boost,c='r',ls='')


#MRC1138 (Dannerbauer+2014)
ax.scatter(6,6188,s=60,marker='o',facecolors='c',edgecolors='c',label='MRC1138 protocluster')
ax.errorbar(6,6188,xerr = 4,yerr=6188.0/np.sqrt(58),c='c',ls='')

#SA22 protocluster (Steidel+1998)
ax.scatter(0.8,1300,s=60,marker='s',facecolors='g',edgecolors='g',label='SA22 protocluster')
ax.errorbar(0.8,1300,xerr = 0.4,yerr=800,c='g',ls='')

#HDF protocluster (Chapman+2009)
ax.scatter(0.6765,2300,s=60,marker='s',facecolors='k',edgecolors='k',label='HDF protocluster')
ax.errorbar(0.6765,2300,xerr=0.3622,yerr = 300,c='k',ls='')

#PCL1002+0222 (Casey+2015)
ax.scatter(0.8,5100,s=60,marker='s',facecolors='b',edgecolors='b',label='PCL1002+0222')
ax.errorbar(0.8,5100,xerr=0.3,yerr = 500,c='b',ls='')

#
#mspace = np.linspace(0,15,100)
#mall = np.array([0.6765,0.8,0.8,6,1,1.5*1])
#sfrall = np.array([2300,1300,5100,6188,SFR_1700_boost,SFR_1549_boost])
#
#def power(x,a,b):
#    return a*x**b
#params,cov = cf(power,mall,sfrall)
#a = params[0]
#b = params[1]
#
#ax.plot(mspace,power(mspace,a,b),'k')
#


pl.xlabel(r'$M_{cl} /10^{14}M_{\odot}$')
pl.ylabel(r'$\sum SFR (M_{\odot} yr^{-1})$')
ax.legend(fontsize=11,loc=4,ncol=2,frameon=False,numpoints=1,scatterpoints=1)
pl.savefig('../Figures/SFR_Mcl.pdf',bbox_inches='tight')
#pl.show()
pl.close()
'''
