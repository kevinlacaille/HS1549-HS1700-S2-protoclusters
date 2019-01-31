import aplpy as ap
import numpy as np
import pylab as pl
import coords as co

#open file
cat_850 = open('../h1700/flux_4sig.gaia','r')
cat_MIPS = open('../h1700/MIPS/MIPS_ID.gaia','r')
cat_op = open('../../Protocluster_paper/H1700/Data/Catalogues_new/q1700_all_radec_python.gaia','r')
cat_IRAC = open('../h1700/IRAC/IRAC_ID.gaia','r')

#sources without multi-lambda data
data_24_b = ['*','*','*',4,'*','*','*','*',9,'*','*','*','*','*',15,'*','*',18,19,'*','*',22,'*','*','*',26,'*']
data_ch4_b = ['*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*']
data_ch2_b = ['*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*','*']
data_J_b = ['*','*','*',4,'*','*','*','*',9,10,'*','*','*','*',15,'*','*','*',19,'*','*',22,'*',24,25,26,'*']

#import 850um data
ID = []
RA = []
DEC = []
flux_850 = []
for line in cat_850.readlines()[2:]:
    tmp = line.split()
    ID.append(tmp[0])
    RA.append(co.convHMS(tmp[1]))
    DEC.append(co.convDMS(tmp[2]))
    flux_850.append(float(tmp[3]))

#import MIPS IDs
ID_MIPS = []
RA_MIPS = []
DEC_MIPS = []
for line in cat_MIPS.readlines()[2:]:
    tmp = line.split()
    if tmp[1][0:2] == '17':
        ID_MIPS.append(tmp[0])
        RA_MIPS.append(co.convHMS(tmp[1]))
        DEC_MIPS.append(co.convDMS(tmp[2]))
    else:
        pass

#import IRAC IDs
ID_IRAC = []
RA_IRAC = []
DEC_IRAC = []
for line in cat_IRAC.readlines()[2:]:
    tmp = line.split()
    if tmp[1][0:2] == '17':
        ID_IRAC.append(tmp[0])
        RA_IRAC.append(co.convHMS(tmp[1]))
        DEC_IRAC.append(co.convDMS(tmp[2]))
    else:
        pass

#import "optical" data
ID_op = []
RA_op = []
DEC_op = []
z = []
for line in cat_op.readlines()[2:]:
    tmp = line.split()
    ID_op.append(tmp[0])
    RA_op.append(co.convHMS(tmp[1]))
    DEC_op.append(co.convDMS(tmp[2]))
    z.append(float(tmp[4]))
RA_op = np.array(RA_op)
DEC_op = np.array(DEC_op)
z = np.array(z)

#set max flux so we can have correct contrast
max_850 = flux_850 #12
max_450 = [27.5471,16.5369,12.3198,26.8498,19.3488,8.85379,20.3172,13.5599,32.1374,19.5339,27.4439,13.2162,14.433,17.3893,18.4549,27.5171,14.9662,17.7167,20.5034,14.0107,14.6879,20.5114,12.858,19.8971,20.2303,16.603,13.221]
max_24 = [0.126542,0.170227,0.198449,0,0.234537,0.0514561,0.413828,0.623477,0,0.335913,0.600445,12.4331,0.0423905,0.333705,0,1.16658,0.074127,0,0,0.566733,0.0636076,0,0.264334,0.0667247,0.127113,0,0.184137]
max_ch2 = [0.00321828,0.00327927,0.00755771,0.001,0.00473219,0.00415315,0.00220657,0.00584779,0.005,0.00747939,0.01,0.55387,0.0005,0.00722975,0.00514863,0.0567032,0.0017,0.00987087,0.00231791,0.0264548,0.00075,0.00174437,0.00370191,0.0042353,0.0100492,0.0108639,0.00154]
max_J = [3453.04,3453.04,4000.0,0,5000.0,4500.0,7000.0,6000.0,0,0,6000.0,10000.0,2600.0,2500.0,0,15000.0,4500.0,5000.0,0,18000.0,2500,0,2200.0,0,0,0,2500.0]#13

#data
loc = '../../Protocluster_paper/H1700/Data/'
h_850 = '../h1700/h1700_850.fits'
h_450 = '../h1700/h1700_450.fits'
h_850_snr = '../h1700/h1700_850_snr.fits'
h_450_snr = '../h1700/h1700_450_snr.fits'
h_24 = loc + 'q1700_24um.fits'
h_ch4 = loc + 'h1700_ch4_Rs_new.fits'
h_ch3 = loc + 'h1700_ch3_Rs_new.fits'
h_ch2 = loc + 'h1700_ch2_Rs_fixed.fits'
h_ch1 = loc + 'h1700_ch1_Rs_new.fits'
h_J = loc + 'q1700Jreg.fits'



#whole field - 850um
fig = pl.figure(figsize=(8,8),frameon=False)
pl.rc('font',size=22)
pl.rc('mathtext', default='regular')

RA_c = co.convHMS('17:01:06.362') - 1.0/60
DEC_c = co.convDMS('64:12:44.19') - 0.25/60

gc = ap.FITSFigure(h_850_snr, figure = fig)
#gc.tick_labels.set_xformat('hh:mm:ss')
#gc.tick_labels.set_yformat('dd:mm:ss')
gc.hide_xaxis_label()
gc.hide_xtick_labels()
gc.hide_yaxis_label()
gc.hide_ytick_labels()
gc.recenter(RA_c, DEC_c,width = 15./60, height= 15./60)
gc.show_colorscale(cmap='gray_r',vmax=4,vmin=-2,interpolation='nearest')
#gc.add_colorbar(axis_label_text=r'S$_{850\mu m}$ (mJy)',axis_label_rotation=-90,axis_label_pad=10.5)

for i in range(len(RA)):
    if ID[i] == '3' or ID[i] == '5' or ID[i] == '8' or ID[i] == '14':
        gc.show_circles(RA[i], DEC[i],radius = 15.0/3600,color='red',linewidth=3)
    elif ID[i] == '1' or ID[i] == '12' or ID[i] == '15' or ID[i] == '16':
        gc.show_circles(RA[i], DEC[i],radius = 15.0/3600,color='blue',linewidth=3)
    else:
        gc.show_circles(RA[i], DEC[i],radius = 15.0/3600,color='green',linewidth=3)
    if ID[i] != '19':
        gc.add_label(RA[i],DEC[i]+25.0/3600,ID[i],color = 'black',size=18)
    else:
        gc.add_label(RA[i],DEC[i]+25.0/3600,ID[i],color = 'black',size=18)

for i in range(len(RA_op)):
    if 2.25 <= z[i] <= 2.35:
        gc.show_markers(RA_op[i], DEC_op[i],marker='.', facecolor='darkorange', edgecolor='darkorange',s=80)

gc.show_circles(np.mean(RA_op[np.where(abs(z-2.3)<=.05)[0]]),np.mean(DEC_op[np.where(abs(z-2.3)<=.05)[0]]),radius=1.5/60,linewidth=3,color='black',linestyle='--')
#gc.show_circles(np.median(RA_op[np.where(abs(z-2.3)<=.05)[0]]),np.median(DEC_op[np.where(abs(z-2.3)<=.05)[0]]),radius=1.5/60,linewidth=3,color='cyan')
#
#print 'offset = ' + str(3600*(np.median(DEC_op[np.where(abs(z-2.3)<=.05)[0]]) - np.mean(DEC_op[np.where(abs(z-2.3)<=.05)[0]])))

#gc.show_contour(h_850,hdu=1,levels = [4], colors = 'limegreen',linewidths=3)
#gc.show_contour(h_ch2,levels = [0.001], colors = 'red',linewidths=2)
#gc.show_contour(h_24,levels = [0.0001], colors = 'blue',linewidths=2)

gc.add_scalebar(2.0/60,linewidth=3, corner='bottom left') #2'
gc.scalebar.set_label("2 arcmin")
gc.scalebar.set_color('black')
gc.save('../Figures/h1700_850_map_gray.pdf',dpi=250)
#pl.show()
pl.close()
gc.close()



#core
fig = pl.figure(figsize=(8,8),frameon=False)
pl.rc('font',size=22)
pl.rc('mathtext', default='regular')

RA_c = np.median(RA_op)#co.convHMS('17:01:01.00')
DEC_c = np.median(DEC_op)#co.convDMS('64:12:38.00')

gc = ap.FITSFigure(h_ch2, figure = fig)
#gc.tick_labels.set_xformat('hh:mm:ss')
#gc.tick_labels.set_yformat('dd:mm:ss')
gc.hide_xaxis_label()
gc.hide_xtick_labels()
gc.hide_yaxis_label()
gc.hide_ytick_labels()
gc.recenter(RA_c, DEC_c,width = 4.0/60, height= 4.0/60)
gc.show_colorscale(cmap='gray_r',vmax=0.0076/2,vmin=-0.00025*3,interpolation='nearest')
#gc.add_colorbar(axis_label_text=r'S$_{850\mu m}$ (mJy)',axis_label_rotation=-90,axis_label_pad=10.5)

#gc.show_circles(co.convHMS('17:01:00.6'),co.convDMS('64:12:06.4'),radius = 15.0/3600,color='red')
gc.show_contour(h_850_snr,levels = [4,7,11], colors = 'limegreen',linewidths=2)
gc.show_contour(h_450_snr,levels = [3], colors = 'magenta',linewidths=2)

for j in range(len(ID_IRAC)):
    if ID_IRAC[j][5:].split('_')[0] == '3' or ID_IRAC[j][5:].split('_')[0] == '5' or ID_IRAC[j][5:].split('_')[0] == '8' or ID_IRAC[j][5:].split('_')[0] == '14' or ID_IRAC[j][5:].split('_')[0] == '20':
        gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 2.5/3600,color='green',linewidth=3)
    if ID_IRAC[j][5:] == '8_1':
        gc.add_label(RA_IRAC[j]+5.0/3600,DEC_IRAC[j]+15./3600,ID_IRAC[j][5:],color = 'black',size=22)
#    if ID_IRAC[j][5:] == '8_2':
#        gc.add_label(RA_IRAC[j]-20.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'white',size=12)
    if ID_IRAC[j][5:] == '3_1':
            gc.add_label(RA_IRAC[j]-10./3600,DEC_IRAC[j]+15.0/3600,ID_IRAC[j][5:],color = 'black',size=22)
#    if ID_IRAC[j][5:] == '3_2':
#        gc.add_label(RA_IRAC[j]+20.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'white',size=12)
#    if ID_IRAC[j][5:] == '3_3':
#        gc.add_label(RA_IRAC[j]-20.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'white',size=12)
    if ID_IRAC[j][5:] == '5_2':
        gc.add_label(RA_IRAC[j]+40.0/3600,DEC_IRAC[j]+10./3600,ID_IRAC[j][5:],color = 'black',size=22)
#    if ID_IRAC[j][5:] == '5_1':
#        gc.add_label(RA_IRAC[j],DEC_IRAC[j]-10.0/3600,ID_IRAC[j][5:],color = 'white',size=12)
#    if ID_IRAC[j][5:] == '5_3':
#        gc.add_label(RA_IRAC[j],DEC_IRAC[j]+10.0/3600,ID_IRAC[j][5:],color = 'white',size=12)
    if ID_IRAC[j][5:] == '14':
        gc.add_label(RA_IRAC[j]+35.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'black',size=22)
#    if ID_IRAC[j][5:] == '12':
#        gc.add_label(RA_IRAC[j],DEC_IRAC[j]+10.0/3600,ID_IRAC[j][5:],color = 'white',size=12)
#    if ID_IRAC[j][5:] == '6_1':
#        gc.add_label(RA_IRAC[j]-20.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'white',size=12)
#    if ID_IRAC[j][5:] == '6_2':
#        gc.add_label(RA_IRAC[j]+20.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'white',size=12)
#    if ID_IRAC[j][5:] == '2_1':
#        gc.add_label(RA_IRAC[j]+20.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'white',size=12)
#    if ID_IRAC[j][5:] == '2_2':
#        gc.add_label(RA_IRAC[j]-20.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'white',size=12)

    if ID_IRAC[j][5:] == '3_1' or ID_IRAC[j][5:] == '5_2' or ID_IRAC[j][5:] == '8_1' or ID_IRAC[j][5:] == '14':
        gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 2.5/3600,color='red',linewidth=3)
    if ID_IRAC[j][5:] == '16' or ID_IRAC[j][5:] == '1' or ID_IRAC[j][5:] == '12' or ID_IRAC[j][5:].split('_')[0] == '15':
        gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 2.5/3600,color='blue',linewidth=3)
    if ID_IRAC[j][5:].split('_')[0] == '2' or ID_IRAC[j][5:].split('_')[0] == '6':
        gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 2.5/3600,color='green',linewidth=3)


gc.show_circles(RA_c,DEC_c,radius = 1.5/60,color='black',linewidth=5,linestyle='--')
gc.add_label(RA_c,DEC_c+1.57/60+5.0/3600,'3-arcmin core',color = 'black',size=22)

gc.add_scalebar(0.5/60,linewidth=3)
gc.scalebar.set_label("0.5 arcmin")
gc.scalebar.set_color('black')
gc.save('../Figures/h1700_ch2_core.pdf',dpi=250)
#pl.show()
pl.close()
gc.close()



#whole field - 450um
fig = pl.figure(figsize=(8,8),frameon=False)
pl.rc('font',size=15)
pl.rc('mathtext', default='regular')

RA_c = co.convHMS('17:01:06.362')
DEC_c = co.convDMS('64:12:44.19')

gc = ap.FITSFigure(h_450, figure = fig)
gc.tick_labels.set_xformat('hh:mm:ss')
gc.tick_labels.set_yformat('dd:mm:ss')
gc.recenter(RA_c, DEC_c,width = 12.0/60, height= 12.0/60)
gc.show_colorscale(cmap='gray_r',vmax=22,vmin=-10,interpolation='nearest')
gc.add_colorbar(axis_label_text=r'S$_{850\mu m}$ (mJy)',axis_label_rotation=-90,axis_label_pad=10.5)

for i in range(len(RA)):
    gc.show_circles(RA[i], DEC[i],radius = 15.0/3600,color='blue')
    gc.add_label(RA[i],DEC[i]+25.0/3600,ID[i],color = 'white',size=15)

gc.add_scalebar(2.0/60) #2'
gc.scalebar.set_label("2 arcmin")
gc.scalebar.set_color('white')
gc.save('../Figures/h1700_450_map.pdf',dpi=250)
pl.show()
pl.close()
gc.close()



p = []
p_ID = []
p_RA = []
p_DEC = []
p_theta = []

#850,450,24,ch4,ch2,J-band
for i in range(len(ID)+1)[1:]:
    # print i
    # i=1
    pl.figure()
    fig = pl.figure()
    pl.rc('font',size=25)

    RA_c = RA[i-1]
    DEC_c = DEC[i-1]

    #850
    gc = ap.FITSFigure(h_850, figure=fig, subplot=[0,0,0.75,1])
    gc.hide_xaxis_label()
    gc.hide_xtick_labels()
    gc.hide_yaxis_label()
    gc.hide_ytick_labels()
    gc.tick_labels.set_yformat('dd:mm:ss')
    gc.show_colorscale(cmap='gray_r',vmin = -3,vmax = max_850[i-1])
    gc.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
    gc.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
    for j in range(len(RA_op)):
       if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

           if 2.23 <= z[j] <= 2.37:
               gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='red',linewidth=2) #1''
               #gc.add_label(RA_op[j],DEC_op[j]+4E-4,ID_op[j],color = 'red')

           elif z[j] < 2.37:
               gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''
               #gc.add_label(RA_op[j],DEC_op[j]+4E-4,ID_op[j],color = 'cyan')

           if z[j] > 2.37:
               gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''
               #gc.add_label(RA_op[j],DEC_op[j]+4E-4,ID_op[j],color = 'cyan')

    gc.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
    for j in range(len(ID_IRAC)):
       if int(ID_IRAC[j][5:].split('_')[0]) == i:
           gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 1.0/3600,color='lime',linewidth=2)

    gc.add_scalebar(0.0)
    gc.scalebar.set_label(r'850 $\mu$m')
    gc.scalebar.set_color(color='black')

    #450
    gc = ap.FITSFigure(h_450, figure=fig, subplot=[0.75,0,0.75,1])
    gc.hide_yaxis_label()
    gc.hide_ytick_labels()
    gc.hide_xaxis_label()
    gc.hide_xtick_labels()
    gc.show_colorscale(cmap='gray_r',vmin = -10,vmax = max_450[i-1])
    gc.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
    gc.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
    gc.add_scalebar(0.0)
    gc.scalebar.set_label(r'450 $\mu$m')
    gc.scalebar.set_color(color='black')
    for j in range(len(RA_op)):
       if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

           if 2.23 <= z[j] <= 2.37:
               gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='red',linewidth=2) #1''
           elif z[j] < 2.37:
               gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''
           if z[j] > 2.37:
               gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''


    gc.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
    for j in range(len(ID_IRAC)):
       if int(ID_IRAC[j][5:].split('_')[0]) == i:
           gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 1.0/3600,color='lime',linewidth=2)

    #24um
    if data_24_b[i-1]!=i:
        gc = ap.FITSFigure(h_24, figure=fig, subplot=[1.5,0,0.75,1])
        gc.hide_yaxis_label()
        gc.hide_ytick_labels()
        gc.hide_xaxis_label()
        gc.hide_xtick_labels()
        gc.show_colorscale(cmap='gray_r',vmin = -0.03,vmax = max_24[i-1])
        gc.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
        gc.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
        gc.add_scalebar(0.0)
        gc.scalebar.set_label(r'24 $\mu$m')
        gc.scalebar.set_color(color='black')

        for j in range(len(ID_IRAC)):
            if int(ID_IRAC[j][5:].split('_')[0]) == i:
    #            gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 1.0/3600,color='lime')

                theta = 3600*np.sqrt(((RA_c-RA_IRAC[j])*np.cos(64*np.pi/180.0))**2 + (DEC_c-DEC_IRAC[j])**2)
                n = 1.25E-3 #1/(np.pi*(7.5)**2)
                pval = 1-np.exp(-np.pi*n*theta**2)
                print ID_IRAC[j] + ' ' + str(round(pval,3))
                p.append(pval)
                p_ID.append(ID_IRAC[j])
                p_RA.append(RA_IRAC[j])
                p_DEC.append(DEC_IRAC[j])
                p_theta.append(theta)

        for j in range(len(RA_op)):
           if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

               if 2.23 <= z[j] <= 2.37:
                   gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='red',linewidth=2) #1''
               elif z[j] < 2.37:
                   gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''
               if z[j] > 2.37:
                   gc.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''

        gc.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
        for j in range(len(ID_IRAC)):
          if int(ID_IRAC[j][5:].split('_')[0]) == i:
            gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 1.0/3600,color='lime',linewidth=2)


    else:
       gc = ap.FITSFigure(h_850, figure=fig, subplot=[1.5,0,0.75,1])
       gc.hide_yaxis_label()
       gc.hide_ytick_labels()
       gc.hide_xaxis_label()
       gc.hide_xtick_labels()
       gc.add_scalebar(0.0)
       gc.scalebar.set_label(r'24 $\mu$m')
       gc.scalebar.set_color(color='black')
       gc.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
       gc.add_label(RA_c,DEC_c,'no data',color='black')

    #8um
    if data_ch4_b[i-1]!=i:
       if i == 9 or i == 16:
           image_ch4 = h_ch3
       else:
           image_ch4 = h_ch4
       gc_ch4 = ap.FITSFigure(image_ch4, figure=fig, subplot=[0,-1,0.75,1])
       gc_ch4.tick_labels.set_xformat('hh:mm:ss')
       gc_ch4.tick_labels.set_yformat('dd:mm:ss')
       gc_ch4.hide_xaxis_label()
       gc_ch4.hide_xtick_labels()
       gc_ch4.hide_yaxis_label()
       gc_ch4.hide_ytick_labels()
       for j in range(len(ID_IRAC)):
           if int(ID_IRAC[j][5:].split('_')[0]) == i:
               theta = 3600*np.sqrt(((RA_c-RA_IRAC[j])*np.cos(64*np.pi/180.0))**2 + (DEC_c-DEC_IRAC[j])**2)
               n = 1/(np.pi*(7.5)**2)
               pval = 1-np.exp(-np.pi*n*theta**2)
               print ID_IRAC[j] + ' ' + str(pval)
               p.append(pval)
               p_ID.append(ID_IRAC[j])
               p_RA.append(RA_IRAC[j])
               p_DEC.append(DEC_IRAC[j])
               p_theta.append(theta)


       gc_ch4.show_colorscale(cmap='gray_r',vmin = -0.0005,vmax = max_ch2[i-1])
       gc_ch4.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
       gc_ch4.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
       gc_ch4.add_scalebar(0.0)
       if i == 9 or i == 16:
           gc_ch4.scalebar.set_label(r'5.6 $\mu$m')
           gc_ch4.scalebar.set_color(color='black')
       else:
           gc_ch4.scalebar.set_label(r'8.0 $\mu$m')
           gc_ch4.scalebar.set_color(color='black')
       for j in range(len(RA_op)):
           if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

               if 2.23 <= z[j] <= 2.37:
                   gc_ch4.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='red',linewidth=2) #1''
               elif z[j] < 2.37:
                   gc_ch4.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''
               if z[j] > 2.37:
                   gc_ch4.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''

       gc_ch4.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
       for j in range(len(ID_IRAC)):
           if int(ID_IRAC[j][5:].split('_')[0]) == i:
               gc_ch4.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 1.0/3600,color='lime',linewidth=2)



    else:
       gc_ch4 = ap.FITSFigure(h_850, figure=fig, subplot=[0,-1,0.75,1])
       gc_ch4.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
       gc_ch4.hide_yaxis_label()
       gc_ch4.hide_ytick_labels()
       gc_ch4.hide_xaxis_label()
       gc_ch4.hide_xtick_labels()
       gc_ch4.add_scalebar(0.0)
       gc_ch4.scalebar.set_label(r'8.0 $\mu$m')
       gc.scalebar.set_color(color='black')
       gc_ch4.add_label(RA_c,DEC_c,'no data',color='black')

    #4.5um
    if data_ch2_b[i-1]!=i:
       if i == 9 or i ==16: #and 16?
           image_ch2 = h_ch1
       else:
           image_ch2 = h_ch2
       gc_ch2 = ap.FITSFigure(image_ch2, figure=fig, subplot=[0.75,-1,0.75,1])
       gc_ch2.tick_labels.set_xformat('hh:mm:ss')
       gc_ch2.hide_xaxis_label()
       gc_ch2.hide_xtick_labels()
       gc_ch2.hide_yaxis_label()
       gc_ch2.hide_ytick_labels()
       gc_ch2.show_colorscale(cmap='gray_r',vmin =-0.00025,vmax = max_ch2[i-1])
       gc_ch2.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
       gc_ch2.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
       gc_ch2.add_scalebar(0.0)
       if i == 9 or i ==16:
           gc_ch2.scalebar.set_label(r'3.6 $\mu$m')
           gc_ch2.scalebar.set_color(color='black')
       else:
           gc_ch2.scalebar.set_label(r'4.5 $\mu$m')
           gc_ch2.scalebar.set_color(color='black')
       for j in range(len(RA_op)):
           if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

               if 2.23 <= z[j] <= 2.37:
                   gc_ch2.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='red',linewidth=2) #1''
               elif z[j] < 2.37:
                   gc_ch2.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''
               if z[j] > 2.37:
                   gc_ch2.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''

       gc_ch2.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
       for j in range(len(ID_IRAC)):
           if int(ID_IRAC[j][5:].split('_')[0]) == i:
               gc_ch2.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 1.0/3600,color='lime',linewidth=2)
               gc_ch2.add_label(RA_IRAC[j],DEC_IRAC[j]+4E-4,ID_IRAC[j],color = 'lime')

    else:
       gc_ch2 = ap.FITSFigure(h_850, figure=fig, subplot=[0.75,-1,0.75,1])
       gc_ch2.hide_yaxis_label()
       gc_ch2.hide_ytick_labels()
       gc_ch2.hide_xaxis_label()
       gc_ch2.hide_xtick_labels()
       gc_ch2.add_scalebar(0.0)
       gc_ch2.scalebar.set_label(r'4.5 $\mu$m')
       gc_ch2.scalebar.set_color(color='black')
       gc_ch2.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
       gc_ch2.add_label(RA_c,DEC_c,'no data',color='black')


    #1.2um
    if data_J_b[i-1]!=i:
       gc_J = ap.FITSFigure(h_J, figure=fig, subplot=[1.5,-1,0.75,1])
       gc_J.hide_xaxis_label()
       gc_J.hide_xtick_labels()
       gc_J.hide_yaxis_label()
       gc_J.hide_ytick_labels()
       gc_J.tick_labels.set_xformat('hh:mm:ss')
       gc_J.show_colorscale(cmap='gray_r',vmin = 500,vmax = max_J[i-1])
       gc_J.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
       gc_J.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
       gc_J.add_scalebar(0.0)
       gc_J.scalebar.set_label(r'1.2 $\mu$m')
       gc_J.scalebar.set_color(color='black')
       for j in range(len(RA_op)):
           if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

               if 2.23 <= z[j] <= 2.37:
                   gc_J.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='red',linewidth=2) #1''
               elif z[j] < 2.37:
                   gc_J.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''
               if z[j] > 2.37:
                   gc_J.show_circles(RA_op[j], DEC_op[j],radius = 1.0/3600, color='cyan',linewidth=2) #1''

       gc_J.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
       for j in range(len(ID_IRAC)):
           if int(ID_IRAC[j][5:].split('_')[0]) == i:
               gc_J.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 1.0/3600,color='lime',linewidth=2)

    else:
       gc_J = ap.FITSFigure(h_850, figure=fig, subplot=[1.5,-1,0.75,1])
       gc_J.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
       gc_J.hide_yaxis_label()
       gc_J.hide_ytick_labels()
       gc_J.hide_xaxis_label()
       gc_J.hide_xtick_labels()
       gc_J.add_scalebar(0.0)
       gc_J.scalebar.set_label(r'1.2 $\mu$m')
       gc.scalebar.set_color(color='black')
       gc_J.add_label(RA_c,DEC_c,'no data',color='black')

#    gc.save('../Figures/cutouts/h1700/h1700_%i.png'%(i))
    gc.save('../../Protocluster_paper/Paper_HS1549_HS1700/final_September17_2018_MNRAS/Figures/Cutouts/h1700/h1700_%i.pdf'%(i),dpi=50)

    gc.close()

    if data_24_b[i-1]!=i:
       gc.close()
    if data_ch4_b[i-1]!=i:
       gc_ch4.close()
    if data_ch2_b[i-1]!=i:
       gc_ch2.close()
    if data_J_b[i-1]!=i:
       gc_J.close()


    pl.close()
#    if i > 0:
    # break


p_RA_new = []
for r in p_RA:
    r = co.deg2HMS(r)
    tmp = r.replace(r[-6:],str(round(float(r[-6:]),1)))
    p_RA_new.append(tmp.replace("17:", "$17^{\mathrm{h}}").replace(":","^{\mathrm{m}}")+"^{\mathrm{s}}$")

p_DEC_new = []
for d in p_DEC:
    d = co.deg2DMS(d)
    tmp = d.replace("64:", "$64^{\mathrm{o}}").replace(":","'")+"''$"
    p_DEC_new.append(tmp)

for i in range(len(p)):
    if len(p_ID[i].split('_')) == 2:
        if p[i]<=0.05:
            print "\\textbf{" + p_ID[i].split('_')[1] + "}" + '\t&\t' + p_RA_new[i] + '\t&\t' + p_DEC_new[i] + '\t&\t' + str(round(p_theta[i],2)) + '\t&\t' + str(round(p[i],3)) + '\\\\'
        elif 0.05<p[i]<=0.1:
            print "\\emph{" + p_ID[i].split('_')[1] + "}" + '\t&\t' + p_RA_new[i] + '\t&\t' + p_DEC_new[i] + '\t&\t' + str(round(p_theta[i],2)) + '\t&\t' + str(round(p[i],3)) + '\\\\'
        else:
            print p_ID[i].split('_')[1] + '\t&\t' + p_RA_new[i] + '\t&\t' + p_DEC_new[i] + '\t&\t' + str(round(p_theta[i],2)) + '\t&\t' + str(round(p[i],3)) + '\\\\'
    else:
        if p[i]<=0.05:
            print "\\textbf{" + p_ID[i].split('_')[1] + "}" + '\\textbf{\_' + p_ID[i].split('_')[2] + '}\t&\t' + p_RA_new[i] + '\t&\t' + p_DEC_new[i] + '\t&\t' + str(round(p_theta[i],2)) + '\t&\t' + str(round(p[i],3)) + '\\\\'
        elif 0.05<p[i]<=0.1:
            print "\\emph{" + p_ID[i].split('_')[1] + "}" + '\_' + p_ID[i].split('_')[2] + '\t&\t' + p_RA_new[i] + '\t&\t' + p_DEC_new[i] + '\t&\t' + str(round(p_theta[i],2)) + '\t&\t' + str(round(p[i],3)) + '\\\\'
        else:
            print p_ID[i].split('_')[1] + '\_' + p_ID[i].split('_')[2] + '\t&\t' + p_RA_new[i] + '\t&\t' + p_DEC_new[i] + '\t&\t' + str(round(p_theta[i],2)) + '\t&\t' + str(round(p[i],3)) + '\\\\'
count = 0
print "IN: "
for i in range(len(p)):
    if p[i]<=0.05:
        count+=1
        print p_ID[i],round(p[i],3)
print str(count) + '/' + str(len(p))

count = 0
print "TENTITIVE: "
for i in range(len(p)):
    if 0.05<p[i]<=0.1:
        count+=1
        print p_ID[i],round(p[i],3)
print str(count) + '/' + str(len(p))


count = 0
for i in range(len(ID_MIPS)):
    if p[i] <=0.2160247106210712:
        inn = 'yes'
        count +=1
    else:
        inn = 'no'
    print ID_MIPS[i] + '\t' + str(round(p[i],3)) + '\t' + inn
print "there are " + str(count) + "MIPS correlations"

pl.figure()
pl.hist(p,np.concatenate([np.linspace(0,0.3,11),[max(p)]]))
pl.axvline(x=0.05,c='r',ls='--')
pl.xlabel('p')
pl.ylabel('N/bin')
pl.rc('font',size=15)
pl.rc('mathtext', default='regular')
pl.savefig('p_hist.png')
pl.show()
