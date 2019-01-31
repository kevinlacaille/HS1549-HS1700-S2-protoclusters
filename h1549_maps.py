import aplpy as ap
import numpy as np
import matplotlib.pyplot as pl
import coords as co

#open file
cat_850 = open('../h1549/flux_4sig.gaia','r')
cat_op = open('../../Protocluster_paper/H1549/Data/q1549_pdb_lae_lbg.gaia','r')
cat_IRAC = open('../h1549/IRAC/IRAC_ID.gaia','r')

#sources without multi-lambda data
data_24_b = ['*','*','*','*',5,6,7,'*','*','*','*','*','*',14,'*','*','*',18,19,20,'*','*','*','*','*','*','*','*','*']#29???
data_ch4_b = ['*','*','*','*',5,6,7,'*','*','*','*','*','*',14,'*','*','*',18,19,20,'*','*','*','*','*','*','*','*',29]
data_ch2_b = ['*','*','*','*',5,6,7,'*','*','*','*','*','*',14,'*','*','*',18,19,20,'*','*','*','*','*','*','*','*',29]
data_J_b = ['*','*','*','*',5,6,7,8,'*','*','*',12,'*',14,'*',16,17,18,19,20,'*','*',23,'*',25,'*','*',28,29]

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

#import IRAC IDs
ID_IRAC = []
RA_IRAC = []
DEC_IRAC = []
for line in cat_IRAC.readlines()[2:]:
    tmp = line.split()
    if tmp[1][0:2] == '15':
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
    z.append(float(tmp[3]))

#set max flux so we can have correct contrast
max_850 = flux_850 #12
max_450 = [26.8383,17.9662,6.56305,16.1776,17.9415,17.2804,34.9245,28.9266,13.3674,9.85836,5.0,12.0,8.5,28.255,14.2235,11.1451,13.9354,27.5608,14.531,27.3275,10.6884,13.5813,20.2269,9.08099,11.1413,10.5075,12.9884,12.9235,9.03872]
max_24 = [0.25,0.223352,0.0643796,0.0643796,0,0,0,0.0643796,0.0484784,0.2,0.0396578,0.195471,0.109745,0,0.104586,0.770801,0.172739,0,0,0,0.158146,0.221373,0.1378,0.0441169,0.12027,0.260649,0.062402,0.0707372,0.195316]#5,6,7,14,18,19,20,29???
max_ch2 = [0.3,0.125,0.05,0.1,0,0,0,0.422203,0.0779297,0.339852,0.0868136,0.0542203,0.03,0,0.0587655,0.206,0.337176,0,0,0,0.139199,0.2,0.15,0.0398506,0.223425,0.144,0.07,0.07,0]#5,6,7,14,18,19,20,29
max_ch4 = [1.6,1.7,1.55,1.6,0,0,0,1.9,1.55,1.65,1.52,1.5,1.5,0,1.5,2.0,1.8,0,0,0,1.6,1.7,1.7,1.5,1.6,1.7,1.525,1.6,0]#5,6,7,14,18,19,20,29
max_J = [15.0,15.0,7.5,7.5,0,0,0,0,6.0,15.0,10.0,0,12.0,0,10.0,0,0,0,0,0,8,11,0,6.5,0,20.0,7.0,0,0]#5,6,7,8,12,14,16,17,18,19,20,23,25,28,29

#data
loc = '../../Protocluster_paper/H1549/Data/'
h_850 = '../h1549/h1549_850.fits'
h_450 = '../h1549/h1549_450.fits'
h_850_snr = '../h1549/h1549_850_snr.fits'
h_450_snr = '../h1549/h1549_450_snr.fits'
h_24 = loc + 'h1549_24um.fits'
h_ch4 = loc + 'h1549_irac_ch4_mosaic_Rot.fits'
h_ch2 = loc + 'h1549_irac_ch2_mosaic_Rot.fits'
h_J = loc + 'JK/q1549_J_Rot.fits'
h_SMA = loc + 'SMA/H1549COMBOmJy.fits'

'''24:4,8,12,'''


#whole field - 850um
fig = pl.figure(figsize=(8,8),frameon=False)
pl.rc('font',size=22)
pl.rc('mathtext', default='regular')

RA_c = co.convHMS('15:51:52.604')
DEC_c = co.convDMS('19:11:22.90')

gc = ap.FITSFigure(h_850_snr, figure = fig)
#gc.tick_labels.set_xformat('hh:mm:ss')
#gc.tick_labels.set_yformat('dd:mm:ss')
gc.hide_xaxis_label()
gc.hide_xtick_labels()
gc.hide_yaxis_label()
gc.hide_ytick_labels()
gc.recenter(RA_c, DEC_c,width=15./60,height=15./60)#width = 11.3/60, height= 11.3/60)
gc.show_colorscale(cmap='gray_r',vmax=4,vmin=-2,interpolation='nearest')
#gc.add_colorbar(axis_label_text=r'S$_{850\mu m}$ (mJy)',axis_label_rotation=-90,axis_label_pad=10.5)


for i in range(len(RA[1:])):
    if ID[i+1] == '2':
        gc.show_circles(RA[i+1], DEC[i+1],radius = 15.0/3600,color='red',linewidth=3)
    elif ID[i+1] == '3':
        gc.show_circles(RA[i+1], DEC[i+1],radius = 15.0/3600,color='blue',linewidth=3)
    elif ID[i+1] == '4':
        gc.show_circles(RA[i+1], DEC[i+1],radius = 15.0/3600,color='blue',linewidth=3)
    elif ID[i+1] == '10':
        gc.show_circles(RA[i+1], DEC[i+1],radius = 15.0/3600,color='blue',linewidth=3)
    else:
        gc.show_circles(RA[i+1], DEC[i+1],radius = 15.0/3600,color='green',linewidth=3)
    if i != 22:
        gc.add_label(RA[i+1],DEC[i+1]+25.0/3600,ID[i+1],color = 'black',size=18)
    else:
        gc.add_label(RA[i+1],DEC[i+1]+27.0/3600,ID[i+1],color = 'black',size=18)

gc.show_circles(co.convHMS('15:51:53.8'), co.convDMS('+19:11:09.860'),radius = 4.5/3600,color='red',linewidth=3)
#gc.add_label(co.convHMS('15:51:53.8')+18.0/3600, co.convDMS('+19:11:09.860')+2.0/3600,'1_1',color='black',size=18)

gc.show_circles(co.convHMS('15:51:53.2'), co.convDMS('+19:10:59.100'),radius = 4.5/3600,color='red',linewidth=3)
gc.add_label(co.convHMS('15:51:53.2')+0./3600, co.convDMS('+19:10:59.100')+27.0/3600,'1_(1,2,3)',color='black',size=18)

gc.show_circles(co.convHMS('15:51:52.5'), co.convDMS('+19:11:03.860'),radius = 4.5/3600,color='red',linewidth=3)
#gc.add_label(co.convHMS('15:51:52.5')-1.0/3600, co.convDMS('+19:11:03.860')-15.0/3600,'1_3',color='black',size=18)

for i in range(len(RA_op)):
    if 2.82 <= z[i] <= 2.88:
        gc.show_markers(RA_op[i], DEC_op[i],marker='.', facecolor='darkorange', edgecolor='darkorange',s=80)

gc.show_circles(co.convHMS('15:51:53.2'),co.convDMS('+19:11:03.860'),radius=1.5/60,linewidth=3, linestyle='--', color='black')

#gc.show_contour(h_850,hdu=1,levels = [4], colors = 'limegreen',linewidths=3)
#gc.show_contour(h_ch2,levels = [0.01], colors = 'red',linewidths=2)
#gc.show_contour(h_24,levels = [0.01], colors = 'blue',linewidths=2)


gc.add_scalebar(2.0/60,linewidth=3, corner='bottom left') #2'
gc.scalebar.set_label("2 arcmin")
gc.scalebar.set_color('black')
gc.save('../Figures/h1549_850_map_gray.pdf',dpi=250)
#pl.show()
pl.close()
gc.close()



#core
fig = pl.figure(figsize=(8,8),frameon=False)
pl.rc('font',size=22)
pl.rc('mathtext', default='regular')

RA_c = co.convHMS('15:51:53.2')
DEC_c = co.convDMS('+19:11:03.860')

gc = ap.FITSFigure(h_ch2, figure = fig)
#gc.tick_labels.set_xformat('hh:mm:ss')
#gc.tick_labels.set_yformat('dd:mm:ss')
gc.hide_xaxis_label()
gc.hide_xtick_labels()
gc.hide_yaxis_label()
gc.hide_ytick_labels()
gc.recenter(RA_c, DEC_c,width = 4.0/60, height= 4.0/60)
gc.show_colorscale(cmap='gray_r',vmax=0.1,vmin=-0.00025,interpolation='nearest') #vmax=0.0076,vmin=-0.00025,
#gc.add_colorbar(axis_label_text=r'S$_{850\mu m}$ (mJy)',axis_label_rotation=-90,axis_label_pad=10.5)

#gc.show_circles(co.convHMS('17:01:00.6'),co.convDMS('64:12:06.4'),radius = 15.0/3600,color='red')
gc.show_contour(h_850_snr,levels = [4,7,11], colors = 'limegreen',linewidths=2)
gc.show_contour(h_450_snr,levels = [3], colors = 'magenta',linewidths=2)

gc.show_circles(co.convHMS('15:51:53.8'), co.convDMS('+19:11:09.860'),radius = 2.5/3600,color='red',linewidth=3)
gc.add_label(co.convHMS('15:51:53.8')+22.0/3600, co.convDMS('+19:11:09.860'),'1_1',color='black',size=22)

gc.show_circles(co.convHMS('15:51:53.2'), co.convDMS('+19:10:59.100'),radius = 2.5/3600,color='red',linewidth=3)
gc.add_label(co.convHMS('15:51:53.2'), co.convDMS('+19:10:59.100')-12.0/3600,'1_2',color='black',size=22)

gc.show_circles(co.convHMS('15:51:52.5'), co.convDMS('+19:11:03.860'),radius = 2.5/3600,color='red',linewidth=3)
gc.add_label(co.convHMS('15:51:52.5')-17.0/3600, co.convDMS('+19:11:03.860'),'1_3',color='black',size=22)

for j in range(len(ID_IRAC)):
    if ID_IRAC[j][5:].split('_')[0] == '2':
        gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 2.5/3600,color='red',linewidth=3)
    if ID_IRAC[j][5:] == '3_1':
        gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 2.5/3600,color='blue',linewidth=3)
#        gc.add_label(RA_IRAC[j]-10.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'black',size=15)
#    if ID_IRAC[j][5:] == '10':
#        gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 2.5/3600,color='blue',linewidth=3)
#        gc.add_label(RA_IRAC[j]-10.0/3600,DEC_IRAC[j],ID_IRAC[j][5:],color = 'black',size=15)
    if ID_IRAC[j][5:] == '3_2' or ID_IRAC[j][5:].split('_')[0] == '24' or ID_IRAC[j][5:].split('_')[0] == '15' or ID_IRAC[j][5:].split('_')[0] == '21' or ID_IRAC[j][5:].split('_')[0] == '11' or ID_IRAC[j][5:].split('_')[0] == '9' or ID_IRAC[j][5:].split('_')[0] == '13' or ID_IRAC[j][5:] == '26_1':
        gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 2.5/3600,color='green',linewidth=3)


gc.show_circles(RA_c,DEC_c,radius = 1.5/60,color='black',linewidth=5, linestyle='--')
gc.add_label(RA_c,DEC_c+1.57/60+5.0/3600,'3-arcmin core',color = 'black',size=22)

gc.add_scalebar(0.5/60,linewidth=3, corner='bottom right') #0.5''
gc.scalebar.set_label("0.5 arcmin")
gc.scalebar.set_color('black')
gc.save('../Figures/h1549_ch2_core.pdf',dpi=250)
#pl.show()
pl.close()
gc.close()


#whole field - 450um
fig = pl.figure(figsize=(8,8),frameon=False)
pl.rc('font',size=15)
pl.rc('mathtext', default='regular')

RA_c = co.convHMS('15:51:54.604')
DEC_c = co.convDMS('19:11:15.90')

gc = ap.FITSFigure(h_450, figure = fig)
gc.tick_labels.set_xformat('hh:mm:ss')
gc.tick_labels.set_yformat('dd:mm:ss')
gc.recenter(RA_c, DEC_c,width = 12.0/60, height= 12.0/60)
gc.show_colorscale(cmap='gray_r',vmax=22,vmin=-10,interpolation='nearest')
gc.add_colorbar(axis_label_text=r'S$_{850\mu m}$ (mJy)',axis_label_rotation=-90,axis_label_pad=10.5)

for i in range(len(RA)):
    gc.show_circles(RA[i], DEC[i],radius = 15.0/3600,color='blue')
    gc.add_label(RA[i],DEC[i]+25.0/3600,ID[i],color = 'black',size=15)

gc.add_scalebar(2.0/60) #2'
gc.scalebar.set_label("2 arcmin")
gc.scalebar.set_color('black')
gc.save('../Figures/h1549_450_map.pdf',dpi=250)
pl.show()
pl.close()
gc.close()



ID_in = []
RA_in = []
DEC_in = []
z_in = []

p = []
p_ID = []
p_RA = []
p_DEC = []
p_theta = []

#850,450,24,ch4,ch2,J-band
for i in range(len(ID)+1)[1:]:

    # i=1

    pl.figure()
    fig = pl.figure()
    pl.rc('font',size=25)
    RA_c = RA[i-1]
    DEC_c = DEC[i-1]


    id850 = []
    idin = []
    rain = []
    decin = []
    zin = []

    #850
    gc = ap.FITSFigure(h_850, figure=fig, subplot=[0,0,0.75,1])
    gc.hide_xaxis_label()
    gc.hide_xtick_labels()
    gc.hide_yaxis_label()
    gc.hide_ytick_labels()
    gc.show_colorscale(cmap='gray_r',vmin = -3,vmax = max_850[i-1])
    gc.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
    gc.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
    for j in range(len(RA_op)):
        if abs(RA_op[j] - RA_c) < (30.0/3600)/2 and abs(DEC_op[j] - DEC_c) < (30.0/3600)/2:

            if 2.82 <= z[j] <= 2.92:
                gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='red',linewidth=2) #1''
                #gc.add_label(RA_op[j],DEC_op[j]+4E-4,ID_op[j],color = 'red')

                idin.append(ID_op[j])
                rain.append(RA_op[j])
                decin.append(DEC_op[j])
                zin.append(z[j])

            if z[j]<2.82:
                gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan', linewidth=2) #1''
                #gc.add_label(RA_op[j],DEC_op[j]+4E-4,ID_op[j],color = 'cyan')

            if z[j]>2.92:
                gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''
                #gc.add_label(RA_op[j],DEC_op[j]+4E-4,ID_op[j],color = 'cyan')


    gc.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w', linewidths = 2)
    if i!=1:
        for j in range(len(ID_IRAC)):
            if int(ID_IRAC[j][5:].split('_')[0]) == i:
                gc.show_circles(RA_IRAC[j], DEC_IRAC[j],radius = 1.0/3600,color='lime',linewidth=2)
    if i == 1:
        gc.show_contour(h_SMA,levels = [2.9,3.5,4.626], colors = ['yellow','yellow','yellow'],linewidths = 2)

    gc.add_scalebar(0.0)
    gc.scalebar.set_label(r'850 $\mu$m')
    gc.scalebar.set_color(color='black')


    ID_in.append(idin)
    RA_in.append(rain)
    DEC_in.append(decin)
    z_in.append(zin)

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
            if 2.82 <= z[j] <= 2.92:
                gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='red',linewidth=2) #1''

            elif z[j]<2.82:
                gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''

            elif z[j]>2.92:
                gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''

    gc.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
    if i == 1:
        gc.show_contour(h_SMA,levels = [2.9,3.5,4.626], colors = ['yellow','yellow','yellow'],linewidths = 2)
    else:
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
                theta = 3600*np.sqrt(((RA_c-RA_IRAC[j])*np.cos(19*np.pi/180.0))**2 + (DEC_c-DEC_IRAC[j])**2)
                n = 1.25E-3 #1/(np.pi*(7.5)**2)
                pval = 1-np.exp(-np.pi*n*theta**2)
                print ID_IRAC[j] + ' ' + str(round(theta,3)) + ' ' +str(round(pval,3))
                p.append(pval)
                p_ID.append(ID_IRAC[j])
                p_RA.append(RA_IRAC[j])
                p_DEC.append(DEC_IRAC[j])
                p_theta.append(theta)

        for j in range(len(RA_op)):
            if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:
                if 2.82 <= z[j] <= 2.92:
                    gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='red',linewidth=2) #1''

                elif z[j]<2.82:
                    gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''

                elif z[j]>2.92:
                    gc.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''


        gc.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
        if i == 1:
            gc.show_contour(h_SMA,levels = [2.9,3.5,4.626], colors = ['yellow','yellow','yellow'],linewidths = 2)
        else:
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
        image_ch4 = h_ch4
        gc_ch4 = ap.FITSFigure(image_ch4, figure=fig, subplot=[0,-1,0.75,1])
        gc_ch4.hide_xaxis_label()
        gc_ch4.hide_xtick_labels()
        gc_ch4.hide_yaxis_label()
        gc_ch4.hide_ytick_labels()
        gc_ch4.tick_labels.set_xformat('hh:mm:ss.s')
        gc_ch4.show_colorscale(cmap='gray_r',vmin = 1.41,vmax = max_ch4[i-1])
        gc_ch4.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
        gc_ch4.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
        gc_ch4.add_scalebar(0.0)
        gc_ch4.scalebar.set_label(r'8.0 $\mu$m')
        gc_ch4.scalebar.set_color(color='black')
        for j in range(len(RA_op)):
            if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

                if 2.82 <= z[j] <= 2.92:
                    gc_ch4.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='red',linewidth=2) #1''

                elif z[j]<2.82:
                    gc_ch4.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''

                elif z[j]>2.92:
                    gc_ch4.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''



        gc_ch4.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
        if i == 1:
            gc_ch4.show_contour(h_SMA,levels = [2.9,3.5,4.626], colors = ['yellow','yellow','yellow'],linewidths = 2)
        else:
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
        gc_ch4.scalebar.set_color(color='black')
        gc_ch4.add_label(RA_c,DEC_c,'no data',color='black')

    #4.5um
    if data_ch2_b[i-1]!=i:
        image_ch2 = h_ch2
        gc_ch2 = ap.FITSFigure(image_ch2, figure=fig, subplot=[0.75,-1,0.75,1])
        gc_ch2.tick_labels.set_xformat('hh:mm:ss.s')
        gc_ch2.hide_xaxis_label()
        gc_ch2.hide_xtick_labels()
        gc_ch2.hide_yaxis_label()
        gc_ch2.hide_ytick_labels()
        gc_ch2.show_colorscale(cmap='gray_r',vmin =-0.00025,vmax = max_ch2[i-1])
        gc_ch2.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
        gc_ch2.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
        gc_ch2.add_scalebar(0.0)
        gc_ch2.scalebar.set_label(r'4.5 $\mu$m')
        gc_ch2.scalebar.set_color(color='black')
        for j in range(len(RA_op)):
            if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

                if 2.82 <= z[j] <= 2.92:
                    gc_ch2.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='red',linewidth=2) #1''

                elif z[j]<2.82:
                    gc_ch2.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''

                elif z[j]>2.92:
                    gc_ch2.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''

#                gc_ch2.add_label(RA_op[j],DEC_op[j]+4E-4,ID_op[j],color = 'black',size=12)
        gc_ch2.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
        if i == 1:
            gc_ch2.show_contour(h_SMA,levels = [2.9,3.5,4.626], colors = ['yellow','yellow','yellow'],linewidths = 2)
            for j in range(len(ID_IRAC)):
                if int(ID_IRAC[j][5:].split('_')[0]) == i:
                    gc_ch2.add_label(RA_IRAC[j]+3e-4,DEC_IRAC[j]+4E-4,ID_IRAC[j],color = 'lime')
        else:
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
        gc_J.tick_labels.set_xformat('hh:mm:ss.s')
        gc_J.hide_xaxis_label()
        gc_J.hide_xtick_labels()
        gc_J.hide_yaxis_label()
        gc_J.hide_ytick_labels()
        gc_J.show_colorscale(cmap='gray_r',vmin = -2.5,vmax = max_J[i-1])
        gc_J.recenter(RA_c, DEC_c,width = 30.0/3600, height= 30.0/3600)
        gc_J.show_circles(RA_c, DEC_c,radius = 7.5/3600, linewidth = 3, linestyle='--', color='black')
        gc_J.add_scalebar(0.0)
        gc_J.scalebar.set_label(r'1.2 $\mu$m')
        gc_J.scalebar.set_color(color='black')
        for j in range(len(RA_op)):
            if RA_op[j] - RA_c < (30.0/3600)/2 and DEC_op[j] - DEC_c < (30.0/3600)/2:

                if 2.82 <= z[j] <= 2.92:
                    gc_J.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='red',linewidth=2) #1''

                elif z[j]<2.82:
                    gc_J.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''

                elif z[j]>2.92:
                    gc_J.show_circles(RA_op[j], DEC_op[j],radius = 0.75/3600, color='cyan',linewidth=2) #1''

        gc_J.show_contour(h_850_snr,levels = [4,5,6,8,10,12], colors = 'w',linewidths = 2)
        if i == 1:
            gc_J.show_contour(h_SMA,levels = [2.9,3.5,4.626], colors = ['yellow','yellow','yellow'],linewidths = 2)
        else:
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
        gc_J.scalebar.set_color(color='black')
        gc_J.add_label(RA_c,DEC_c,'no data',color='black')

    gc.save('../../Protocluster_paper/Paper_HS1549_HS1700/final_September17_2018_MNRAS/Figures/Cutouts/h1549/h1549_%i.pdf'%(i),dpi=50)

    gc.close()

    if data_24_b[i-1]!=i:
        gc.close()
    if data_ch4_b[i-1]!=i:
        gc_ch4.close()
    if data_ch2_b[i-1]!=i:
        gc_ch2.close()
    if data_J_b[i-1]!=i:
        gc_J.close()

#    pl.show()
    pl.close()
#    if i >0:


    # break


p_RA_new = []
for r in p_RA:
    r = co.deg2HMS(r)
    tmp = r.replace(r[-6:],str(round(float(r[-6:]),1)))
    p_RA_new.append(tmp.replace("15:", "$15^{\mathrm{h}}").replace(":","^{\mathrm{m}}")+"^{\mathrm{s}}$")

p_DEC_new = []
for d in p_DEC:
    d = co.deg2DMS(d)
    tmp = d.replace("19:", "$19^{\mathrm{o}}").replace(":","'")+"''$"
    p_DEC_new.append(tmp)

for i in range(len(p)):
    if len(p_ID[i].split('_')) == 2:
        print p_ID[i].split('_')[1] + '\t&\t' + p_RA_new[i] + '\t&\t' + p_DEC_new[i] + '\t&\t' + str(round(p_theta[i],2)) + '\t&\t' + str(round(p[i],3)) + '\\\\'
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


#for i in range(len(ID_in)):
#    if len(ID_in) == 0:
#        pass
#    else:
#        print '1549.' + ID[i]
#        for j in range(len(ID_in[i])):
#            print str(ID_in[i][j]) + '\t' + str(round(RA_in[i][j],3)) + '\t' + str(round(DEC_in[i][j],3)) + '\t' + str(round(z_in[i][j],3))
#    print
