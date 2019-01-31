import numpy as np
import pylab as pl
from astropy.io import fits
import scipy
from scipy import ndimage
import coords as co
from aplpy import FITSFigure
from progress.bar import FillingCirclesBar

def cutout(x,y,pix):
    return x-15/float(pix),x+15/float(pix),y-15/float(pix),y+15/float(pix) #30arcsec X 30arcsec cutout, +/-15arcsec about centre point (x,y), 1pixel=4arcsec (850um), 2arcsec (450um)


wavelength = ['850']#['850','450']

for w in wavelength:
    
    cat = open('../../Data/Catalogues_new/q1700_all_radec_python_' + w + '.gaia','r')

    ID = []
    RA = []
    DEC = []
    z = []

    for line in cat.readlines()[2:]:
        tmp = line.split()
        ID.append(tmp[0])
        RA.append(co.convHMS(tmp[1]))
        DEC.append(co.convDMS(tmp[2]))
        z.append(float(tmp[4]))

    ID = np.array(ID)
    RA = np.array(RA)
    DEC = np.array(DEC)
    z = np.array(z)


    #protocluster z-analysis
    n_histo = 75
    bins = np.histogram(z,n_histo)[1]
    z_min = np.histogram(z,n_histo)[1][28]
    z_max = np.histogram(z,n_histo)[1][32]

    N_proto = np.sum(np.histogram(z,n_histo)[0][28:32])
    N = np.sum(np.histogram(z,n_histo)[0])
    print str(N_proto) + ' galaxies (' + str(round(float(N_proto)/N*100,2)) + '%) within ' + str(round(z_min,3)) + '<z<' + str(round(z_max,3))

    ID_proto = ID[np.where((z<=z_max) & (z>=z_min))[0]]
    RA_proto = RA[np.where((z<=z_max) & (z>=z_min))[0]]
    DEC_proto = DEC[np.where((z<=z_max) & (z>=z_min))[0]]
    z_proto = z[np.where((z<=z_max) & (z>=z_min))[0]]

    cat_proto = open('../../Data/Catalogues_new/proto_' + w + '_z.gaia','w')
    cat_proto.write('ID\tRA\tDEC\tz\n--\t--\t---\t-\n')
    for i in range(len(ID_proto)):
        cat_proto.write(str(ID_proto[i]) + '\t' + str(RA_proto[i]) + '\t' + str(DEC_proto[i]) + '\t' + str(z_proto[i]) + '\n')
    cat.close()
    cat_proto.close()
    



    print w + 'um'

    #850um data and coordinates
    f = fits.open('../../Data/SCUBA2/h1700_' + w + '.fits')
    scidata = np.nan_to_num(f[0].data[0])

    rotated = scipy.ndimage.interpolation.rotate(scidata,180)
    image = scipy.fliplr(rotated)

    gc = FITSFigure('../../Data/SCUBA2/h1700_' + w + '.fits')
    (X_proto,Y_proto) = gc.world2pixel(RA_proto,DEC_proto)[0]-1,len(image)-gc.world2pixel(RA_proto,DEC_proto)[1] #X-1 because array shifted from 0, height-Y to flip, no effect from 180deg rotation
    gc.close()


    #30''x30'' cutouts
    if w == '850':
        cutoutsize = 7 #28x28 cutout
    elif w == '450':
        cutoutsize = 15 #30x30 cutout
    total = np.array(cutoutsize*[np.zeros(cutoutsize)]) #8 = len(image[y_d:y_u,x_l:x_r])
    total_color = np.array(cutoutsize*[np.zeros(cutoutsize)])
    total_NB = np.array(cutoutsize*[np.zeros(cutoutsize)])
    index_bad = []
    index_color = []
    index_NB = []

    bar = FillingCirclesBar('Cutouts',max = len(ID_proto))

    for i in range(len(ID_proto)):

        if ID_proto[i][0] != '*':
            
            bar.next()
            
            if w == '850':
                pix = 4 #number of arcseconds/pixel
            elif w == '450':
                pix = 2
            (x_l,x_r,y_d,y_u) = cutout(X_proto[i],Y_proto[i],pix)
            #need shape to be (8,8) and centred about ID_proto[i]
            if image[y_d:y_u,x_l:x_r].shape[0] != cutoutsize: #(y,x), if y!=8
                if np.mean([y_u,y_d]) <= Y_proto[i]: #if <y_image> <= y_proto, cutout centred too low, move cutout up
                    if w == '850':
                        y_d +=1 #shift the bottom bound lower
                    elif w =='450':
                        y_d -=1 #shift the bottom bound lower
                else: #<y_image> > y_proto, cutout too high, move cutout down
                    if w == '850':
                        y_u +=1 #shift the upper bound down
                    elif w == '450':
                        y_u -=1 #shift the upper bound down
            if image[y_d:y_u,x_l:x_r].shape[1] != cutoutsize: #(y,x), if x!=8
                if np.mean([x_l,x_r]) <= X_proto[i]: #if <x_image> <= x_proto, cutout centred too left, move cutout right
                    if w == '850':
                        x_r -=1 #shift the right bound right
                    elif w == '450':
                        x_r +=1
                else: #<x_image> > x_proto, cutout too right, move cutout left
                    if w == '850':
                        x_l -=1 #shift the left bound right
                    elif w == '450':
                        x_l +=1
            if image[y_d:y_u,x_l:x_r].shape != (cutoutsize,cutoutsize):
                print "cutout shape is still wrong! " + str(image[y_d:y_u,x_l:x_r].shape)
                break

            total+=image[y_d:y_u,x_l:x_r]

            if w == '850':
                centre = 3#?
            elif w == '450':
                centre = 7

            '''
            fig = pl.figure()
            pl.rc('font',size=15)
            pl.rc('mathtext', default='regular')
            ax = fig.add_subplot(1, 1, 1)


            imgplot = pl.imshow(image[y_d:y_u,x_l:x_r],cmap='YlOrBr_r')
            imgplot.set_interpolation('nearest')
            
            circ = pl.Circle((centre,centre), radius=7.5/pix,facecolor='none',edgecolor='k')
            ax.add_patch(circ)

            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_xaxis().set_ticks([])
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_yaxis().set_ticks([])
            
            pl.colorbar()

            pl.savefig('All/Cutouts/' + w + '/' + ID_proto[i] + '.png',bbox_inches='tight')
            #pl.show()
            pl.close()
            '''


            #seperate measurments for color and narrow band selected sources
            if ID_proto[i][0:2] == 'BX' or ID_proto[i][0:2] == 'MD' or ID_proto[i][0:1] == 'C' or ID_proto[i][0:1] == 'D' or ID_proto[i][0:1] == 'C' or ID_proto[i][0:1] == 'M' or ID_proto[i][0:2] == 'BM':
                index_color.append(i)
                total_color+=image[y_d:y_u,x_l:x_r]
                
                fig = pl.figure()
                pl.rc('font',size=15)
                pl.rc('mathtext', default='regular')
                ax = fig.add_subplot(1, 1, 1)
                
                
                imgplot = pl.imshow(image[y_d:y_u,x_l:x_r],cmap='YlOrBr_r')
                imgplot.set_interpolation('nearest')
                
                circ = pl.Circle((centre,centre), radius=7.5/pix,facecolor='none',edgecolor='k')
                ax.add_patch(circ)
                
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_xaxis().set_ticks([])
                ax.axes.get_yaxis().set_visible(False)
                ax.axes.get_yaxis().set_ticks([])
                
                pl.colorbar()
                
                pl.savefig('Colour/Cutouts/' + w + '/' + ID_proto[i] + '.png',bbox_inches='tight')
                #pl.show()
                pl.close()
                

            else:
                index_NB.append(i)
                total_NB+=image[y_d:y_u,x_l:x_r]
                
                fig = pl.figure()
                pl.rc('font',size=15)
                pl.rc('mathtext', default='regular')
                ax = fig.add_subplot(1, 1, 1)
                
                
                imgplot = pl.imshow(image[y_d:y_u,x_l:x_r],cmap='YlOrBr_r')
                imgplot.set_interpolation('nearest')
                
                circ = pl.Circle((centre,centre), radius=7.5/pix,facecolor='none',edgecolor='k')
                ax.add_patch(circ)
                
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_xaxis().set_ticks([])
                ax.axes.get_yaxis().set_visible(False)
                ax.axes.get_yaxis().set_ticks([])
                
                pl.colorbar()
                
                pl.savefig('NB/Cutouts/' + w + '/' + ID_proto[i] + '.png',bbox_inches='tight')
                #pl.show()
                pl.close()
                
        else:
            bar.next()
            index_bad.append(i)

    bar.finish()


    #map
    fig = pl.figure()
    pl.rc('font',size=15)
    pl.rc('mathtext', default='regular')
    ax = fig.add_subplot(1, 1, 1)

    if w == '850':
        imgplot = pl.imshow(image,vmin=-1.9,vmax = 4.0,cmap='YlOrBr_r')
    elif w == '450':
        imgplot = pl.imshow(image,vmin=-17.0,vmax = 20.0,cmap='YlOrBr_r')

    imgplot.set_interpolation('nearest')

    for i in index_color:
        circ = pl.Circle((X_proto[i],Y_proto[i]), radius=cutoutsize/4.0,facecolor='none',edgecolor='b')
        ax.add_patch(circ)
    for i in index_NB:
        circ = pl.Circle((X_proto[i],Y_proto[i]), radius=cutoutsize/4.0,facecolor='none',edgecolor='c')
        ax.add_patch(circ)
    for i in index_bad:
        circ = pl.Circle((X_proto[i],Y_proto[i]), radius=cutoutsize/4.0,facecolor='none',edgecolor='r')
        ax.add_patch(circ)

    pl.colorbar()

    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_yaxis().set_ticks([])


    pl.savefig(w + '_field.png',dpi=250,bbox_inches='tight')
    #pl.show()
    pl.close()



    #stack
    len_all = len(ID_proto) - len(index_bad)
    len_color = len(index_color)
    len_NB = len(index_NB)
    len_bad = len(index_bad)

    total/=(len_all)
    total_color/=(len_color)
    total_NB/=(len_NB)

    #z-histo
    fig = pl.figure()
    fig.subplots_adjust(hspace=0,wspace=0)
    pl.rc('font',size=15)
    pl.rc('mathtext', default='regular')

    pl.hist(z,bins=bins,color='w',hatch = "//",label='all')
    pl.hist(z[index_color],bins=bins,color='b',alpha=0.5,label='color')
    pl.hist(z[index_NB],bins=bins,color='c',alpha=0.5,label='NB')
    pl.hist(z[index_bad],bins=bins,color='r',alpha=0.5,label='bad')

    pl.axvline(x = z_min, c='k',ls='--')
    pl.axvline(x = z_max, c='k',ls='--')

    pl.xlabel('z')
    pl.ylabel('N')

    pl.legend(fontsize=13,loc=0)#frameon=False,numpoints=1,fontsize=13,loc=0)
    pl.savefig(w + '_z_proto.pdf',bbox_inches='tight')
    pl.close()



    pic = [total,total_color,total_NB]

    #850um stack cutout
    for i in range(len(pic)):
        
        fig = pl.figure()
        pl.rc('font',size=15)
        pl.rc('mathtext', default='regular')
        ax = fig.add_subplot(1, 1, 1)

        imgplot = pl.imshow(pic[i],cmap='YlOrBr_r')
        imgplot.set_interpolation('nearest')

        circ = pl.Circle((centre,centre), radius=7.5/pix,facecolor='none',edgecolor='k')
        ax.add_patch(circ)

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_yaxis().set_ticks([])

        pl.colorbar()

        if i == 0:
            pl.savefig('All/' + w + '_all_stack.png',bbox_inches='tight')
        elif i == 1:
            pl.savefig('Colour/' + w + '_colour_stack.png',bbox_inches='tight')
        elif i == 2:
            pl.savefig('NB/' + w + '_NB_stack.png',bbox_inches='tight')
        #pl.show()
        pl.close()


        #850um stack profile
        fig = pl.figure(figsize = (10,10))
        fig.subplots_adjust(hspace=0,wspace=0)
        pl.rc('font',size=15)
        pl.rc('mathtext', default='regular')

        if w == '850':
            sigma = 0.8
            ylabel = r'S$_{850\mu m}$ (mJy)'
            xlabel = r'S$_{850\mu m}$ (mJy)'
        elif w == '450':
            sigma = 10.0 #? 6mJy in ~centre 4''
            ylabel = r'S$_{450\mu m}$ (mJy)'
            xlabel = r'S$_{450\mu m}$ (mJy)'

        ax = fig.add_subplot(2, 2, 1)
        ax.plot(np.arange(-int(cutoutsize/2),cutoutsize/2+1),pic[i][cutoutsize/2:][0]) #plot 5th row down, left to right
        if i == 0:
            ax.axhline(y = sigma/np.sqrt(len_all),c='r',ls='-') #sigma~RMS/sqrt(N)
        elif i == 1:
            ax.axhline(y = sigma/np.sqrt(len_color),c='r',ls='-')
        elif i == 2:
            ax.axhline(y = sigma/np.sqrt(len_NB),c='r',ls='-')

        ax.set_xlabel('x-offset (left to right)')
        ax.set_ylabel(ylabel)


        ax2 = fig.add_subplot(2, 2,2)
        ax2.plot(pic[i][:,cutoutsize/2],np.arange(-int(cutoutsize/2),cutoutsize/2+1)) #plot 5th from the left, top down
        if i == 0:
            ax2.axvline(x = sigma/np.sqrt(len_all),c='r',ls='-')
        if i == 1:
            ax2.axvline(x = sigma/np.sqrt(len_color),c='r',ls='-')
        if i == 2:
            ax2.axvline(x = sigma/np.sqrt(len_NB),c='r',ls='-')

        ax2.set_ylabel('y-offset (top down)')
        ax2.set_xlabel(xlabel)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        if i == 0:
            pl.savefig('All/' + w + '_profile_all.pdf',bbox_inches='tight')
        if i == 1:
            pl.savefig('Colour/' + w + '_profile_colour.pdf',bbox_inches='tight')
        if i == 2:
            pl.savefig('NB/' + w + '_profile_NB.pdf',bbox_inches='tight')
        #pl.show()
        pl.close()


