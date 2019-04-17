import numpy as np
import matplotlib.pyplot as pl

# moving average function
def movingaverage(interval, window_size):
	window = np.ones(int(window_size))/float(window_size)
	return np.convolve(interval, window, 'same')

# Temp. ID names
IDs = np.array(['2bH-1d','2bK_1d','4aH_1d','4aJ2_1d','4aK_1d','16aK_1d'])
# Temp. spec. redshifts
z = np.array([2.31,2.31,2.31,2.31,2.31,2.31])

# Import Gemini spectra
wavelength = []
flux = []
for i in IDs:
    x,y = np.loadtxt('../../KBSS/h1700/Gemini/h1700_' + i + '.txt',unpack=True)
    wavelength.append(x)
    flux.append(y/max(y))
wavelength = np.array(wavelength)
flux = np.array(flux)


##############
# Plot spectra
##############
for i in range(len(IDs)):
    ax1 = pl.subplot(111)
    fig = pl.figure(figsize=(8,5))
    pl.rc('font',size=18)
    pl.rc('mathtext', default='regular')

    pl.plot(wavelength[i]/(1+z[i]), flux[i], c='b', ls='steps', label='1700.'+IDs[i][0:3]+'(z='+str(z[i])+')')

    pl.ylabel('Relative flux')
    pl.xlabel(r'Rest wavelength ($\AA$)')

    pl.xlim(6400,6900)
    pl.ylim(-0.4,1.2)
    pl.legend(fontsize=14,loc=0,ncol=1,frameon=True,numpoints=1,scatterpoints=1)

    pl.savefig('../Figures/HS1700_spectra/1700.'+IDs[i][0:3]+'.pdf', bbox_inches='tight')
    pl.close()
