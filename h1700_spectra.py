import numpy as np
import matplotlib.pyplot as pl

# moving average function
def movingaverage(interval, window_size):
	window = np.ones(int(window_size))/float(window_size)
	return np.convolve(interval, window, 'same')

'''
old -> new ID names w/ GNIRS spec-z:

ld		new		z
---		---		-
2b		4.1		2.318
4a		7.1		2.313
5.X		5.2		2.303
13.Y	16.1	1.575
16a		17.1	2.306
'''

# Temp. ID names
IDs = np.array(['4','5_2','7_1','16','17'])
# Temp. spec. redshifts
z = np.array([2.318,2.286,2.313,1.575,2.306])

# Import Gemini spectra
wavelength = []
flux = []
for i in IDs:
    x,y = np.loadtxt('../../KBSS/h1700/Gemini/h1700_' + i + '.txt',unpack=True)
    y = movingaverage(y,5) # bin data
    wavelength.append(x)
    flux.append(y/max(y)) # relative flux to maximum flux
wavelength = np.array(wavelength)
flux = np.array(flux)


lines_wavelength = [6302.046,6365.536,6529.03, 6549.86, 6564.61, 6585.27, 6718.29, 6732.67]
lines_names = ['[OI]','[OI]','[NI]','[NII]', r'H$\alpha$','[NII]','[SII]','[SII]']
name_height = np.array([1.05,1.05,1.05,0.93,0.85,1.05,1.05,0.90])+0.28
name_offset = np.array([-10,-10,-15,-15,-5,-10,-15,-5])

##############
# Plot spectra
##############
for i in range(len(IDs)):

	ax1 = pl.subplot(111)
	fig = pl.figure(figsize=(8,5))
	pl.rc('font',size=18)
	pl.rc('mathtext', default='regular')

	pl.plot(wavelength[i]/(1+z[i]), flux[i], c='k',ls='steps')
	pl.plot(-100,-100,c='w',label='1700.'+IDs[i][0:3]+' (z='+str(z[i])+')')
	# Filled bellow spectra
	#d = np.zeros(len(y))
	pl.fill_between(wavelength[i]/(1+z[i]), flux[i], interpolate=True, color='k',alpha=0.2)
	#pl.fill_between(xs, ys, where=ys<=d, interpolate=True, color='red')

	pl.axhline(y=0, ls=':', c='k')

	for j in range(len(lines_wavelength)):
		pl.axvline(x=lines_wavelength[j],c='r',ls='--')
		pl.text(lines_wavelength[j]+name_offset[j],name_height[j],lines_names[j],fontsize=10,rotation=90,bbox=dict(facecolor='white', edgecolor='none',boxstyle='round',pad=0))



	if abs(z[i]-2.3)<0.1:
		pl.xlim(6000,7200)

	elif abs(z[i]-1.5)<0.1:
		pl.xlim(6000,7200)

	elif abs(z[i]-2.8)<0.1:
		pl.xlim(6000,6500)

	pl.ylim(-0.2,1.4)

	pl.ylabel('Relative flux')
	pl.xlabel(r'Rest wavelength ($\AA$)')

	pl.rcParams['legend.handlelength'] = 0
	pl.rcParams['legend.numpoints'] = 1

	pl.legend(fontsize=12,loc=0,ncol=1,scatterpoints=0,frameon=True)

	pl.savefig('../Figures/HS1700_spectra/1700_'+IDs[i][0:3]+'.pdf', bbox_inches='tight')
	pl.close()
