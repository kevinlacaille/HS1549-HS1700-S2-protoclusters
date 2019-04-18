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
16a		17.1	2.306
5.X		5.2		2.303
13.Y	16.1	1.575
'''

# Temp. ID names
IDs = np.array(['2bH-1d','2bK_1d','4aH_1d','4aJ2_1d','4aK_1d','16aK_1d'])
# Temp. spec. redshifts
z = np.array([None,2.318,None,None,2.31,2.3105])

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

lines_wavelength = [6529.03, 6549.86, 6564.61, 6585.27, 6718.29, 6732.67]
lines_names = ['[N I]','[N II]', r'H$\alpha$','[N II]','[S II]','[S II]']
##############
# Plot spectra
##############
for i in range(len(IDs)):
	if IDs[i] == '2bH-1d' or IDs[i] == '4aH_1d' or IDs[i] == '4aJ2_1d':
		pass
	else:
		ax1 = pl.subplot(111)
		fig = pl.figure(figsize=(8,5))
		pl.rc('font',size=18)
		pl.rc('mathtext', default='regular')

		pl.plot(wavelength[i]/(1+z[i]), flux[i], c='b', label='1700.'+IDs[i][0:3]+'(z='+str(z[i])+')')

		for j in range(len(lines_wavelength)):
			pl.axvline(x=lines_wavelength[j],c='k',ls='--')
			# pl.text(lines_wavelength[j]+5,1,lines_names[j])


		pl.ylabel('Relative flux')
		pl.xlabel(r'Rest wavelength ($\AA$)')

		if abs(z[i]-2.3)<0.1:
			pl.xlim(6400,6900)
		# elif abs(z[i]-1.5)<0.1:
		#     pl.xlim(4700,5100)

		pl.ylim(-0.4,1.2)
		pl.legend(fontsize=14,loc=0,ncol=1,frameon=True,numpoints=1,scatterpoints=1)

		pl.savefig('../Figures/HS1700_spectra/1700_'+IDs[i][0:3]+'.pdf', bbox_inches='tight')
		pl.close()
