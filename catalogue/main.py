# catalogue/main.py
# > ipython -pylab
# ALS 04/07/2015

"""
PURPOSE: cross matching ALPAKA (Mullaney+13) to external catalogue WISE, IRAF. 
		 And infer important quantities like Lbol.

		 Currently done for the Luminous (LOIII>5.e41) Nearby (z<0.2) objects
		 and output as Mullaney_allLN.fits

CAUTION: catalogue matching currently is a time consuming process
"""

import numpy as np
# utility functions
import catalogue_util
import crossmatch_util
import conversion_util

# from astropy.io import fits
# import numpy.lib.recfunctions as rfn
from astropy.table import Table, hstack, Column, join

# from astroquery.irsa import Irsa
# from astropy.coordinates import SkyCoord
import astropy.units as u

# from astropy.cosmology import luminosity_distance
# from astropy.cosmology import FlatLambdaCDM
# cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


def main():
	compileallLN() # done
	compileallvLOIIIr() # done
	compileMagellanObservable() # done
	compileall() # requires 3-4 days to run

	# T2 
	compileT2vLOIII() # done

	makeALPAKA_extended()  

def compileT2vLOIII():
	"""
	PURPOSE: compile table Mullaney_T2vLOIII.fits 
			 containing very Luminous (LOIII>1.e42) Type 2 of all z
			 that is Lbol_OIII > 2.6 e+45 erg / s

	Columns: 
		all the ALPAKA param, and the inferred OIII_5007T_LUM, OIII_5007AVG_FWHM
		WISE 
		IRAS
		inferred Luminosities: Lbol_OIII, LMIR, Lbol(MIRs), LFIR, LCO10EST, etc 
		consists of 360 objects
	"""
	# setups
	filein='ALPAKA_v1.fits'	
	fileout='Mullaney_T2vLOIII.fits'
	# read in Mullaney ALPAKA
	table=Table.read(filein,format='fits')

	# adding inferred OIII quantities
	conversion_util.inferALLinTable(table)

	#reduce table to the Luminous (LOIII>5.e41) and Near (z<0.2) - 935 objects
	selection=np.all([table['OIII_5007T_LUM']>1.e+42,table['AGN_TYPE']==2],axis=0)
	table=table[selection]

	print 'runnning ' + str(sum(selection))+' objects'
	# matching with external catalogue
	crossmatch_util.matchWISE(table)
	crossmatch_util.matchIRAS(table)

	# adding MIR, FIR infered quantities
	conversion_util.inferALLinTable(table)
	# add SDSSNAME
	catalogue_util.addcolSDSSNAME(table)

	table.write(fileout,format='fits')


def compileT2LOIII():
	"""
	PURPOSE: compile table Mullaney_T2LOIII.fits 
			 containing Luminous (LOIII>5.e41) Type 2 of all z

	Columns: 
		all the ALPAKA param, and the inferred OIII_5007T_LUM, OIII_5007AVG_FWHM
		WISE 
		IRAS
		inferred Luminosities: Lbol_OIII, LMIR, Lbol(MIRs), LFIR, LCO10EST, etc 
		consists of 909 objects
	"""
	# setups
	filein='ALPAKA_v1.fits'	
	fileout='Mullaney_T2LOIII.fits'
	# read in Mullaney ALPAKA
	table=Table.read(filein,format='fits')

	# adding inferred OIII quantities
	conversion_util.inferALLinTable(table)

	#reduce table to the Luminous (LOIII>5.e41) and Near (z<0.2) - 935 objects
	selection=np.all([table['OIII_5007T_LUM']>5.e+41,table['AGN_TYPE']==2],axis=0)
	table=table[selection]

	print 'runnning ' + str(sum(selection))+' objects'
	# matching with external catalogue
	crossmatch_util.matchWISE(table)
	crossmatch_util.matchIRAS(table)

	# adding MIR, FIR infered quantities
	conversion_util.inferALLinTable(table)
	# add SDSSNAME
	catalogue_util.addcolSDSSNAME(table)

	table.write(fileout,format='fits')


def compileT2mOIII():
	"""
	PURPOSE: compile table Mullaney_T2LOIII.fits 
			 containing Luminous (LOIII>1.e41) Type 2 of all z

	Columns: 
		all the ALPAKA param, and the inferred OIII_5007T_LUM, OIII_5007AVG_FWHM
		WISE 
		IRAS
		inferred Luminosities: Lbol_OIII, LMIR, Lbol(MIRs), LFIR, LCO10EST, etc 
		consists of 4460 objects
	"""
	# setups
	filein='ALPAKA_v1.fits'	
	fileout='Mullaney_T2mOIII.fits'
	# read in Mullaney ALPAKA
	table=Table.read(filein,format='fits')

	# adding inferred OIII quantities
	conversion_util.inferALLinTable(table)

	#reduce table to the Luminous (LOIII>5.e41) and Near (z<0.2) - 935 objects
	selection=np.all([table['OIII_5007T_LUM']>1.e+41,table['AGN_TYPE']==2],axis=0)
	table=table[selection]

	print 'runnning ' + str(sum(selection))+' objects'
	# matching with external catalogue
	crossmatch_util.matchWISE(table)
	crossmatch_util.matchIRAS(table)

	# adding MIR, FIR infered quantities
	conversion_util.inferALLinTable(table)
	# add SDSSNAME
	catalogue_util.addcolSDSSNAME(table)

	table.write(fileout,format='fits')



def compileallLN():
	"""
	compile table Mullaney_allLN.fits 
	containing all the Luminous (LOIII>5.e41) Nearby (z<0.2) object in ALPAKA

	Columns: 
		all the ALPAKA param, and the inferred OIII_5007T_LUM, OIII_5007AVG_FWHM
		WISE 
		IRAS
		inferred Luminosities: Lbol_OIII, LMIR, Lbol(MIRs), LFIR, LCO10EST, etc 
	"""
	# setups
	filein='ALPAKA_v1.fits'	
	fileout='Mullaney_allLN.fits'
	# read in Mullaney ALPAKA
	table=Table.read(filein,format='fits')

	# adding inferred OIII quantities
	conversion_util.inferALLinTable(table)

	#reduce table to the Luminous (LOIII>5.e41) and Near (z<0.2) - 935 objects
	selection=all([table['Z']<0.2, table['OIII_5007T_LUM']>5.e+41],axis=0)
	table=table[selection]

	# matching with external catalogue
	crossmatch_util.matchWISE(table)
	crossmatch_util.matchIRAS(table)

	# adding MIR, FIR infered quantities
	conversion_util.inferALLinTable(table)

	table.write(fileout,format='fits')


def compileallvLOIIIr():
	"""
	PURPOSE: compile table Mullaney_allvLOIIIr.fits 
			 containing all the very Luminous (LOIII>1.e42) ones with [OIII] lie well in r-band,
			 which is  0.1346 < z < 0.3434
			 contains 969 objects

	Columns: 
		all the ALPAKA param, and the inferred OIII_5007T_LUM, OIII_5007AVG_FWHM
		WISE 
		IRAS
		inferred Luminosities: Lbol_OIII, LMIR, Lbol(MIRs), LFIR, LCO10EST, etc 
	"""
	# setups
	filein='ALPAKA_v1.fits'	
	fileout='Mullaney_allvLOIIIr.fits'
	# read in Mullaney ALPAKA
	table=Table.read(filein,format='fits')

	# adding inferred OIII quantities
	conversion_util.inferALLinTable(table)

	#reduce table to the Luminous (LOIII>5.e41) and Near (z<0.2) - 935 objects
	selection=np.all([table['Z']>0.1346,table['Z']< 0.3434, table['OIII_5007T_LUM']>1.e+42],axis=0)
	table=table[selection]

	# matching with external catalogue
	crossmatch_util.matchWISE(table)
	crossmatch_util.matchIRAS(table)

	# adding MIR, FIR infered quantities
	conversion_util.inferALLinTable(table)

	table.write(fileout,format='fits')


def compileMagellanObservable():
	"""
	compile table Mullaney_Magellan_obsable.fits 
	containing all objects in ALPAKA observable by Magellan in 2014 June (RA < 2 or >10, DEC<15)

	Columns: 
		all the ALPAKA param, and the inferred OIII_5007T_LUM, OIII_5007AVG_FWHM
		WISE 
		IRAS
		inferred Luminosities: Lbol_OIII, LMIR, Lbol(MIRs), LFIR, LCO10EST, etc 
	"""
	# setups
	filein='ALPAKA_v1.fits'	
	fileout='Mullaney_Magellan_obsable.fits'
	# read in Mullaney ALPAKA
	table=Table.read(filein,format='fits')

	# adding inferred OIII quantities
	conversion_util.inferALLinTable(table)

	#reduce table Magellan observable on June 2014 (RA < 2 or >10, DEC<15) - 7280 objects
	selRA=any([table['RA']<2.*15., table['RA']>10.*15.],axis=0)
	selection=all([selRA, table['DEC']<15.],axis=0)
	table=table[selection]

	# matching with external catalogue
	crossmatch_util.matchWISE(table)
	crossmatch_util.matchIRAS(table)

	# adding MIR, FIR infered quantities
	conversion_util.inferALLinTable(table)

	table.write(fileout,format='fits')




def compileall():
	"""
	compile table Mullaney_all.fits 
	containing all objects in ALPAKA

	Columns: 
		all the ALPAKA param, and the inferred OIII_5007T_LUM, OIII_5007AVG_FWHM
		WISE 
		IRAS
		inferred Luminosities: Lbol_OIII, LMIR, Lbol(MIRs), LFIR, LCO10EST, etc 
	"""
	# setups
	filein='ALPAKA_v1.fits'	
	fileout='Mullaney_all.fits'
	filetemp='Mullaney_all_temp.fits'
	# read in Mullaney ALPAKA
	table=Table.read(filein,format='fits')

	# adding inferred OIII quantities
	conversion_util.inferALLinTable(table)

	# matching with external catalogue
	crossmatch_util.matchWISE(table)
	table.write(filetemp,format='fits')
	crossmatch_util.matchIRAS(table)
	table.write(filetemp,format='fits')
	# adding MIR, FIR infered quantities

	conversion_util.inferALLinTable(table)

	table.write(fileout,format='fits')



def makeALPAKA_extended():
	"""
	add in infered OIII quantities and write as ALPAKA_extended
	"""
	filein='ALPAKA_v1.fits'	
	fileout='ALPAKA_extended.fits'
	# read in Mullaney ALPAKA
	table=Table.read(filein,format='fits')

	# adding inferred OIII quantities
	conversion_util.inferALLinTable(table)
	table.write(fileout,format='fits')



if __name__ == '__main__':
	main()

