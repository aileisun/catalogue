# crossmatch_utils.py
# ALS 04/07/2015
"""
Contains functions for cross-matching external catalogues

CONTENT:
    addWISE_Lum_rf(table, w_rf=8.0)
    
    correctWISE(table, selection=False)
    
    matchHerschelData_manual(table, data)
    
    matchIRAS(table, selection=False)
    
    matchNEDNAME(table, selection=False)
    
    matchWISE(table, selection=False)
    
    writeFIRSTsearch(filename, table)
"""

from pylab import *
import os
import conversion_util
reload(conversion_util)
import catalogue_util
reload(catalogue_util)

from astropy.table import Table, hstack, Column


from astroquery.irsa import Irsa
from astroquery.ned import Ned
from astropy.coordinates import SkyCoord
import astropy.units as u

from astropy import constants as const
# from astropy.cosmology import luminosity_distance
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)



#========= small functions
def lambda2nu(wavelength):
	return (const.c/wavelength).to('Hz')

#========= cross-matching with other catalogue

# def matchWISE(table,selection=False,searchradius=5.,catalog='wise_allwise_p3as_psd'):
	# """
	# PURPOSE: This funciton inquery WISE photometries given a table of targets with RA, DEC,  and attach the result to the table. 

	# DESCRIPTION: 
	# 	Columns created:
	# 		Qquery columns
	# 		['w1mpro','w2mpro','w3mpro','w4mpro',]
	# 		['w1sigmpro','w2sigmpro','w3sigmpro','w4sigmpro',]
	# 		['w1snr','w2snr','w3snr','w4snr',]


	# INPUT: 
	# 	table: a table with at least two columns 'RA', 'DEC'
	# 	selection: a bool array with length equals to table specifying what rows to match
	# 	searchradius: radius of search in unit of arcsec
	# 	catalog: which wise catalog to use
	# 		wise_allwise_p3as_psd: 		ALLWISE data release Nov. 2013 (default)
	# 		wise_allsky_4band_p3as_psd: ALLSKY data release March 2012

	# NOTE: 
	# 	1. Vega mag -> Fnu

	# 		Conversion from WISE catalog Vega mag to Fnu in cgs unit is based on 
	# 		"WISE All-Sky release explanatory supplement: Data Processing",
	# 		http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#example

	# 		adapting conversion of flat spectrum (F_nu proportional to nu^0),
	# 		which gives Vega zero-magnitude flux density:  (see Tab. 1 2nd column)

	# 			F_nu0		W1		W2		W3		W4
	# 			[Jy]		309.540	171.787	31.674	8.363

	# 		Conversion is done using formula: 

	# 			F_nu = F_nu0 * 10^(-magVG/2.5)

	# 		Adoption of different spectral shape causes differences in the zero-point mag. 
	# 		For example, if F_nu is proportional to nu^-2 then F_nu0 would differ by 9% in W3. 

	# 	2. Vega mag -> AB mag

	# 		Conversion from WISE Vega mag to AB mag, assuming flat spectrum, follows: 

	# 		m_ab = m_vega + deltam 

	# 			deltam	W1		W2		W3		W4
	# 					2.699	3.339	5.174	6.620

	# """

	# # list of columns to be accessed from query
	# query_colnames_mpro=['w1mpro','w2mpro','w3mpro','w4mpro',]
	# query_colnames_sigmpro=['w1sigmpro','w2sigmpro','w3sigmpro','w4sigmpro',]
	# query_colnames_snr=['w1snr','w2snr','w3snr','w4snr',]
	# WISE_colnames_others=['WISE_sep','WISE_nsource']


	# # list of columns to be created
	# WISE_colnames_magVG=['W1_magVG','W2_magVG','W3_magVG','W4_magVG']
	# WISE_colnames_magAB=['W1_magAB','W2_magAB','W3_magAB','W4_magAB']
	# WISE_colnames_Fnu=['W1_Fnu','W2_Fnu','W3_Fnu','W4_Fnu']
	# WISE_colnames_Lnu=['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']
	# WISE_colnames_mag_err=['W1_mag_err','W2_mag_err','W3_mag_err','W4_mag_err']
	# WISE_colnames_dex_err=['W1_dex_err','W2_dex_err','W3_dex_err','W4_dex_err']
	# WISE_colnames_snr=['W1_snr','W2_snr','W3_snr','W4_snr']


	# #==== selection
	# if type(selection)==bool:
	# 	if selection ==False: selection=ones(len(table))
	# #==== check if table has WISE columns, if not creat them
	# # 
	# for WISE_colname in WISE_colnames_magVG+WISE_colnames_magAB+WISE_colnames_Fnu+WISE_colnames_Lnu+WISE_colnames_mag_err+WISE_colnames_dex_err+WISE_colnames_snr+WISE_colnames_others:
	# 	if WISE_colname not in table.dtype.names:
	# 		print "creating column ", WISE_colname
	# 		table[WISE_colname]=0.
	# 		if WISE_colname in  WISE_colnames_Fnu:
	# 			table[WISE_colname].unit="erg s-1 cm-2 Hz-1"
	# 		if WISE_colname in WISE_colnames_Lnu:
	# 			table[WISE_colname].unit="erg s-1"
	# table['WISE_sep'].unit='arcsec'


	# #==== matching with ALLWISE
	# for i in range(len(table)):
	# 	if selection[i]:
	# 		magVGs=zeros(4)
	# 		magABs=zeros(4)
	# 		Fnus=zeros(4)*u.erg/u.cm**2/u.s/u.Hz
	# 		# calculate distance modulo
	# 		lum_dist=cosmo.luminosity_distance(table['Z'][i]).to(u.pc)
	# 		# dist_mod=-5.*log10(lum_dist/(10*u.pc)).value
	# 		# retrieving magAB
	# 		c = SkyCoord(table['RA'][i], table['DEC'][i], 'icrs', unit='deg')
	# 		tabwise=Irsa.query_region(c,catalog=catalog, spatial='Cone',radius=searchradius*u.arcsec)[['designation','ra','dec','w1mpro','w1sigmpro','w1snr','w2mpro','w2sigmpro','w2snr','w3mpro','w3sigmpro','w3snr','w4mpro','w4sigmpro','w4snr']]
	# 		if len(tabwise)==1:
	# 			for j in range(4):
	# 				# converting magVega to magAB and then Fnu
	# 				magVGs[j]=tabwise[query_colnames_mpro[j]]
	# 				magABs[j]=magVGs[j]+deltam[j]
	# 				Fnus[j]=(F_nu0[j]*10**(-magVGs[j]/2.5)).cgs
	# 				# Fnus[j]=conversion_util.magAB2Fnu(magABs[j])
	# 				# logging quantities to table
	# 				table[WISE_colnames_magVG[j]][i]=magVGs[j]
	# 				table[WISE_colnames_magAB[j]][i]=magABs[j]
	# 				table[WISE_colnames_Fnu[j]][i]=Fnus[j].value
	# 				table[WISE_colnames_Lnu[j]][i]=(nus[j]*Fnus[j]*4.*pi*lum_dist**2).to('erg/s').value
	# 				table[WISE_colnames_mag_err[j]][i]=tabwise[query_colnames_sigmpro[j]]
	# 				table[WISE_colnames_dex_err[j]][i]=tabwise[query_colnames_sigmpro[j]]/2.5
	# 				table[WISE_colnames_snr[j]][i]=tabwise[query_colnames_snr[j]]

	# 			# WISE_sep: distance between SDSS and WISE coordinates
	# 			cwise = SkyCoord(tabwise['ra'], tabwise['dec'], 'icrs', unit='deg')
	# 			table['WISE_sep'][i] = c.separation(cwise).arcsec
	# 			table['WISE_nsource'][i]=1.
	# 			print i, table['OBJID'][i], " matched"
	# 		elif len(tabwise)>=2:
	# 			# calculate distance modulo
	# 			# determining the closest object
	# 			cwise1 = SkyCoord(tabwise['ra'][0], tabwise['dec'][0], 'icrs', unit='deg')
	# 			cwise2 = SkyCoord(tabwise['ra'][1], tabwise['dec'][1], 'icrs', unit='deg')
	# 			dist1=c.separation(cwise1).arcsec
	# 			dist2=c.separation(cwise2).arcsec
	# 			k= dist1>dist2
	# 			for j in range(4):
	# 				# converting magVega to magAB and then Fnu
	# 				magVGs[j]=tabwise[query_colnames_mpro[j]][k]
	# 				magABs[j]=magVGs[j]+deltam[j]
	# 				Fnus[j]=(F_nu0[j]*10**(-magVGs[j]/2.5)).cgs
	# 				# Fnus[j]=conversion_util.magAB2Fnu(magABs[j])
	# 				# logging quantities to table
	# 				table[WISE_colnames_magVG[j]][i]=magVGs[j]
	# 				table[WISE_colnames_magAB[j]][i]=magABs[j]
	# 				table[WISE_colnames_Fnu[j]][i]=Fnus[j].value
	# 				table[WISE_colnames_Lnu[j]][i]=(nus[j]*Fnus[j]*4.*pi*lum_dist**2).to('erg/s').value
	# 				table[WISE_colnames_mag_err[j]][i]=tabwise[query_colnames_sigmpro[j]][k]
	# 				table[WISE_colnames_dex_err[j]][i]=tabwise[query_colnames_sigmpro[j]][k]/2.5
	# 				table[WISE_colnames_snr[j]][i]=tabwise[query_colnames_snr[j]][k]
	# 			table['WISE_sep'][i] = min(dist1,dist2)
	# 			table['WISE_nsource'][i]=len(tabwise)
	# 			print i, table['OBJID'][i], " confusion"
	# 		elif len(tabwise)==0:
	# 			table['WISE_sep'][i]=-1
	# 			table['WISE_nsource'][i]=0
	# 			print i, table['OBJID'][i], " no match"
	# 		else:
	# 			print i, table['OBJID'][i], " exception"

def matchWISE(table,selection=False,searchradius=5.,catalog='ALLWISE', updatelocal=False):
	"""
	PURPOSE: Infer WISE relevent columns from a table with WISE inquery and redshift columns
			Inferred columns: 
			['W1_magVG','W2_magVG','W3_magVG','W4_magVG']
			['W1_magAB','W2_magAB','W3_magAB','W4_magAB']
			['W1_Fnu','W2_Fnu','W3_Fnu','W4_Fnu']
			['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']
			['W1_mag_err','W2_mag_err','W3_mag_err','W4_mag_err']
			['W1_dex_err','W2_dex_err','W3_dex_err','W4_dex_err']
			['W1_snr','W2_snr','W3_snr','W4_snr']
			Additional Columns: 
			['WISE_sep','WISE_nsource']

	INPUT: table
		a table of objects with WISE inquery columns
	"""
	#==== setup 
	F_nu0=array([309.540,171.787,31.674,8.363])*u.Jy
	deltam=array([2.699,3.339,5.174,6.620])
	ws=array([3.4,4.6,12.,22.])*u.micron
	nus=lambda2nu(ws)

	#==== selection
	if type(selection)==bool:
		if selection ==False: selection=ones(len(table),dtype='bool')

	#==== creating columns
	colinfersuffixs=['_magVG','_mag_err','_dex_err','_snr','_magAB','_Fnu','_nuLnu_rf']
	coladditional=['WISE_sep','WISE_nsource']

	for suffix in colinfersuffixs:
		for BAND in ['W1','W2','W3','W4']:
			colname=BAND+suffix
			if colname not in table.dtype.names:
				print "creating column ", colname
				table[colname]=0.
				# table[colname].mask=True
				if suffix == '_Fnu': table[colname].unit="erg s-1 cm-2 Hz-1"
				if suffix == '_nuLnu_rf': table[colname].unit="erg s-1"
	for colname in coladditional:
		if colname not in table.dtype.names:
			print "creating column ", colname
			table[colname]=0.
			# table[colname].mask=True
			if colname == 'WISE_sep': table[colname].unit='arcsec'

	#==== Operation - writing WISE columns
	for i in arange(len(table))[selection]:
		if 'OBJID' in table.colnames:
			print i, table['OBJID'][i]
		else: 
			print i
		# do WISE query
		tabwise=queryWISE(table[i],catalog=catalog,searchradius=searchradius, updatelocal=updatelocal)
		if len(tabwise)!=1: raise ValueError("Length of tabwise not 1")
		# write the content
		for j in range(4):
			# name the bands
			band, BAND='w'+str(j+1), 'W'+str(j+1)
			# writing table
			table[BAND+'_magVG'][i]		=tabwise[band+'mpro']	
			table[BAND+'_mag_err'][i]	=tabwise[band+'sigmpro']
			table[BAND+'_dex_err'][i]	=tabwise[band+'sigmpro']/2.5
			table[BAND+'_snr'][i]		=tabwise[band+'snr']
			table[BAND+'_magAB'][i]		=table[BAND+'_magVG'][i]+deltam[j]
			table[BAND+'_Fnu'][i]		=(F_nu0[j]*10**(-table[BAND+'_magVG'][i]/2.5)).cgs.value
			table[BAND+'_nuLnu_rf'][i]	=conversion_util.L_from_F(nus[j]*table[BAND+'_Fnu'].quantity[i],table['Z'][i]).to('erg/s').value
		table['WISE_sep'][i] = tabwise['WISE_sep']
		table['WISE_nsource'][i]=tabwise['WISE_nsource']



def queryWISE(row, catalog='ALLWISE', searchradius=5., forcequery=False, updatelocal=True, verbos=True):
	"""
	PURPOSE: given a row with RA, DEC, return the WISE quiry either from local repository or by astroquery

	INPUT:  row : a row with RA DEC
			catalog: which wise catalog to use
				wise_allwise_p3as_psd: 		ALLWISE data release Nov. 2013 (default)
				wise_allsky_4band_p3as_psd: ALLSKY data release March 2012
			searchradius: radius of search in unit of arcsec
			forcequery=False: when True, force to do astroquery
			updatelocal=True: when True, update the local repository after doing the astroquery

	RETURN OUTPUT: 

	"""
	# sanity check
	if 'RA' not in row.columns or 'DEC' not in row.columns: raise ValueError("Row does not contain RA, DEC")

	# setup
	catalogkey={'ALLSKY':'allsky_4band_p3as_psd', 'ALLWISE':'allwise_p3as_psd'}
	filelocal='/Users/aisun/Documents/astro/projects/feedback/sample/algorithm/catalogue/wisequery_local/'+catalog+'.fits'
	# setup for object
	c = SkyCoord(row['RA'], row['DEC'], 'icrs', unit='deg')

	#== run local query
	runWISEquery=False
	if os.path.isfile(filelocal) and not forcequery: 
		if verbos: print "Running local  WISE quiry ", '%.2f'%row['RA'], '%.2f'%row['DEC']
		tablocal=Table.read(filelocal)		
		# sanity check of local table
		if len(unique(tablocal))!=len(tablocal): raise ValueError("Local WISE table contains duplicated rows")
		# match row with local table
		tabwise=catalogue_util.selectRADEC(tablocal,[row],searchradius,verbos=False)
		if len(tabwise)==0: runWISEquery=True
	else: runWISEquery=True

	#== run external query if local query failed
	if runWISEquery or forcequery:
		if verbos: print "Running remote WISE quiry ", '%.2f'%row['RA'], '%.2f'%row['DEC']
		#== WISE query
		tabwise=Irsa.query_region(c,catalog=catalogkey[catalog], spatial='Cone',radius=searchradius*u.arcsec)
		#== clean up 
		# del columns with dtype='object', which cannot be stored in fits
		tabwise['objid']=tabwise['designation'].astype('str')
		for column in tabwise.columns:
			if tabwise[column].dtype==dtype('O'):
				del tabwise[column]
		# rename ra, dec to capitalized
		tabwise.rename_column('ra','RA')
		tabwise.rename_column('dec','DEC')

		if updatelocal and len(tabwise)>=1:
			tabwise['WISE_querydate']=[datetime.date.today().isoformat()]
			if not os.path.isfile(filelocal): # if local file does not exist then create it with current tabwise
				tabwise.write(filelocal,format='fits',overwrite=True)
			else: # if local file does exist then append the new rows
				newtablocal=tablocal.copy()
				for i in range(len(tabwise)):
					# if there is no current row in the local table
					if not any((logical_and(newtablocal['RA']==tabwise[i]['RA'],newtablocal['DEC']==tabwise[i]['DEC']))):
						newtablocal.add_row(tabwise[i])
				# if there is new thing added then save it
				if len(newtablocal)>len(tablocal): newtablocal.write(filelocal,format='fits',overwrite=True)

	#== processing tabwise for output - deal with duplication and add on info of nsource and sep
	if len(tabwise)==0:
		# in tabwise create a masked row with only info of RA, DEC
		if verbos: print "NO MATCH"
		tabwise.add_row()
		for col in tabwise.colnames:
			tabwise[col].mask=True
		tabwise['RA']=row['RA']
		tabwise['DEC']=row['DEC']
		# create tabout and attach WISE_nsource and WISE_sep to it
		tabout=tabwise.copy()
		tabout['WISE_nsource']=[0]
		tabout['WISE_sep']=[nan]
	elif len(tabwise)==1:
		if verbos: print "MATCHED"		
		# create tabout and attach WISE_nsource and WISE_sep to it
		tabout=tabwise.copy()
		tabout['WISE_nsource']=[1]
		tabout['WISE_sep'] = c.separation(SkyCoord(tabwise['RA'], tabwise['DEC'], 'icrs', unit='deg')).arcsec[0]
	elif len(tabwise)>1:
		if verbos: print "CONFUSION"		
		# find the closest object
		distances=zeros(len(tabwise))
		for i in range(len(tabwise)):
			distances[i]=c.separation(SkyCoord(tabwise['RA'][i], tabwise['DEC'][i], 'icrs', unit='deg')).arcsec
		i_closest=argmin(distances)
		# create tabout to be the closest object and attach WISE_nsource and WISE_sep to it
		tabout=Table(tabwise[i_closest])
		tabout['WISE_nsource']=[len(tabwise)]
		tabout['WISE_sep']= [amin(distances)]

	return tabout



def correctWISE(table,selection=False):
	# correct the WISE MagAB to real magAB and F_nu
	if type(selection)==bool:
		if selection ==False: selection=ones(len(table))

	WISE_colnames_Mag=['W1_LUM_MagAB','W2_LUM_MagAB','W3_LUM_MagAB','W4_LUM_MagAB']
	WISE_ncolnames=['W1_Fnu','W2_Fnu','W3_Fnu','W4_Fnu','W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']
	WISE_colnames_Fnu=['W1_Fnu','W2_Fnu','W3_Fnu','W4_Fnu']
	WISE_colnames_Lnu=['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']
	for WISE_colname in WISE_ncolnames:
		if WISE_colname not in table.dtype.names:
			print "creating column ", WISE_colname
			table[WISE_colname]=0.
		if WISE_colname in  WISE_colnames_Fnu:
			table[WISE_colname].unit="erg /cm2 /Hz/s"
		if WISE_colname in WISE_colnames_Fnu:
			table[WISE_colname].unit="erg /s"

	ws=array([3.4,4.6,12.,22.])*u.micron
	nus=lambda2nu(ws)
	magABs=zeros(4)
	Fnus=zeros(4)*u.erg/u.cm**2/u.s/u.Hz
	for i in range(len(table)):
		if selection[i]:
			lum_dist=cosmo.luminosity_distance(table['Z'][i]).to(u.pc)
			dist_mod=-5.*log10(lum_dist/(10*u.pc)).value
			for j in range(len(WISE_colnames_Mag)):
				magABs[j]=table[WISE_colnames_Mag[j]][i]-dist_mod

				Fnus[j]=conversion_util.magAB2Fnu(magABs[j])
				table[WISE_colnames_Fnu[j]][i]=Fnus[j].value
				table[WISE_colnames_Lnu[j]][i]=(nus[j]*Fnus[j]*4.*pi*lum_dist**2).to('erg/s').value




def matchIRAS(table,selection=False):
	#==== selection
	if type(selection)==bool:
		if selection ==False: selection=ones(len(table))
	colnames=['IRAS_Fnu12', 'IRAS_Fnu25', 'IRAS_Fnu60', 'IRAS_Fnu100', 'IRAS_CATALOGUE','IRAS_NSOURCE']
	colinits=[0.,0.,0.,0.,'               ',0]
	colunits=['Jy','Jy','Jy','Jy',None,None]
	for i in range(len(colnames)):
		colname=colnames[i]
		if colname not in table.dtype.names:
			print "creating column ", colname
			table[colname]=colinits[i]
			table[colname].unit=colunits[i]
	for i in range(len(table)):
		if selection[i]:
			c = SkyCoord(table['RA'][i], table['DEC'][i], 'icrs', unit='deg')
			tabquery=Irsa.query_region(c,catalog='iraspsc', spatial='Cone',radius=1*u.arcmin)[['ra','dec','fnu_12','fnu_25','fnu_60','fnu_100']]
			if len(tabquery)==1:
				table['IRAS_Fnu12'][i]=tabquery['fnu_12']
				table['IRAS_Fnu25'][i]=tabquery['fnu_25']
				table['IRAS_Fnu60'][i]=tabquery['fnu_60']
				table['IRAS_Fnu100'][i]=tabquery['fnu_100']
				table['IRAS_CATALOGUE'][i]='point source'
				table['IRAS_NSOURCE'][i]=1
				print i, table['OBJID'][i], " psc match"
			elif len(tabquery)>2:
				table['IRAS_Fnu12'][i]=tabquery['fnu_12'][0]
				table['IRAS_Fnu25'][i]=tabquery['fnu_25'][0]
				table['IRAS_Fnu60'][i]=tabquery['fnu_60'][0]
				table['IRAS_Fnu100'][i]=tabquery['fnu_100'][0]
				table['IRAS_CATALOGUE'][i]='point source'
				table['IRAS_NSOURCE'][i]=len(tabquery)
				print i, table['OBJID'][i], " psc confusion"
			elif len(tabquery)==0:
				tabquery=Irsa.query_region(c,catalog='irasfsc', spatial='Cone',radius=1*u.arcmin)[['ra','dec','fnu_12','fnu_25','fnu_60','fnu_100']]
				if len(tabquery)==1:
					table['IRAS_Fnu12'][i]=tabquery['fnu_12']
					table['IRAS_Fnu25'][i]=tabquery['fnu_25']
					table['IRAS_Fnu60'][i]=tabquery['fnu_60']
					table['IRAS_Fnu100'][i]=tabquery['fnu_100']
					table['IRAS_CATALOGUE'][i]='faint source'
					table['IRAS_NSOURCE'][i]=1
					print i, table['OBJID'][i], " fsc match"
				elif len(tabquery)>2:
					table['IRAS_Fnu12'][i]=tabquery['fnu_12'][0]
					table['IRAS_Fnu25'][i]=tabquery['fnu_25'][0]
					table['IRAS_Fnu60'][i]=tabquery['fnu_60'][0]
					table['IRAS_Fnu100'][i]=tabquery['fnu_100'][0]
					table['IRAS_CATALOGUE'][i]='faint source'
					table['IRAS_NSOURCE'][i]=len(tabquery)
					print i, table['OBJID'][i], " fsc confusion"
				elif len(tabquery)==0:
					print i, table['OBJID'][i], " no match"




def matchNEDNAME(table,selection=False):
	colname='NEDNAME'
	if type(selection)==bool:
		if selection ==False: selection=ones(len(table))
	if colname not in table.dtype.names:
		table[colname]='                              '

	for i in range(len(table)):
		if selection[i]:
			c = SkyCoord(table['RA'][i], table['DEC'][i], 'icrs', unit='deg')

			tabNED = Ned.query_region(c,radius=2.*u.arcsec)

			if len(tabNED)==1:
				table[colname][i]=tabNED['Object Name'][0]
				print i, table['OBJID'][i], " matched"
			if len(tabNED)>=2:
				table[colname][i]=tabNED['Object Name'][0]
				print i, table['OBJID'][i], " confusion"
			if len(tabNED)==0:
				print i, table['OBJID'][i], " no match"

#=========== manual operations 

def writeFIRSTsearch(filename,table):
	file=open(filename,'w')
	for i in range(len(table)):
		c = SkyCoord(ra=table[i]['RA'], dec=table[i]['DEC'], unit=(u.degree, u.degree))
		file.write('%02d'%c.ra.hms[0]+' '+'%02d'%c.ra.hms[1]+' '+'%02d'%c.ra.hms[2]
				   +'.'+str(c.ra.hms[2]*10%10)[0]+' ')
		file.write('%+03d'%c.dec.dms[0]+' '+'%02d'%abs(c.dec.dms[1])+' '+'%02d'%abs(c.dec.dms[2])+' '+str(table[i]['OBJID'])+'\n')
	file.close()


def matchHerschelData_manual(table,data):
	table['PACSDATA']=False
	table['SPIREDATA']=False
	for objid in unique(data[data['instrument']=='PACS']['objid_01']):
		table['PACSDATA'][table['OBJID']==objid]=True
		print objid, " PACS"
	for objid in unique(data[data['instrument']=='SPIRE']['objid_01']):
		table['SPIREDATA'][table['OBJID']==objid]=True
		print objid, " SPIRE"
