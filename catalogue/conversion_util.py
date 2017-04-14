# conversion_utils.py
# ALS 04/07/2015

"""
Contain functions for conversions from some values to other values. 
"""
from pylab import *
import astropy.units as u
from astropy.table import Table, Column
import astropy.constants as const

from scipy import optimize
from scipy.interpolate import UnivariateSpline

# import catalogue_util
import crossmatch_util

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

kmps = u.def_unit('kmps', u.km / u.s)
u.add_enabled_units([kmps])
u.kmps=kmps

CO10nurest=115.27*u.GHz

def inferALLinTable(table):
	"""
	PURPOSE: make all possible inference  in a table
	SYNTAX: convertALLinTable(table)
	OUTPUT: write columns in the table, including: 
		'SDSSNAME'
		 'OIII_5007T_LUM','OIII_5007AVG_FWHM','Lbol_OIII',
		 ('W1_Fnu','W2_Fnu','W3_Fnu','W4_Fnu','W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf') if requires conversions
		 'WISE_nuLnu_rf12um','WISE_nuLnu_rf22um','WISE_nuLnu_rf30um',
		 'Lbol_12um','Lbol_22um','Lbol_30um'
		 'LMIR'
		 'LFIR'
		 'LpCO10EST'
	"""


	#==== infer quantities
	# OIII
	if set(['OIII_5007_LUM','OIII_5007B_LUM','OIII_5007_FLUX','OIII_5007B_FLUX','OIII_5007_FWHM','OIII_5007B_FWHM']) < set(table.colnames) and \
		not (set(['OIII_5007T_LUM','OIII_5007AVG_FWHM']) < set(table.colnames)):
		inferOIII_5007T_LUM(table)
		inferOIII_5007AVG_FWHM(table)

	if 'OIII_5007T_LUM' in table.colnames:
		table['Lbol_OIII']=Lbol_from_LOIII(table['OIII_5007T_LUM'])

	# WISE
	if set(['W1_LUM_MagAB','W2_LUM_MagAB','W3_LUM_MagAB','W4_LUM_MagAB']) < set(table.colnames) and not (set(['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']) < set(table.colnames)):
		crossmatch_util.correctWISE(table)

	if set(['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']) < set(table.colnames):
		spline_w_rf=[8.,12.,15.,22.,30.]
		for w_rf in spline_w_rf:
			addWISE_Lum_rf_spline(table,w_rf=w_rf,splineorder=2,touseWband=array([False,True,True,True]))
		# Lbol
		table['Lbol_8um']=Lbol_from_nuLnu8um(table['WISE_nuLnu_rf8um'])
		table['Lbol_12um']=Lbol_from_nuLnu12um(table['WISE_nuLnu_rf12um'])
		table['Lbol_15um']=Lbol_from_nuLnu15um(table['WISE_nuLnu_rf15um'])
		table['Lbol_22um']=Lbol_from_nuLnu22um(table['WISE_nuLnu_rf22um'])
		table['Lbol_30um']=Lbol_from_nuLnu30um(table['WISE_nuLnu_rf30um'])
		# LMIR
		table['LMIR']=LMIR_from_nuLnu25um(table['WISE_nuLnu_rf22um'])

	# LIR
	if set(['IRAS_Fnu60','IRAS_Fnu100','Z']) < set(table.colnames):
		table['LFIR']=LFIRest_from_Fnu60Fnu100(table['IRAS_Fnu60'],table['IRAS_Fnu100'],cosmo.luminosity_distance(table['Z']))
		table['LpCO10EST']=LpCO10_from_LFIR(table['LFIR'])
		table['LCO10EST']=L_from_Lp(table['LpCO10EST'],CO10nurest)
		addcolSdvCOs(table)


	# # RADIO
	# for suffix in ['_N','_IDS', '_FLUX', '_LUM']:
	# 	try: table.rename_column('FIR'+suffix,'FIRST'+suffix)
	# 	except: pass

	

#======= luminosity and flux
def magAB2Fnu(magAB):
	"""
	Converting AB magnitude to Fnu [erg / (cm2 Hz s)]
	"""
	return 10.**(-(magAB+48.6)/2.5)*u.erg/u.cm**2/u.s/u.Hz

def F_from_L(L,z):
	"""
	Converting from intrinsic luminosity to observed flux given redshift. 
	"""
	Dl=cosmo.luminosity_distance(z)	
	return L/(4.*pi*Dl**2)

def L_from_F(F,z):
	"""
	Converting from observed flux to intrinsic luminosity given redshift. 
	"""
	Dl=cosmo.luminosity_distance(z)	
	return F*(4.*pi*Dl**2)

def L_from_Lp(Lp,nurest):
	"""
	PURPOSE: convert L' [K kmps pc2] to L [Lsun] given line frequency
	SYNTAX: L_from_Lp(Lp,nurest)
	INPUT: 
			Lp: L' in units of [K kmps pc2]
			nurest: rest fram frequency in units of GHz
	"""
	return 3.20e-11*u.Lsun*(nurest/u.GHz)**3*Lp/(u.K*u.km/u.s*u.pc**2)

# def Fnudv_from_F(F,nurest,z):
# 	Fnudv=F/(nurest/(1.+z)/const.c)
# 	return Fnudv.to(u.Jy*u.km/u.s)

def Fnu_from_F_and_dv(F,nurest,z,dv):
	"""
	Calculate observed flux density (Fnu) in [Jy] given observed flux (F), 
	given rest frequency of line (nurest), redshift (z), and line width (dv)
	"""
	dnuobs=nurest/(1.+z)*dv/const.c
	Fnu=F/dnuobs
	return Fnu.to(u.Jy)


#======= OIII

def inferOIII_5007T_LUM(table):
	"""
	infer and add a column of total luminosity of OIII 5007 line (OIII_5007T_LUM)
	in the ALPAKA table from the other two columns (OIII_5007_LUM) and (OIII_5007B_LUM)
	"""
	OIII_5007T_LUM=table['OIII_5007_LUM']+table['OIII_5007B_LUM']
	table.add_column(Column(name='OIII_5007T_LUM',data=OIII_5007T_LUM,unit=u.erg/u.s))

def inferOIII_5007AVG_FWHM(table):
	"""
	infer and add a column of luminosity weighted averaged OIII line width (OIII_5007AVG_FWHM) 
	in the ALPAKA table
	"""
	table['OIII_5007AVG_FWHM']=0.
	totalflux=table['OIII_5007_FLUX']+table['OIII_5007B_FLUX']
	mask=totalflux!=0
	table['OIII_5007AVG_FWHM'][mask]=((table['OIII_5007_FWHM'][mask]*table['OIII_5007_FLUX'][mask]/totalflux[mask])**2+(table['OIII_5007B_FWHM'][mask]*table['OIII_5007B_FLUX'][mask]/totalflux[mask])**2)**0.5

def Lbol_from_LOIII(LOIII):
	"""
	estimate Lbol from LOIII
	from Liu+09 eq. 1, with +/- 0.5 dex error
	"""
	return 10.**3.5*(LOIII/(3.846e+33*u.erg/u.s))**0.99*u.Lsun.to(u.erg/u.s)*u.erg/u.s

#====== MIR

def addWISE_Lum_rf_power(table,w_rf=8.,touseWband=array([False,True,True,True])):
	"""
	Interpolate/Extrapolate MIR rest frame nuLnu from WISE band-2, 3, 4.
	A power law is fitted to the three bands and interpolation is evaluated using the power law. 
	The result will be stored in a new column of the table named e.g 'WISE_nuLnu_rf8um'

	INPUT: 
	w_rf=8. desires rest frame wavelength in um at which nuLnu is evaluated via interpolation. 

	splineorder=3 
		order of spline, notice order of spline has to be smaller than number of points
	touseWband=array([False,True,True,True])
		whether to use one WISE band or not [W1, W2, W3, W4]

	"""
	# define new columns name
	WISE_colnames_Lnu=['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']
	colname='WISE_nuLnu_rf'+'%.0f'%w_rf+'um' # e.g. 'WISE_nuLnu_rf8um'

	# create new column if not existed
	if colname not in table.dtype.names:
		table[colname]=0.
		table[colname].unit=u.erg/u.s
	
	# define fit function - linear in log-log space (power law)
	fitfunc = lambda p, x: p[0] + p[1] * x
	errfunc = lambda p, x, y: (y - fitfunc(p, x))
	ws=array([3.4,4.6,12.,22.])
	nuLnus=zeros(4)

	# iterate through rows in table
	for i in range(len(table)):
		for j in range(4):
			nuLnus[j]=table[WISE_colnames_Lnu[j]][i]
		# deredshift and take log of the three bandds
		wsrf=ws/(1.+table[i]['Z'])
		logx=log10(wsrf[touseWband])
		logy=log10(nuLnus[touseWband])
		
		# fit power law
		pinit = [1.0, -1.0]
		out = optimize.leastsq(errfunc, pinit, args=(logx, logy), full_output=1)
		pfinal = out[0]

		table[colname][i]=10**fitfunc(pfinal, log10(w_rf))

def addWISE_Lum_rf_spline(table,w_rf=8.,splineorder=2,touseWband=array([False,True,True,True])):
	"""
	Interpolate/Extrapolate MIR rest frame nuLnu at given wavelgnths from specified WISE band-2, 3, 4 using specified order spline.
	The result will be stored in a new column of the table named e.g 'WISE_nuLnu_rf8um'
	
	INPUT: 
		w_rf=8. desires rest frame wavelength in um at which nuLnu is evaluated via interpolation. 

		splineorder=2
			order of spline, notice order of spline has to be smaller than number of points used. 
			e.g. if use W2, 3, 4 then splineorder has to be =<2. 
		touseWband=array([False,True,True,True])
			whether to use one WISE band or not [W1, W2, W3, W4]
	"""
	# define new columns name
	WISE_colnames_Lnu=['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']
	colname='WISE_nuLnu_rf'+'%.0f'%w_rf+'um' # e.g. 'WISE_nuLnu_rf8um'

	# create new column if not existed
	if colname not in table.dtype.names:
		table[colname]=0.
		table[colname].unit=u.erg/u.s
	
	# iterate through rows in table
	for i in range(len(table)):
		nuLnus=zeros(4)
		for j in range(4):
			nuLnus[j]=table[WISE_colnames_Lnu[j]][i]
		table[colname][i]=splineWISE_Lum_rf(nuLnus,z=table[i]['Z'],wrf_target=w_rf,splineorder=splineorder,touseWband=touseWband)


def splineWISE_Lum_rf(nuLnus,z,wrf_target=8.,splineorder=2,touseWband=array([False,True,True,True])):
	"""
	the spline kernel inside addWISE_Lum_rf_spline

	INPUT:
			nuLnus: the source intrinsic nuLnu observed at 4 WISE bands
			z     : source redshift
			wsrf  : the source rest frame wavelgnths corresponding to 4 WISE bands
			wrf_target=8.: the wavlength [micron] to interpolate
			splineorder=2
				order of spline, notice order of spline has to be smaller than number of points used. 
				e.g. if use W2, 3, 4 then splineorder has to be =<2. 
			touseWband=array([False,True,True,True])
				whether to use one WISE band or not [W1, W2, W3, W4]

	"""
	# setup
	ws=array([3.4,4.6,12.,22.])
	wsrf=ws/(1.+z)
	# spline
	logx=log10(wsrf[touseWband])
	logy=log10(nuLnus[touseWband])
	
	spl = UnivariateSpline(logx, logy, k=splineorder)
	return 10**spl(log10(wrf_target))


def Lbol_from_nuLnu8um(nuLnu):
	"""
	AGN bolometric luminosity from nuLnu at rest frame 8 micron
	with a bolometric correction of 9 from Richards et al. (2006), see also Liu et al. (2013b)
	"""
	return 8.8*nuLnu

def Lbol_from_nuLnu12um(nuLnu):
	"""
	AGN bolometric luminosity from nuLnu at rest frame 12 micron
	with a bolometric correction of 9 from Richards et al. (2006), see also Liu et al. (2013b)
	"""
	return 9.*nuLnu

def Lbol_from_nuLnu15um(nuLnu):
	"""
	AGN bolometric luminosity from nuLnu at rest frame 15 micron
	with a bolometric correction of 9 from Richards et al. (2006). 
	"""
	return 9.*nuLnu

def Lbol_from_nuLnu22um(nuLnu):
	"""
	AGN bolometric luminosity from nuLnu at rest frame 22 micron
	with a bolometric correction of 11 from Richards et al. (2006)
	"""
	return 11.*nuLnu

def Lbol_from_nuLnu30um(nuLnu):
	"""
	AGN bolometric luminosity from nuLnu at rest frame 30 micron
	with a bolometric correction of 13 from Richards et al. (2006), see also Liu et al. (2013b)
	"""
	return 13.*nuLnu

def LMIR_from_nuLnu25um(nuLnu):
	"""
	rough estimate of Mid-infrared Luminosity from nuLnu at rest frame 25 micron 
	from Krips+12 Table 3. originally from Dale & Helou (2002)
	WARNING: it's a rough esimtate for z ~0.1
	"""
	return 1.75*nuLnu


#======== FIR

def LFIRest_from_Fnu60Fnu100(Fnu60,Fnu100,Dl):
	"""
	estimate LFIR using IRAS Fnu60, Fnu100
	from Solomon 97 eq. 4, assuming r(S60/S100) = 1.8
	"""
	return 3.94e5*u.Lsun*(2.58*(Fnu60/u.Jy)+(Fnu100/u.Jy))*1.8*(Dl/u.Mpc)**2

def LpCO10_from_LFIR(LFIR):
	"""
	DESCRIPTION:
		estimate L'CO10 from LFIR 
		using Villar-Martin+13 Fig. 4 top panel regression line
	"""
	return 80.*(u.K*u.kmps*u.pc**2)*(LFIR/u.Lsun)**0.65


#======== molecular


def Mmol_from_LpCO10(LpCO10, sourcetype='ULIRG'):
	"""
	DESCRIPTION: 
		estimate mass of molecular from L'CO(1-0)
		assuming alpha_CO value in Carilli and Walter 2013

	PARAMS
	---------
	LpCO10: quantity
		Lp in unit of (u.K*u.kmps*u.pc**2)
	sourcetype: 'string'
		either 'Milky Way' or 'ULIRG'
	"""
	if sourcetype == 'Milky Way':
		alpha_CO = 4. *u.M_sun/(u.K*u.km/u.s*u.pc**2)

	elif sourcetype == 'ULIRG': 
		alpha_CO = 0.8 *u.M_sun/(u.K*u.km/u.s*u.pc**2)
	else: 
		raise NameError('source type not recognized')

	return (LpCO10 * alpha_CO).to(u.M_sun)


def Mmdense_from_LpHCN10(LpHCN10):
	"""
	DESCRIPTION: 
		estimate mass of dense molecular from L'HCN(1-0)
		assuming alpha_HCN value in Gao & Solomon 2004a
	"""
	return u.Msun*10.*(LpHCN10/(u.K*u.kmps*u.pc**2))


def Lp_from_Sdv(Sdv,nurest,z):
	"""
	PURPOSE: 
		Calculate molecular line luminosity Lp [K km s-1 pc2] from flux [Jy km s-1],
		given redshift z and line frequency nurest

	DESCRIPTION: 
		Following Carilli & Walter 2013 Section 2.4, but presenting frequency in nurest:

			Lp/ (K kmps pc^2) = 3.25e7 * Sdv/(Jy kmps) * (1+z)^-1 * (nurest/GHz)^-2 * (Dl/Mpc)^2

		note: Luminosity distance will be automatically calculated from z

	INPUT: 
			Sdv [Jy km s-1] : line flux
			nurest [GHz]    : line rest frequency
			z   []          : redshift

	OUTPUT: 
			Lp  [K kmps pc^2]: Luminosity
	"""
	Dl=cosmo.luminosity_distance(z)
	Lp=(u.K*u.kmps*u.pc**2) * 3.25e7 * Sdv/(u.Jy*u.kmps) * (1.+z)**-1 * (nurest/u.GHz)**-2 * (Dl/u.Mpc)**2
	return Lp


def Sdv_from_Lp(Lp,nurest,z):
	"""
	PURPOSE: 
		Calculate molecular line flux [Jy km s-1] from luminosity Lp [K km s-1 pc2]
		given redshift z and line frequency nurest

	DESCRIPTION: 
		Following Carilli & Walter 2013 Section 2.4, but presenting frequency in nurest:

			Sdv/(Jy kmps) = 3.077e-8 * Lp/(K kmps pc^2) * (1+z) * (nurest/GHz)^2 * (Dl/Mpc)^-2

		note: Luminosity distance will be automatically calculated from z

	INPUT: 
			Lp  [K kmps pc^2]: Luminosity
			nurest [GHz]    : line rest frequency
			z   []          : redshift

	OUTPUT: 
			Sdv [Jy km s-1] : line flux
	"""
	Dl=cosmo.luminosity_distance(z)
	Sdv= (u.Jy*u.kmps) * 3.077e-8 * Lp/(u.K*u.kmps*u.pc**2) * (1.+z) * (nurest/u.GHz)**2 * (Dl/u.Mpc)**-2
	return Sdv


def getnu_molline(line='CO10'):
	"""
	PURPOSE: 
		return the rest frame ferquency in GHz of the specified molecular lines

	DESCRIPTION:
		available lines are: CO10, CO21, CO32, CO43 
	"""
	linenurest={'CO10':115.2712018*u.GHz,'CO21':230.538000*u.GHz,'CO32':345.7959899*u.GHz,'CO43':461.0407682*u.GHz,}

	return linenurest[line]



def addcolSdvCOs(table, sourcetype='SMG'):
	"""
	PURPOSE: add new column to table with estimated Sdv of lines CO10, CO21, CO32, CO43. 
			 given columns "LpCO10EST", "Z"

	DESCRIPTION: 
			 all quantities are inferred from LpCO10EST column of the table
			 Sdv is inferred from Lp using function conversion_util.Sdv_from_Lp(), which is using Carilli & Walter 13 formula. 

			 Higher excitation lines are estimated from CO10 using specified line ratios (SMG or QSO):
				 SMG: {"rCO21CO10":0.85,"rCO32CO10":0.66,"rCO43CO10":0.46,}		 
				 QSO: {"rCO21CO10":0.99,"rCO32CO10":0.97,"rCO43CO10":0.87,}

			 Note: rCO32CO10 of SMG are similar to that of J1356, which can be used as a good guess for the rest of the bubbles

	INPUT:   table (astropy Table): 
				a table containing columns "LpCO10EST", "Z"

			 sourcetype (string)='SMG' : 
			 	a string specifying which typical line ratio from Carilli & Walter 13
			 	options are 'SMG', 'QSO'


	OUTPUT: 
			 new columns in table: 
			 "SdvCO10EST", "SdvCO21EST", "SdvCO32EST", "SdvCO43EST"
	"""
	#=== set up
	LpCO10=table['LpCO10EST']
	z=table['Z']

	#=== calculate Sdv for lines
	nurest = getnu_molline('CO10')
	table['SdvCO10EST'] = Sdv_from_Lp(LpCO10,nurest,z)

	nurest = getnu_molline('CO21')
	lineratio = getCO_lineratio(ratiotype = "rCO21CO10", sourcetype=sourcetype)
	table['SdvCO21EST'] = Sdv_from_Lp(LpCO10*lineratio, nurest, z)

	nurest = getnu_molline('CO32')
	lineratio = getCO_lineratio(ratiotype = "rCO32CO10", sourcetype=sourcetype)
	table['SdvCO32EST'] = Sdv_from_Lp(LpCO10*lineratio, nurest, z)

	nurest = getnu_molline('CO43')
	lineratio = getCO_lineratio(ratiotype = "rCO43CO10", sourcetype=sourcetype)
	table['SdvCO43EST'] = Sdv_from_Lp(LpCO10*lineratio, nurest, z)


def getCO_lineratio(ratiotype = "rCO21CO10", sourcetype='SMG'): 
	""" return CO line ratio of specific transitions and source type"""
	# Assumption of line ratios
	if sourcetype == 'SMG':
		 # typical SMG line ratio (Carilli & Walter 13), which has rCO32CO10 similar to J1356 
		print "Assuming typical SMG CO line ratios from Carilli & Walter 13"
		lineratio = {"rCO21CO10":0.85,"rCO32CO10":0.66,"rCO43CO10":0.46,}
	elif sourcetype == 'QSO':
		 # typical QSO line ratio (Carilli & Walter 13)
		print "Assuming typical QSO CO line ratios from Carilli & Walter 13"
		lineratio = {"rCO21CO10":0.99,"rCO32CO10":0.97,"rCO43CO10":0.87,}
	return lineratio[ratiotype]

def calc_nbeam(dtheta, beam_maj, beam_min): 
	"""
	calculate the number of beams in the area of the source. assuming the source flux is distributed in a square region of diameter dtheta and area dtheta^2 observed by a beam size of area 1.13*beam_maj*beam_min. 

	"""
	nbeam = max(1, dtheta**2/(1.13 * beam_maj * beam_min))
	return nbeam


def calcSNR(Sdv, dtheta, dv, S_obs_RMS_perchanbeam, beam_maj, beam_min, channelw):
	"""
	PURPOSE: calculate the SNR per channel per beam and total SNR given source and observation parameters

	SYNTAX:  calcSNR(Sdv, dtheta, dv, S_obs_RMS_perchanbeam, channel, beam_maj, beam_min)

	INPUT:  
			Sdv [Jy km/s]: source flux 
			dtheta  [arcsec] : source angular size
			dv  [km/s]   : source line width
			S_obs_RMS_perchanbeam [Jy]     : observation RMS flux density per channel per beam
			beam_maj[arcsec]: observation beam size major axis
			beam_min[arcsec]: observation beam size minor axis
			channelw	[km/s]	: observation channel width

	DESCRIPTION:

			Here I assume the source flux is distributed in a square region of diameter dtheta and area dtheta^2, 
			and spectrally top hat line shape of width dv, observed by a beam size of area
			1.13*beam_maj*beam_min and a channel of width dv.

	RETURN OUTPUT: 

			SNR_perchanbeam: SNR per channel per beam
			SNR_toal: total SNR
	"""

	# calculate the number of beams in the source area
	nbeam = calc_nbeam(dtheta, beam_maj, beam_min)
	nchannel = dv/channelw

	# flux density in a beam area
	S_src_perbeam = Sdv/dv/nbeam

	# SNR per beam per channel is the ratio between S and S_obs_RMS_perchanbeam in one beam and channel
	SNR_perchanbeam = (S_src_perbeam/S_obs_RMS_perchanbeam).to(u.dimensionless_unscaled)

	# SNR total is SNR per beam per channel times sqrt(# of beam*channel)
	SNR_total = SNR_perchanbeam * np.sqrt(nbeam) * np.sqrt(nchannel)

	return SNR_perchanbeam, SNR_total


def calcSdv_RMS_sourcetotal(dtheta, dv, S_obs_RMS_perchanbeam, beam_maj, beam_min, channelw): 
	"""
	calculate the RMS of the source Sdv of a given the observational RMS of S per channel per beam (S_obs_RMS_perchanbeam). The source diameter (dtheta), line width (dv), and the observation beam size (beam_maj, beam_min), and channel width (channelw) need to be given. 

	PARAMS
	----------
	dtheta, dv, S_obs_RMS_perchanbeam, beam_maj, beam_min, channelw
	see calcSNR
	"""

	# calculate the number of beams in the source area
	nbeam = calc_nbeam(dtheta, beam_maj, beam_min)
	nchannel = dv/channelw

	S_RMS_total = S_obs_RMS_perchanbeam / np.sqrt(nchannel*nbeam)
	Sdv_RMS_total = S_RMS_total*dv
	return Sdv_RMS_total


def calcM_mol_RMS(dtheta, dv, S_obs_RMS_perchanbeam, beam_maj, beam_min, channelw, z, line='CO32', sourcetype_alphaCO='ULIRG', sourcetype_COlratio='QSO'): 
	""" 
	calculate the molecular mass RMS detection limit of a given source and observation setup 
	"""

	Sdv_RMS_sourcetotal = calcSdv_RMS_sourcetotal(dtheta, dv, S_obs_RMS_perchanbeam, beam_maj, beam_min, channelw)

	# get source Lp of the observed line
	nurest = getnu_molline(line=line)
	Lp = Lp_from_Sdv(Sdv_RMS_sourcetotal, nurest, z)

	# infer Lp of CO10 assuming line ratio
	lineratio = getCO_lineratio(ratiotype = "r"+line+"CO10", sourcetype=sourcetype_COlratio)
	LpCO10 = Lp/lineratio

	# infer M_mol assuming alphaCO
	M_mol_RMS = Mmol_from_LpCO10(LpCO10, sourcetype=sourcetype_alphaCO)

	return M_mol_RMS

#======== SFR

def SFR_from_LnuFF(Lnu, nu, Te):
	"""
	PURPOSE: 
		Estimate SFR from free-free emission in the radio
	DESCRIPTION:
		Equation (6) in Murphy+2012, http://adsabs.harvard.edu/abs/2012ApJ...761...97M
		Used by Alatalo+12 for SFR of N1266, http://adsabs.harvard.edu/abs/2015ApJ...798...31A
	SYNTAX: SFR_from_FF(Lnu, nu, Te)
	"""

	return u.Msun/u.yr*4.6e-28*(Te/1.e4*u.K)**-0.45*(nu/u.GHz)**0.1*(Lnu/(u.erg/u.s/u.Hz))


def SFR_from_LIR(LIR):
	"""
	DESCRIPTION: convert  LIR to SFR 
	using Kennicutt+98
	"""
	return u.Msun/u.yr*1.7e-10*(LIR/u.Lsun)
