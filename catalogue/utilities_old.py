# handy tools
from pylab import *

from astropy.table import Table, hstack, Column

from astroquery.irsa import Irsa
from astroquery.ned import Ned
from astropy.coordinates import SkyCoord
import astropy.units as u

from astropy import constants as const
# from astropy.cosmology import luminosity_distance
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

import os

from astropy import constants as const

from scipy import optimize


kms=u.km/u.s

def SDSSName(c):
	# WARNING! not tested
	return "SDSSJ"+"%02d"%c.ra.hms[0]+"%02d"%c.ra.hms[1]+("%019.16f"%c.ra.hms[2])[0:5]+"%+03d"%c.dec.dms[0]+"%02d"%c.dec.dms[1]+("%019.16f"%3.49)[0:4]


def getSDSSName(c):
	# return "SDSSJ"+"%02d"%c.ra.hms[0]+"%02d"%c.ra.hms[1]+("%019.16f"%c.ra.hms[2])[0:5]+"%+03d"%c.dec.dms[0]+"%02d"%c.dec.dms[1]+("%019.16f"%c.dec.dms[2])[0:4]
	return "SDSSJ"+"%02d"%c.ra.hms[0]+"%02d"%c.ra.hms[1]+"%+03d"%c.dec.dms[0]+"%02d"%abs(c.dec.dms[1])

def addcolSDSSName(table):
	colname='SDSSNAME'
	if colname not in table.dtype.names:
		print "creating column ", colname
		table[colname]='                    '
	for i in range(len(table)):
		c = SkyCoord(ra=table[i]['RA'], dec=table[i]['DEC'], unit=(u.degree, u.degree))
		table[colname][i]=getSDSSName(c)


def gethhmmsss(c,delimiter=':'):
	return "%02d"%c.ra.hms[0]+delimiter+"%02d"%c.ra.hms[1]+delimiter+("%019.16f"%c.ra.hms[2])[0:4]

def getsddmmss(c,delimiter=':'):
	if c.dec >=0:
		sign='+'
	else:
		sign='-'
	return sign+"%02d"%c.dec.dms[0]+delimiter+"%02d"%abs(c.dec.dms[1])+delimiter+"%02d"%abs(round(c.dec.dms[2]))

def addcolRAhmsDECdms(table):
	if 'RA_hms' not in table.dtype.names:
		table['RA_hms']='               '
	if 'DEC_dms' not in table.dtype.names:
		table['DEC_dms']='               '
	for i in range(len(table)):
		c = SkyCoord(ra=table[i]['RA'], dec=table[i]['DEC'], unit=(u.degree, u.degree))
		table[i]['RA_hms']=gethhmmsss(c)
		table[i]['DEC_dms']=getsddmmss(c)

#=========== write JSky object list
def writeJSkylist(filename,table):
	file=open(filename,'w')
	if 'Name' in table.dtype.names:
		name=table['Name']
	elif 'name' in table.dtype.names:
		name=table['name']
	elif 'NAME'  in table.dtype.names:
		name=table['NAME']
	else:
		name=arange(len(table))
	for i in range(len(table)):
		file.write(str(name[i])+' ')
		c = SkyCoord(ra=table[i]['RA'], dec=table[i]['DEC'], unit=(u.degree, u.degree))
		file.write('%02d'%c.ra.hms[0]+' '+'%02d'%c.ra.hms[1]+' '+'%02d'%c.ra.hms[2]
				   +'.'+str(c.ra.hms[2]*10%10)[0]+' ')
		file.write('%+03d'%c.dec.dms[0]+' '+'%02d'%abs(c.dec.dms[1])+' '+'%02d'%abs(c.dec.dms[2])+' ')
		file.write('2000 \n')
	file.close()


#=========== write FIRST search list
def writeFIRSTsearch(filename,table):
	file=open(filename,'w')
	for i in range(len(table)):
		c = SkyCoord(ra=table[i]['RA'], dec=table[i]['DEC'], unit=(u.degree, u.degree))
		file.write('%02d'%c.ra.hms[0]+' '+'%02d'%c.ra.hms[1]+' '+'%02d'%c.ra.hms[2]
				   +'.'+str(c.ra.hms[2]*10%10)[0]+' ')
		file.write('%+03d'%c.dec.dms[0]+' '+'%02d'%abs(c.dec.dms[1])+' '+'%02d'%abs(c.dec.dms[2])+' '+str(table[i]['OBJID'])+'\n')
	file.close()



#=========== write SDSS Image Search
# http://skyserver.sdss3.org/public/en/tools/chart/list.aspx
def writeSDSSSImgearch(filename,table,header=''):
	if 'Name' in table.dtype.names:
		name=table['Name']
		savetxt(filename, array([name,table['RA'],table['DEC']]).transpose(),delimiter=' ', comments='#', fmt=['%s','%.3f','%.3f'],header=header)
	if 'name' in table.dtype.names:
		name=table['name']
		savetxt(filename, array([name,table['RA'],table['DEC']]).transpose(),delimiter=' ', comments='#', fmt=['%s','%.3f','%.3f'],header=header)
	else:
		name=arange(len(table)).astype(int)
		savetxt(filename, array([name,table['RA'],table['DEC']]).transpose(),delimiter=' ', comments='#', fmt=['%i','%.3f','%.3f'],header=header)


#=========== write finding chart wget script
def writeSDSSfchartwget(filename,table):
	# require:  table['RA'], table['DEC']
	# optional: table['Name'] or table['name']
	if 'Name' in table.dtype.names:
		name=table['Name']
	if 'name' in table.dtype.names:
		name=table['name']
	if 'NAME' in table.dtype.names:
		name=table['NAME']
	else:
		name=arange(len(table)).astype('str')
	file=open(filename,'w')
	for i in range(len(table)):
		file.write('wget -O '+str(name[i])+'.jpg \'http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx?ra='+'%.5f'%table[i]['RA']+'&dec='+'%.5f'%table[i]['DEC']+'&scale=0.4&width=850&height=850&opt=GLI&query=R(10,20)\'\n')
	file.close()

def makefchart(directory,table):
	mycwd = os.getcwd()
	os.system('mkdir '+directory+'findcharts')
	os.chdir(directory+'findcharts/')
	writeSDSSfchartwget('wgetscript.txt',table)
	os.system('chmod 755 wgetscript.txt')
	os.system('./wgetscript.txt')
	os.chdir(mycwd)


#=========== write stamp wget script
def writeSDSSstampwget(filename,table):
	# require:  table['RA'], table['DEC']
	file=open(filename,'w')
	if 'Name' in table.dtype.names:
		name=table['Name']
	elif 'name' in table.dtype.names:
		name=table['name']
	elif 'NAME' in table.dtype.names:
		name=table['NAME']
	else:
		name=arange(len(table)).astype('str')
	for i in range(len(table)):
		#file.write('wget -O '+str(i)+'_'+'%.5f'%table['RA'][i]+'_'+'%.5f'%table['DEC'][i]+'.jpg \'http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx?ra='+'%.5f'%table[i]['RA']+'&dec='+'%.5f'%table[i]['DEC']+'&scale=0.1&width=600&height=600&opt=GL&query=R(10,20)\'\n')
		# zoomed in
		file.write('wget -O stamp_'+name[i]+'.jpg \'http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx?ra='+'%.5f'%table[i]['RA']+'&dec='+'%.5f'%table[i]['DEC']+'&scale=0.1&width=300&height=300&opt=G&query=R(10,20)\'\n')
	file.close()

def makestamp(directory,table):
	mycwd = os.getcwd()
	os.system('mkdir '+directory+'stamps')
	os.chdir(directory+'stamps/')
	writeSDSSstampwget('wgetscript.txt',table)
	os.system('chmod 755 wgetscript.txt')
	os.system('./wgetscript.txt')
	os.chdir(mycwd)


#============ Magellan catalogue
def writeMagellanCatalogue(filename,table):
	if any(table.dtype.names=='Name'): name=table['Name']
	elif any(table.dtype.names=='name'): name=table['name']
	else: name=arange(table.size).astype('str')

	if any(table.dtype.names=='PA'): PA=table['PA']
	else: PA=zeros(table.size)
	
	file=open(filename,'w')
	file.write('#                   RA         Dec       equinox RApm Decpm offset rot  RA_probe1  Dec_probe1 equinox RA_probe2  Dec_probe2 equinox pm_epoch\n### name            hh:mm:ss.s sdd:mm:ss yyyy.0  s.ss s.ss  angle  mode hh:mm:ss.s sdd:mm:ss  yyyy.0  hh:mm:ss.s sdd:mm:ss  yyyy.0  yyyy.0 #\n')

	for i in range(table.size):
		file.write(name[i]+' ')
		c = ICRS(ra= table[i]['RA'], dec= table[i]['DEC'], unit=(u.degree, u.degree))
		file.write(getSDSSName(c)+' ')
		file.write(gethhmmsss(c)+' '+getsddmmss(c)+' ')
		file.write('2000.0  0.00 0.00  '+'%3.2f'%((PA[i]+180.-46.15)%360))
		file.write(' EQU  00:00:00.0 +00:00:00  2000.0  00:00:00.0 +00:00:00  2000.0  2000.0 # \n')
	file.close()


#============ write APO Catalogue
def writeAPOCatalogue(filename,table,header="APO catalogue"):
	# table['PA'] in degree east of north
	# catalogue RotAng = PA - 90
	#
	file=open(filename,'w')
	file.write("# "+header+"\n")
	file.write("# Name	Pos1 (e.g. RA)	Pos2 (e.g. Dec)	Options in keyword-value format\n")

	for i in range(len(table)):
		file.write(str(table[i]['NAME'])+'	')
		c = SkyCoord(ra= table[i]['RA'], dec= table[i]['DEC'], unit=(u.degree, u.degree))
		file.write(gethhmmsss(c)+'	'+getsddmmss(c)+'	'+'RotAng='+str(table[i]['PA']-90)+'\n')
	file.close()


#============= write summary
def writeSummary(filename,table):
	if 'N' not in table.dtype.names:
		table['N']=arange(len(table))
	colnames=['N','OBJID','RA','DEC','Z','AGN_TYPE','OIII_5007T_LUM','OIII_5007AVG_FWHM','W4_LUM_MagAB']
	table[colnames].write(filename,format='ascii.fixed_width',delimiter='')




def matchWISE(table,selection=False):
	# table: a table with at least two columns 'RA', 'DEC'
	# selection: an array of bool entried with length equals to table
	#==== selection
	if type(selection)==bool:
		if selection ==False: selection=ones(len(table))
	#==== check if table has WISE columns, if not creat them
	WISE_colnames_Fnu=['W1_Fnu','W2_Fnu','W3_Fnu','W4_Fnu']
	WISE_colnames_Lnu=['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']
	for WISE_colname in WISE_colnames_Fnu+WISE_colnames_Lnu+['WISE_Sep','WISE_nsource']:
		if WISE_colname not in table.dtype.names:
			print "creating column ", WISE_colname
			table[WISE_colname]=0.
			if WISE_colname in  WISE_colnames_Fnu:
				table[WISE_colname].unit="erg /cm2 /Hz/s"
			if WISE_colname in WISE_colnames_Fnu:
				table[WISE_colname].unit="erg /s"
	table['WISE_Sep'].unit='arcsec'

	#==== set up
	query_colname=['w1mpro','w2mpro','w3mpro','w4mpro',]
	deltam=array([2.683,3.319,5.242,6.604])
	ws=array([3.4,4.6,12.,22.])*u.micron
	nus=lambda2nu(ws)

	#==== matching with ALLWISE
	for i in range(len(table)):
		if selection[i]:
			magABs=zeros(4)
			Fnus=zeros(4)*u.erg/u.cm**2/u.s/u.Hz
			# calculate distance modulo
			lum_dist=cosmo.luminosity_distance(table['Z'][i]).to(u.pc)
			dist_mod=-5.*log10(lum_dist/(10*u.pc)).value
			c = SkyCoord(table['RA'][i], table['DEC'][i], 'icrs', unit='deg')
			tabwise=Irsa.query_region(c,catalog='wise_allwise_p3as_psd', spatial='Cone',radius=5*u.arcsec)[['designation','ra','dec','w1mpro','w2mpro','w3mpro','w4mpro']]
			if len(tabwise)==1:
				# retrieving magAB
				for j in range(4):
					magABs[j]=tabwise[query_colname[j]]+deltam[j]
					Fnus[j]=magAB2Fnu(magABs[j])
					table[WISE_colnames_Fnu[j]][i]=Fnus[j].value
					table[WISE_colnames_Lnu[j]][i]=(nus[j]*Fnus[j]*4.*pi*lum_dist**2).to('erg/s').value

				# WISE_Sep: distance between SDSS and WISE coordinates
				cwise = SkyCoord(tabwise['ra'], tabwise['dec'], 'icrs', unit='deg')
				table['WISE_Sep'][i] = c.separation(cwise).arcsec
				table['WISE_nsource'][i]=1.
				print i, table['OBJID'][i], " matched"
			elif len(tabwise)>=2:
				# calculate distance modulo
				# determining the closest object
				cwise1 = SkyCoord(tabwise['ra'][0], tabwise['dec'][0], 'icrs', unit='deg')
				cwise2 = SkyCoord(tabwise['ra'][1], tabwise['dec'][1], 'icrs', unit='deg')
				dist1=c.separation(cwise1).arcsec
				dist2=c.separation(cwise2).arcsec
				k= dist1>dist2
				for j in range(4):
					magABs[j]=tabwise[query_colname[j]][k]+deltam[j]
					Fnus[j]=magAB2Fnu(magABs[j])
					table[WISE_colnames_Fnu[j]][i]=Fnus[j].value
					table[WISE_colnames_Lnu[j]][i]=(nus[j]*Fnus[j]*4.*pi*lum_dist**2).to('erg/s').value
				table['WISE_Sep'][i] = min(dist1,dist2)
				table['WISE_nsource'][i]=len(tabwise)
				print i, table['OBJID'][i], " confusion"
			elif len(tabwise)==0:
				table['WISE_Sep'][i]=-1
				table['WISE_nsource'][i]=0
				print i, table['OBJID'][i], " no match"
			else:
				print i, table['OBJID'][i], " exception"
	for WISE_colname in WISE_colnames_Fnu:
		table[WISE_colname].unit="erg /cm2 /Hz/s"
	for WISE_colname in WISE_colnames_Lnu:
		table[WISE_colname].unit="erg /s"



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

				Fnus[j]=magAB2Fnu(magABs[j])
				table[WISE_colnames_Fnu[j]][i]=Fnus[j].value
				table[WISE_colnames_Lnu[j]][i]=(nus[j]*Fnus[j]*4.*pi*lum_dist**2).to('erg/s').value


def addWISE_Lum_rf(table,w_rf=8.):
	# use WISE2, 3, 4 bands to fit an powerlow and to interpolate a resf frame nuLnu
	WISE_colnames_Lnu=['W1_nuLnu_rf','W2_nuLnu_rf','W3_nuLnu_rf','W4_nuLnu_rf']
	colname='WISE_nuLnu_rf'+'%.0f'%w_rf+'um' # e.g. 'WISE_nuLnu_rf8um'
	if colname not in table.dtype.names:
		table[colname]=0.
		table[colname].unit=u.erg/u.s
	
	fitfunc = lambda p, x: p[0] + p[1] * x
	errfunc = lambda p, x, y: (y - fitfunc(p, x))
	ws=array([3.4,4.6,12.,22.])
	nuLnus=zeros(4)
	for i in range(len(table)):
		for j in range(4):
			nuLnus[j]=table[WISE_colnames_Lnu[j]][i]
		wsrf=ws/(1.+table[i]['Z'])
		logx=log10(wsrf[1:])
		logy=log10(nuLnus[1:])
		
		pinit = [1.0, -1.0]
		out = optimize.leastsq(errfunc, pinit, args=(logx, logy), full_output=1)
		pfinal = out[0]

		table[colname][i]=10**fitfunc(pfinal, log10(w_rf))




def lambda2nu(wavelength):
	return (const.c/wavelength).to('Hz')


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


def matchHerschelData_manual(table,data):
	table['PACSDATA']=False
	table['SPIREDATA']=False
	for objid in unique(data[data['instrument']=='PACS']['objid_01']):
		table['PACSDATA'][table['OBJID']==objid]=True
		print objid, " PACS"
	for objid in unique(data[data['instrument']=='SPIRE']['objid_01']):
		table['SPIREDATA'][table['OBJID']==objid]=True
		print objid, " SPIRE"




#=============== match with NED
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


def matchSDSSNAME(table,selection=False):
	colname='NEDSDSSNAME'
	if type(selection)==bool:
		if selection ==False: selection=ones(len(table))
	if colname not in table.dtype.names:
		table[colname]='                              '
	
	for i in range(len(table)):
		if selection[i]:
			c = SkyCoord(table['RA'][i], table['DEC'][i], 'icrs', unit='deg')
			
			tabNED = Ned.query_region(c,radius=2.*u.arcsec)
			for j in range(len(tabNED)):
				if tabNED['Object Name'][j][:4]=='SDSS' and table[colname][i][:4]!='SDSS':
					table[colname][i]=tabNED['Object Name'][j]
					print i, table['OBJID'][i], " matched", tabNED['Object Name'][j]
				elif tabNED['Object Name'][j][:4]=='SDSS' and table[colname][i][:4]=='SDSS':
					print i, table['OBJID'][i], " confusion", tabNED['Object Name'][j]
			if table[colname][i][:4]!='SDSS':
					print i, table['OBJID'][i], " no match"


#=========== plotting

def plotzOIIIW4(table,marker='.',toclf=True,labelsuffix='',plotFWHM500=True):
	if toclf: clf()
	markersize=ceil((-table['W4_LUM_MagAB']-23)).astype(int)
	if plotFWHM500 and 'OIII_5007AVG_FWHM' in table.dtype.names:
		selFWHM500=table['OIII_5007AVG_FWHM']>500
		plot(table['Z'][selFWHM500],log10(table['OIII_5007T_LUM'][selFWHM500]),ls='',marker='+',ms=6,color='black',label='OIII FWHM >500'+labelsuffix)
	for i in range(max(markersize)+1):
		selection= markersize==i
		plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker=marker,ms=i,color='black')
		if i ==7:
			plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker=marker,ms=i,color='red',label='W4 = -29~-30'+labelsuffix)
			if plotFWHM500 and 'OIII_5007AVG_FWHM' in table.dtype.names:
				selection=all([selection,selFWHM500],axis=0)
				plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker='+',ms=i,color='red')
		if i ==6:
			plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker=marker,ms=i,color='blue',label='W4 = -28~-29'+labelsuffix)
		if i ==5:
			plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker=marker,ms=i,color='green',label='W4 = -27~-28'+labelsuffix)
	xlabel('z')
	ylabel('log10 LOIII')
	legend(loc='upper left')


def plotzOIIIWXum(table,col_WISE='WISE_nuLnu_rf8um',matchMagellan=True):
	clf()
	selection=table[col_WISE]<1.e43
	plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker='.',color='black',label='log '+col_WISE+' < 43')
	
	selection=all([table[col_WISE]>1.e43,table[col_WISE]<1.e44],axis=0)
	plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker='.',color='green',label='log '+col_WISE+' = 43-44')
	
	selection=all([table[col_WISE]>1.e44,table[col_WISE]<1.e46],axis=0)
	plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker='.',color='blue',label='log '+col_WISE+' = 44-45',ms=8)
	
	selection=table[col_WISE]>1.e45
	plot(table['Z'][selection],log10(table['OIII_5007T_LUM'][selection]),ls='',marker='.',color='red',label='log '+col_WISE+' > 45',ms=10)
	
	xlabel('z')
	ylabel('L_OIII (Mullaney) [erg/s]')
	legend(loc='upper left')

	if matchMagellan:
		filename='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/Magellan_2014June/obsprep/targets_observed.txt'
		targets_Magellan=Table.read(filename,format='ascii')

		for i in range(len(targets_Magellan)):
			match=table[absolute(table['RA']-targets_Magellan[i]['RA'])<0.001]
			match=match[absolute(match['DEC']-targets_Magellan[i]['DEC'])<0.001]
			id=targets_Magellan[i]['number']
			if len(match)==1: text(match['Z'],log10(match['OIII_5007T_LUM']),str(id),color='black')
			else: print id, 'not matched',len(match)




def overplotMagellan(table):
	#============ read targets Magellan
	filename='/Users/aisun/Documents/Astro/Thesis/bbselection/sample/Mullaney/Magellan_2014June/obsprep/targets_observed.txt'
	targets_Magellan=Table.read(filename,format='ascii')
	#========
	for i in range(len(targets_Magellan)):
		match=table[absolute(table['RA']-targets_Magellan[i]['RA'])<0.001]
		match=match[absolute(match['DEC']-targets_Magellan[i]['DEC'])<0.001]
		id=targets_Magellan[i]['number']
		if len(match)==1:
			if match['W4_LUM_MagAB']<-29: color='red'
			elif match['W4_LUM_MagAB']<-28 and match['W4_LUM_MagAB']>-29: color='blue'
			elif match['W4_LUM_MagAB']<-27 and match['W4_LUM_MagAB']>-28: color='green'
			else: color='black'
			text(match['Z'],log10(match['OIII_5007T_LUM']),str(id),color=color)
			if targets_Magellan[i]['Magellan_FUN'] == 0:
				text(match['Z'],log10(match['OIII_5007T_LUM']),'====',color=color)
		else:
			print id, 'not matched',len(match)


#========== RADEC select


def selectRADEC(table,list):
	imatches=array([])
	for i in range(len(list)):
		imatch=where(all([[absolute(table['RA']-list[i]['RA'])<0.001],
						  [absolute(table['DEC']-list[i]['DEC'])<0.001]],axis=0)[0])[0]
		if imatch.size ==1:
			imatches=append(imatches,[imatch])
			print 'MATCHED',list[i]['RA'],list[i]['DEC']
		if imatch.size >1:
			imatches=append(imatches,[imatch[0]])
			print 'WARNING: CONFUSION',list[i]['RA'],list[i]['DEC']
		if imatch.size ==0:
			print 'WARNING: NO MATCH',list[i]['RA'],list[i]['DEC']
	imatches=unique(imatches)
	return table[imatches.tolist()]


def selectionRADEC(table,list):
	if any(isinf(table['RA'])): raise NameError('Error: infinite RA')
	selection=isinf(table['RA'])
	for i in range(len(list)):
		newsel=all([[absolute(table['RA']-list[i]['RA'])<0.001],
						  [absolute(table['DEC']-list[i]['DEC'])<0.001]],axis=0)[0]
		selection=any([selection,newsel],axis=0)
	return selection



def RADECmatchtables(table,list):
	# extra columns in list will be added to table
	extra_col=array([])
	common_col=array([])
	for i in range(len(list.dtype.names)):
		if not list.dtype.names[i] in table.dtype.names:
			col=Column(name=list.dtype.names[i],dtype=list.dtype[i],length=len(table))
			table.add_column(col)
			extra_col=append(extra_col,[list.dtype.names[i]])
		elif list.dtype.names[i]!='RA' and list.dtype.names[i] != 'DEC':
			common_col=append(common_col,[list.dtype.names[i]])
	print "extra column", extra_col
	print "common column",common_col
	for i in range(len(list)):
		imatch=where(all([[absolute(table['RA']-list[i]['RA'])<0.001],
					   [absolute(table['DEC']-list[i]['DEC'])<0.001]],axis=0)[0])[0]
		if imatch.size ==1:
			for j in range(len(extra_col)):
				table[extra_col[j]][imatch]=list[i][extra_col[j]]
			for j in range(len(common_col)):
				table[common_col[j]][imatch]=list[i][common_col[j]]
			print "MATCHED",list[i]['RA'],list[i]['DEC']
		if imatch.size >1:
			for j in range(len(extra_col)):
				table[extra_col[j]][imatch[0]]=list[i][extra_col[j]]
			for j in range(len(common_col)):
				table[common_col[j]][imatch[0]]=list[i][common_col[j]]
			print 'WARNING: CONFUSION',list[i]['RA'],list[i]['DEC']
		if imatch.size ==0:
			print 'WARNING: NO MATCH',list[i]['RA'],list[i]['DEC']
	return table




#============== Luminosities


def LFIRest(Fnu60,Fnu100,Dl):
	# estimate LFIR using
	# Solomon 97 eq. 4, assuming r(S60/S100) = 1.8
	return 3.94e5*u.Lsun*(2.58*(Fnu60/u.Jy)+(Fnu100/u.Jy))*1.8*(Dl/u.Mpc)**2

def LCO10_FIR(LFIR):
	# estimated L'CO10 using Villar-Martin+13 Fig. 4 top panel regression line
	return 80.*(u.K*kms*u.pc**2)*(LFIR/u.Lsun)**0.65

def MagAB2Lnu(MagAB):
	return 10.**((MagAB-51.533)/-2.5)*u.erg/u.s/u.Hz

def Lnu25toLMIR(Lnu25):
	# from Krips+12 originally from Dale & Helou (2002)
	# rough esimtate for z ~0.1
	nu = const.c/(25.e-6*u.m)
	return 1.75*nu*Lnu25


def nuLnu25toLMIR(nuLnu25):
	# from Krips+12 originally from Dale & Helou (2002)
	# rough esimtate for z ~0.1
	return 1.75*nuLnu25


def LOIII2Lbol(LOIII):
	# from Liu+09 eq. 1
	# estimates with +/- 0.5 dex error
	return 10.**3.5*(LOIII/(3.846e+33*u.erg/u.s))**0.99*u.Lsun.to(u.erg/u.s)*u.erg/u.s

def Lnu25toLbol_MIR(Lnu25):
	# from Richards et al. (2006),
	# rough esimtate with bol correction of 11
	nu = const.c/(25.e-6*u.m)
	return 11.*nu*Lnu25


def nuLnu25toLbol_MIR(nuLnu25):
	# from Richards et al. (2006),
	# rough esimtate with bol correction of 11
	return 11.*nuLnu25


def magAB2Fnu(magAB):
	return 10.**(-(magAB+48.6)/2.5)*u.erg/u.cm**2/u.s/u.Hz


