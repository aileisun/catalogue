# catalogue_util.py
# ALS 04/07/2015
"""
Contains functions for various uses related to the catalogue

CONTENT:
	Name operations:

		addalRAhmsDECdms()
		

		addcolSDSSName()
		getSDSSName()
		getRADECfromSDSSName()
		gethhmmsss()
		getsddmmss()

	Write up files from catalogue
		writeAPOCatalogue()
		writeFIRSTsearch()
		writeJSkylist()
		writeMagellanCatalogue()
		writeSDSSSImgearch()
		writeSDSSfchartwget()
		writeSDSSstampwget()
		writeSummary()
		makefchart()
		makestamp()

	Table matching operation
		selectRADEC()
		selectionRADEC()
		RADECmatchtables()

	Plotting tools
		plotoverMagellan()
		plotoverMagellan()
		plotzOIIIW4()
		plotzOIIIWXum()


"""

# from pylab import *

import numpy as np
import astropy.table as at
import matplotlib.pyplot as plt

import astropy.units as u

from astropy.coordinates import SkyCoord


#========== name operatoin 

def addcolSDSSNAME(tab):
	"""
	add column SDSSNAME to table
	"""
	tab['SDSSNAME']='              '
	for i in range(len(tab)):
		tab['SDSSNAME'][i]=getSDSSName_fromlist(tab['RA','DEC'][i])

def getSDSSName_fromlist(listin):
	"""
	return SDSSName (string) from a list with 'RA' and 'DEC'
	"""
	from astropy.coordinates import SkyCoord
	c = SkyCoord(listin['RA'], listin['DEC'], 'icrs', unit='deg')
	return getSDSSName(c)

def getSDSSName_fromRADEC(ra, dec):
	"""
	return SDSS name given 
	PRAMS
	------
	ra (float) in deg
	dec (float) in deg
	"""
	from astropy.coordinates import SkyCoord
	c = SkyCoord(ra, dec, 'icrs', unit='deg')
	return getSDSSName(c)


def getSDSSName(c):
	# return "SDSSJ"+"%02d"%c.ra.hms[0]+"%02d"%c.ra.hms[1]+("%019.16f"%c.ra.hms[2])[0:5]+"%+03d"%c.dec.dms[0]+"%02d"%c.dec.dms[1]+("%019.16f"%c.dec.dms[2])[0:4]
	return "SDSSJ"+"%02d"%c.ra.hms[0]+"%02d"%c.ra.hms[1]+"%+03d"%c.dec.dms[0]+"%02d"%abs(c.dec.dms[1])

def	getSkyCoordfromSDSSName(SDSSName):
	"""
	get RA DEC from SDSS of format like 'SDSSJ005621.72+003235.7'
	output return is SkyCoord object c
	"""
	if len(SDSSName)==23:
		SDSSName=SDSSName[4:]
		RAstr=SDSSName[1:3]+":"+SDSSName[3:5]+":"+SDSSName[5:10]
		DECstr=SDSSName[10:13]+":"+SDSSName[13:15]+":"+SDSSName[15:]
		c = SkyCoord(RAstr, DECstr, frame='icrs', unit=(u.hourangle, u.deg))
		return c
	else:
		raise ValueError("invalid SDSSName input")


def addcolSDSSName(table):
	colname='SDSSNAME'
	if colname not in table.dtype.names:
		print "creating column ", colname
		table[colname]='                  '

	if 'ra' in table.colnames:
		ra = table['ra']
		dec = table['dec']
	elif 'RA' in table.colnames:
		ra = table['RA']
		dec = table['DEC']

	for i in range(len(table)):
		c = SkyCoord(ra=ra[i], dec=dec[i], unit=(u.degree, u.degree))
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



#=========== write files for preparation of observation
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
		name=np.arange(len(table))
	for i in range(len(table)):
		file.write(str(name[i])+' ')
		c = SkyCoord(ra=table[i]['RA'], dec=table[i]['DEC'], unit=(u.degree, u.degree))
		file.write('%02d'%c.ra.hms[0]+' '+'%02d'%c.ra.hms[1]+' '+'%02d'%c.ra.hms[2]
				   +'.'+str(c.ra.hms[2]*10%10)[0]+' ')
		file.write('%+03d'%c.dec.dms[0]+' '+'%02d'%abs(c.dec.dms[1])+' '+'%02d'%abs(c.dec.dms[2])+' ')
		file.write('2000 \n')
	file.close()


#=========== write SDSS Image Search
# http://skyserver.sdss3.org/public/en/tools/chart/list.aspx
def writeSDSSSImgearch(filename,table,header=''):
	if 'Name' in table.dtype.names:
		name=table['Name']
		savetxt(filename, np.array([name,table['RA'],table['DEC']]).transpose(),delimiter=' ', comments='#', fmt=['%s','%.3f','%.3f'],header=header)
	if 'name' in table.dtype.names:
		name=table['name']
		savetxt(filename, np.array([name,table['RA'],table['DEC']]).transpose(),delimiter=' ', comments='#', fmt=['%s','%.3f','%.3f'],header=header)
	else:
		name=np.arange(len(table)).astype(int)
		savetxt(filename, np.array([name,table['RA'],table['DEC']]).transpose(),delimiter=' ', comments='#', fmt=['%i','%.3f','%.3f'],header=header)


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
		name=np.arange(len(table)).astype('str')
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
		name=np.arange(len(table)).astype('str')
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
	if np.any(table.dtype.names=='Name'): name=table['Name']
	elif np.any(table.dtype.names=='name'): name=table['name']
	else: name=np.arange(table.size).astype('str')

	if np.any(table.dtype.names=='PA'): PA=table['PA']
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


def writeSummary(filename,table):
	if 'N' not in table.dtype.names:
		table['N']=np.arange(len(table))
	# colnames=['N','OBJID','RA','DEC','Z','AGN_TYPE','OIII_5007T_LUM','OIII_5007AVG_FWHM','W4_LUM_MagAB']
	colnames=['OBJID','SDSSNAME','RA','DEC','Z','AGN_TYPE', 'OIII_5007T_LUM','OIII_5007AVG_FWHM','Lbol_OIII','W4_nuLnu_rf','Lbol_22um','LFIR','FIRST_LUM']
	table[colnames].write(filename,format='ascii.fixed_width',delimiter='')



#========== select and matching

def selectRADEC(table, list, radius=5., verbos=True, include_all_confusions=False):
	"""
	PURPOSE: return the rows in the table that has RA DEC matching the list within search radius
	INPUT:
			radius=5. : search radius in units of arcsecs
	"""
	# find radius in units of degree
	radius_deg=radius/3600.

	imatches=np.array([])
	for i in range(len(list)):
		# cone search
		imatch = np.where(np.sqrt((table['RA']-list[i]['RA'])**2+(table['DEC']-list[i]['DEC'])**2)<radius_deg)[0]
		# square search -- abandoned
		# imatch=np.where(np.all([[np.absolute(table['RA']-list[i]['RA'])<0.001],[np.absolute(table['DEC']-list[i]['DEC'])<0.001]],axis=0)[0])[0]
		if imatch.size == 1:
			imatches = np.append(imatches,[imatch])
			if verbos: print 'MATCHED',list[i]['RA'],list[i]['DEC']
		if imatch.size >1 :
			if not include_all_confusions:
				imatches = np.append(imatches, imatch[0])
				if verbos: print 'CONFUSION only one included',list[i]['RA'],list[i]['DEC']
			else: 
				imatches = np.append(imatches, imatch)
				if verbos: print 'CONFUSION all included',list[i]['RA'],list[i]['DEC']
		if imatch.size == 0:
			if verbos: print 'NO MATCH',list[i]['RA'],list[i]['DEC']
	imatches = np.unique(imatches).astype('int')
	return table[imatches]


def selectionRADEC(table,list):
	if np.any(np.isinf(table['RA'])): raise NameError('Error: infinite RA')
	selection=np.isinf(table['RA'])
	for i in range(len(list)):
		newsel=np.all([[np.absolute(table['RA']-list[i]['RA'])<0.001],
						  [np.absolute(table['DEC']-list[i]['DEC'])<0.001]],axis=0)[0]
		selection=np.any([selection,newsel],axis=0)
	return selection



def RADECmatchtables(table,list):
	# extra columns in list will be added to table
	extra_col=np.array([])
	common_col=np.array([])
	for i in range(len(list.dtype.names)):
		if not list.dtype.names[i] in table.dtype.names:
			col=at.Column(name=list.dtype.names[i],dtype=list.dtype[i],length=len(table))
			table.add_column(col)
			extra_col = np.append(extra_col,[list.dtype.names[i]])
		elif list.dtype.names[i]!='RA' and list.dtype.names[i] != 'DEC':
			common_col = np.append(common_col,[list.dtype.names[i]])
	print "extra column", extra_col
	print "common column",common_col
	for i in range(len(list)):
		imatch=np.where(np.all([[np.absolute(table['RA']-list[i]['RA'])<0.001],
					   [np.absolute(table['DEC']-list[i]['DEC'])<0.001]],axis=0)[0])[0]
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




#=========== plotting

def plotzOIIIW4(table,marker='.',toclf=True,labelsuffix='',plotFWHM500=True):
	if toclf: clf()
	markersize=ceil((-table['W4_LUM_MagAB']-23)).astype(int)
	if plotFWHM500 and 'OIII_5007AVG_FWHM' in table.dtype.names:
		selFWHM500=table['OIII_5007AVG_FWHM']>500
		plt.plot(table['Z'][selFWHM500],np.log10(table['OIII_5007T_LUM'][selFWHM500]),ls='',marker='+',ms=6,color='black',label='OIII FWHM >500'+labelsuffix)
	for i in range(max(markersize)+1):
		selection= markersize==i
		plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker=marker,ms=i,color='black')
		if i ==7:
			plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker=marker,ms=i,color='red',label='W4 = -29~-30'+labelsuffix)
			if plotFWHM500 and 'OIII_5007AVG_FWHM' in table.dtype.names:
				selection=np.all([selection,selFWHM500],axis=0)
				plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker='+',ms=i,color='red')
		if i ==6:
			plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker=marker,ms=i,color='blue',label='W4 = -28~-29'+labelsuffix)
		if i ==5:
			plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker=marker,ms=i,color='green',label='W4 = -27~-28'+labelsuffix)
	plt.xlabel('z')
	plt.ylabel('log10 LOIII')
	plt.legend(loc='upper left')


def plotzOIIIWXum(table,col_WISE='WISE_nuLnu_rf8um',matchMagellan=True):
	plt.clf()
	selection=table[col_WISE]<1.e43
	plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker='.',color='black',label='log '+col_WISE+' < 43')
	
	selection=np.all([table[col_WISE]>1.e43,table[col_WISE]<1.e44],axis=0)
	plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker='.',color='green',label='log '+col_WISE+' = 43-44')
	
	selection=np.all([table[col_WISE]>1.e44,table[col_WISE]<1.e46],axis=0)
	plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker='.',color='blue',label='log '+col_WISE+' = 44-45',ms=8)
	
	selection=table[col_WISE]>1.e45
	plt.plot(table['Z'][selection],np.log10(table['OIII_5007T_LUM'][selection]),ls='',marker='.',color='red',label='log '+col_WISE+' > 45',ms=10)
	
	plt.xlabel('z')
	plt.ylabel('L_OIII (Mullaney) [erg/s]')
	plt.legend(loc='upper left')

	if matchMagellan:
		filename='/Users/aisun/Documents/Astro/Thesis/bbselection/SDSS/sample/Mullaney/Magellan_2014June/obsprep/targets_observed.txt'
		targets_Magellan=Table.read(filename,format='ascii')

		for i in range(len(targets_Magellan)):
			match=table[np.absolute(table['RA']-targets_Magellan[i]['RA'])<0.001]
			match=match[np.absolute(match['DEC']-targets_Magellan[i]['DEC'])<0.001]
			id=targets_Magellan[i]['number']
			if len(match)==1: text(match['Z'],np.log10(match['OIII_5007T_LUM']),str(id),color='black')
			else: print id, 'not matched',len(match)




def plotoverMagellan(table):
	#============ read targets Magellan
	filename='/Users/aisun/Documents/Astro/Thesis/bbselection/sample/Mullaney/Magellan_2014June/obsprep/targets_observed.txt'
	targets_Magellan=Table.read(filename,format='ascii')
	#========
	for i in range(len(targets_Magellan)):
		match=table[np.absolute(table['RA']-targets_Magellan[i]['RA'])<0.001]
		match=match[np.absolute(match['DEC']-targets_Magellan[i]['DEC'])<0.001]
		id=targets_Magellan[i]['number']
		if len(match)==1:
			if match['W4_LUM_MagAB']<-29: color='red'
			elif match['W4_LUM_MagAB']<-28 and match['W4_LUM_MagAB']>-29: color='blue'
			elif match['W4_LUM_MagAB']<-27 and match['W4_LUM_MagAB']>-28: color='green'
			else: color='black'
			text(match['Z'],np.log10(match['OIII_5007T_LUM']),str(id),color=color)
			if targets_Magellan[i]['Magellan_FUN'] == 0:
				text(match['Z'],np.log10(match['OIII_5007T_LUM']),'====',color=color)
		else:
			print id, 'not matched',len(match)
