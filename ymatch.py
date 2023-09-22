import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.coordinates import match_coordinates_sky
import pandas as pd
from astropy.table import Table, join, join_skycoord,vstack
from astropy.io import fits
import matplotlib.pyplot as plt 
import sys, glob, os
import time
import logging
from loguru import logger
import argparse
from locale import ALT_DIGITS
from lib.LogInstance import Logger
from execute.yphotutils import datename 
from lib.phot.ybias import get_current_dir
import traceback
from astropy.stats import sigma_clip,sigma_clipped_stats
import ntpath
from execute.yphotutils import datename,d2hms,mad,wds9reg,read_param,crossmatch,reg,ref_expt_select,readfits,hdr_edit_fits,re_hdr_edit_fits,read_list,timeplus,HJDTrans,get_scamp_head,overlapRect,pointRect,_spherical_to_cartesian,_great_circle_distance,mkdir

def match_gaiadr3_50cm(source,path_cat,match_r = 3.0):
	if not(os.path.exists(path_cat)):
		logger.info('there is no such directory')
		return None
	name_catalog='gaiadr3'
	#Read the center position from header
	fitsname = source.split('_sexcat')[0]+'.fits'
	# logger.info('#########################')
	# logger.info(source)
	#logger.info(fitsname)
	header=fits.getheader(fitsname)
	try:
		ra_cen=header['OBJCTRA']
		dec_cen=header['OBJCTDEC']
	except:
		ra_cen=header['RA']
		dec_cen=header['DEC']
	
	c0 = SkyCoord(ra_cen, dec_cen, unit=(u.hourangle, u.deg),frame='icrs')  
	#Read the position (ra, dec) of stars in the file 
	
	data0=Table.read(source,hdu=2)
	flag=data0['FLAGS']
	idx=np.where(flag==0)[0]
	data=data0[idx]
	ra1=data['ALPHA_J2000']
	dec1=data['DELTA_J2000']
	c = SkyCoord(ra=np.array(ra1)*u.degree, dec=np.array(dec1)*u.degree,frame='icrs')
	#Find the catalog files that contains the center position of the target. 
	
	catinfo=fits.open(path_cat+'gaiadr3catalog.info') #Gaia DR3 catalogue
	cramax=catinfo[1].data['ramax']
	cramin=catinfo[1].data['ramin']
	cdecmax=catinfo[1].data['decmax']
	cdecmin=catinfo[1].data['decmin']
	idxca=np.where((c0.ra.degree>=cramin-1/np.cos(np.radians(c0.dec.degree))) & (c0.ra.degree<=cramax+1/np.cos(np.radians(c0.dec.degree))) & (c0.dec.degree>=cdecmin-1) & (c0.dec.degree<=cdecmax+1))[0]

	catalogname=catinfo[1].data['name'][idxca]

	cattable=[Table()]*len(idxca)
	for i in range(len(idxca)):
		#logger.info(i,len(idxca))
		#cattablei=pd.read_table(path_cat+catalogname[i],header=1000,sep=',')
		#cattablei=Table.from_pandas(cattablei) 
		cattablei=Table.read(path_cat+catalogname[i])
		rai=np.array(cattablei['ra_J2023'])
		deci=np.array(cattablei['dec_J2023'])  
		 
		idxra=np.where(rai!=999.0)[0]
		rai=rai[idxra]
		deci=deci[idxra]  
		cattablei=cattablei[idxra]
		clogi=SkyCoord(ra=rai*u.degree, dec=deci*u.degree)
		region=c0.separation(clogi)<40*u.arcminute
		idxdi=np.where(region==True)[0]
		if len(idxdi)==0:
			continue
		else:
			tablei=cattablei[idxdi]
			cattable[i]=tablei
	catalog=cattable[0]
	for j in range(1,len(cattable)):
		catalog=vstack([catalog, cattable[j]])
	 
	# Matching with catalog 
	if len(catalog)==0:
		logger.info('No match catalog data')
		return None
	else:
		ra2=np.array(catalog['ra_J2023'])
		dec2=np.array(catalog['dec_J2023'])      
		clog=SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
		start2=time.time()
		 
		 
		 
		idx, d2d, d3d = c.match_to_catalog_sky(clog) 

		index=np.where(d2d.arcsec<match_r)[0]

		tablel=data[index]
		tabler=catalog[idx[index]]
		if len(tablel)==0 or len(tabler)==0 :
		
			logger.info('No match data')
			return None
		number=tablel['NUMBER']
		offset=d2d.arcsec[index]
		 
		tabler.add_column(offset,name='sepdis')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,join_type='left',keys='NUMBER')  
		ra_mp=newtable['ALPHA_J2000'].data
		dec_mp=newtable['DELTA_J2000'].data
		ra_gaia=newtable['ra_J2023'].data
		dec_gaia=newtable['dec_J2023'].data
		deltara=(ra_mp*np.cos(np.radians(dec_mp))-ra_gaia*np.cos(np.radians(dec_gaia)))*3600.0
		deltadec=(dec_mp-dec_gaia)*3600.0

		return newtable

def match_ps2_50cm(source,path_cat,match_r = 3.0):
	name_catalog='ps2'
	
	#Read the center position from header
	fitsname = source.split('_sexcat')[0]+'.fits'
	 
	 
	#logger.info(fitsname)
	header=fits.getheader(fitsname)
	try:
		ra_cen=header['OBJCTRA']
		dec_cen=header['OBJCTDEC']
	except:
		ra_cen=header['RA']
		dec_cen=header['DEC']
	
	c0 = SkyCoord(ra_cen, dec_cen, unit=(u.hourangle, u.deg),frame='icrs')  
	#Read the position (ra, dec) of stars in the file 
	
	data0=Table.read(source,hdu=2)
	flag=data0['FLAGS']
	idx=np.where(flag==0)[0]
	data=data0[idx]
	ra1=data['ALPHA_J2000']
	dec1=data['DELTA_J2000']
	c = SkyCoord(ra=np.array(ra1)*u.degree, dec=np.array(dec1)*u.degree,frame='icrs')
	#Find the file containing 40 arcmin near the center position from the catalog files
	
	catinfo=fits.open(path_cat+'ps2catalog.info') #ps2 catalogue
	cramax=catinfo[1].data['ramax']
	cramin=catinfo[1].data['ramin']
	cdecmax=catinfo[1].data['decmax']
	cdecmin=catinfo[1].data['decmin']
	idxca=np.where((c0.ra.degree>=cramin-1/np.cos(np.radians(c0.dec.degree))) & (c0.ra.degree<=cramax+1/np.cos(np.radians(c0.dec.degree))) & (c0.dec.degree>=cdecmin-1) & (c0.dec.degree<=cdecmax+1))[0]

	catalogname=catinfo[1].data['name'][idxca]

	cattable=[Table()]*len(idxca)
	for i in range(len(idxca)):
		cattablei=Table.read(path_cat+catalogname[i])
		rai=np.array(cattablei['raMean_J2023']).reshape(1,-1)[0]
		deci=np.array(cattablei['decMean_J2023']).reshape(1,-1)[0]  
		# rai=np.array(cattablei['ra_J2023'])
		# deci=np.array(cattablei['dec_J2023'])      
		idxra=np.where(rai!=999.0)[0]
		rai=rai[idxra]
		deci=deci[idxra]  
		cattablei=cattablei[idxra]
		clogi=SkyCoord(ra=rai*u.degree, dec=deci*u.degree)
		region=c0.separation(clogi)<40*u.arcminute
		idxdi=np.where(region==True)[0]
		if len(idxdi)==0:
			continue
		else:
			tablei=cattablei[idxdi]
			cattable[i]=tablei
	catalog=cattable[0]
	for j in range(1,len(cattable)):
		catalog=vstack([catalog, cattable[j]])
	time2=time.time()
	
	# Matching with catalog 
	if len(catalog)==0:
		logger.info('No match catalog data')
		return None
	else:
		ra2=np.array(catalog['raMean_J2023']).reshape(1,-1)[0]
		dec2=np.array(catalog['decMean_J2023']).reshape(1,-1)[0] 
		     
		clog=SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
		start2=time.time()
		 
		idx, d2d, d3d = c.match_to_catalog_sky(clog) 

		index=np.where(d2d.arcsec<match_r)[0]

		tablel=data[index]
		tabler=catalog[idx[index]]
		if len(tablel)==0 or len(tabler)==0 :
		
			logger.info('No match data')
			return None
		time3=time.time()
		number=tablel['NUMBER']
		offset=d2d.arcsec[index]
		tabler.add_column(offset,name='sepdis')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,join_type='left',keys='NUMBER')  
		ra_mp=newtable['ALPHA_J2000'].data
		dec_mp=newtable['DELTA_J2000'].data
		ra_ps2=np.array(newtable['raMean_J2023']).reshape(1,-1)[0]
		dec_ps2=np.array(newtable['decMean_J2023']).reshape(1,-1)[0]
		deltara=(ra_mp*np.cos(np.radians(dec_mp))-ra_ps2*np.cos(np.radians(dec_ps2)))*3600.0
		deltadec=(dec_mp-dec_ps2)*3600.0
		#newtable.add_column(deltara,name='deltara')
		#newtable.add_column(deltadec,name='deltadec')
		time4=time.time()
		logger.info(f'Original Num:{len(data)} , Xmatch with ps2 Num:{len(tablel)}')	
			
		return newtable

def match_2mass_50cm(source,path_cat,match_r = 3.0):
	name_catalog='2mass'
	#path_catalog='/media/huohe/hlserver/'
	#catalog=glob.glob(path_catalog+'gaiasource/'+'*.gz')
	#path_file='/home/50cm/50cm_scripts/50cm_pro/reception/'+date+'/sci/'+source+'/'
	#filename=glob.glob(path_file+'*subbkg_sexcat.fits')
	#Read the center position from header
	#header=fits.getheader(source[:-21]+'.fits')
	fitsname = source.split('_sexcat')[0]+'.fits'
 
	header=fits.getheader(fitsname)
    	
	try:
		ra_cen=header['OBJCTRA']
		dec_cen=header['OBJCTDEC']
	except:
		ra_cen=header['RA']
		dec_cen=header['DEC']
	c0 = SkyCoord(ra_cen, dec_cen, unit=(u.hourangle, u.deg),frame='icrs')  
	#Read the position (ra, dec) of stars in the file 
	
	data0=Table.read(source,hdu=2)
	flag=data0['FLAGS']
	idx=np.where(flag==0)[0]
	data=data0[idx]
	ra1=data['ALPHA_J2000']
	dec1=data['DELTA_J2000']
	c = SkyCoord(ra=np.array(ra1)*u.degree, dec=np.array(dec1)*u.degree,frame='icrs')
	#Find the file containing 40 arcmin near the center position from the catalog files
	
	catinfo=fits.open(path_cat+'2masscatalog.info') #2mass catalogue
	cramax=catinfo[1].data['ramax']
	cramin=catinfo[1].data['ramin']
	cdecmax=catinfo[1].data['decmax']
	cdecmin=catinfo[1].data['decmin']
	idxca=np.where((c0.ra.degree>=cramin-1/np.cos(np.radians(c0.dec.degree))) & (c0.ra.degree<=cramax+1/np.cos(np.radians(c0.dec.degree))) & (c0.dec.degree>=cdecmin-1) & (c0.dec.degree<=cdecmax+1))[0]

	catalogname=catinfo[1].data['name'][idxca]

	cattable=[Table()]*len(idxca)
	for i in range(len(idxca)):
		#logger.info(i,len(idxca))
		#cattablei=pd.read_table(path_cat+catalogname[i],header=1000,sep=',')
		#cattablei=Table.from_pandas(cattablei) 
		cattablei=Table.read(path_cat+catalogname[i])
		rai=np.array(cattablei['TMASS_RA_J2023'])
		deci=np.array(cattablei['TMASS_DEC_J2023'])      
		clogi=SkyCoord(ra=rai*u.degree, dec=deci*u.degree)
		region=c0.separation(clogi)<40*u.arcminute
		idxdi=np.where(region==True)[0]
		if len(idxdi)==0:
			continue
		else:
			tablei=cattablei[idxdi]
			cattable[i]=tablei
	catalog=cattable[0]
	for j in range(1,len(cattable)):
		catalog=vstack([catalog, cattable[j]])
	# Matching with catalog 
	if len(catalog)==0:
		logger.info('No match catalog data')
		return None
	else:
		ra2=np.array(catalog['TMASS_RA_J2023'])
		dec2=np.array(catalog['TMASS_DEC_J2023'])      
		clog=SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
		start2=time.time()
		max_sep = 5.0 
	 
		 
		idx, d2d, d3d = c.match_to_catalog_sky(clog) 

		index=np.where(d2d.arcsec<match_r)[0]

		tablel=data[index]
		tabler=catalog[idx[index]]
		if len(tablel)==0 or len(tabler)==0 :
		
			logger.info('No match data')
			return None
		number=tablel['NUMBER']
		offset=d2d.arcsec[index]
		tabler.add_column(offset,name='sepdis')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,join_type='left',keys='NUMBER')  
		ra_mp=newtable['ALPHA_J2000'].data
		dec_mp=newtable['DELTA_J2000'].data
		ra_2mass=newtable['TMASS_RA_J2023'].data
		dec_2mass=newtable['TMASS_DEC_J2023'].data
		deltara=(ra_mp*np.cos(np.radians(dec_mp))-ra_2mass*np.cos(np.radians(dec_2mass)))*3600.0
		deltadec=(dec_mp-dec_2mass)*3600.0
		#newtable.add_column(deltara,name='deltara')
		#newtable.add_column(deltadec,name='deltadec')
		logger.info(f'Original Num:{len(data)} , Xmatch with 2MASS Num:{len(tablel)}')
		return newtable			

def match_wise_50cm(source,path_cat,match_r = 3.0):
	name_catalog='wise'
	
	#Read the center position from header
	fitsname = source.split('_sexcat')[0]+'.fits'
 
	header=fits.getheader(fitsname)
	try:
		ra_cen=header['OBJCTRA']
		dec_cen=header['OBJCTDEC']
	except:
		ra_cen=header['RA']
		dec_cen=header['DEC']
	c0 = SkyCoord(ra_cen, dec_cen, unit=(u.hourangle, u.deg),frame='icrs')  
	#Read the position (ra, dec) of stars in the file 
	
	data0=Table.read(source,hdu=2)
	flag=data0['FLAGS']
	idx=np.where(flag==0)[0]
	data=data0[idx]
	ra1=data['ALPHA_J2000']
	dec1=data['DELTA_J2000']
	c = SkyCoord(ra=np.array(ra1)*u.degree, dec=np.array(dec1)*u.degree,frame='icrs')
	#Find the file containing 40 arcmin near the center position from the catalog files
	
	catinfo=fits.open(path_cat+'wisecatalog.info') #allwise catalogue
	cramax=catinfo[1].data['ramax']
	cramin=catinfo[1].data['ramin']
	cdecmax=catinfo[1].data['decmax']
	cdecmin=catinfo[1].data['decmin']
	idxca=np.where((c0.ra.degree>=cramin-1/np.cos(np.radians(c0.dec.degree))) & (c0.ra.degree<=cramax+1/np.cos(np.radians(c0.dec.degree))) & (c0.dec.degree>=cdecmin-1) & (c0.dec.degree<=cdecmax+1))[0]

	catalogname=catinfo[1].data['name'][idxca]

	cattable=[Table()]*len(idxca)
	for i in range(len(idxca)):
		
		#cattablei=pd.read_table(path_cat+catalogname[i],header=1000,sep=',')
		#cattablei=Table.from_pandas(cattablei) 
		cattablei=Table.read(path_cat+catalogname[i])
		rai=np.array(cattablei['ra'])
		deci=np.array(cattablei['dec'])      
		clogi=SkyCoord(ra=rai*u.degree, dec=deci*u.degree)
		region=c0.separation(clogi)<40*u.arcminute
		idxdi=np.where(region==True)[0]
		if len(idxdi)==0:
			continue
		else:
			tablei=cattablei[idxdi]
			cattable[i]=tablei
	catalog=cattable[0]
	for j in range(1,len(cattable)):
		catalog=vstack([catalog, cattable[j]])
	# Matching with catalog 
	if len(catalog)==0:
		logger.info('No match catalog data')
		return None
	else:
		ra2=np.array(catalog['ra'])
		dec2=np.array(catalog['dec'])      
		clog=SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
		start2=time.time()
		
	 
		 
		idx, d2d, d3d = c.match_to_catalog_sky(clog) 

		index=np.where(d2d.arcsec<match_r)[0]

		tablel=data[index]
		tabler=catalog[idx[index]]
		if len(tablel)==0 or len(tabler)==0 :
		
			logger.info('No match data')
			return None
		number=tablel['NUMBER']
		offset=d2d.arcsec[index]
		tabler.add_column(offset,name='sepdis')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,join_type='left',keys='NUMBER')  
		ra_mp=newtable['ALPHA_J2000'].data
		dec_mp=newtable['DELTA_J2000'].data
		ra_wise=newtable['ra'].data
		dec_wise=newtable['dec'].data
		deltara=(ra_mp*np.cos(np.radians(dec_mp))-ra_wise*np.cos(np.radians(dec_wise)))*3600.0
		deltadec=(dec_mp-dec_wise)*3600.0
		#newtable.add_column(deltara,name='deltara')
		#newtable.add_column(deltadec,name='deltadec')
		logger.info(f'Original Num:{len(data)} , Xmatch with WISE Num:{len(tablel)}')
		return newtable	
 

def match_sdss_50cm(source,path_cat,match_r = 2.0):
	 
	name_catalog='sdssdr17'
	
	#Read the center position from header
	fitsname = source.split('_sexcat')[0]+'.fits'
	header=fits.getheader(fitsname)
	try:
		ra_cen=header['OBJCTRA']
		dec_cen=header['OBJCTDEC']
	except:
		ra_cen=header['RA']
		dec_cen=header['DEC']
	try:
		ra_cen=float(ra_cen)
		dec_cen=float(dec_cen)
		c0 = SkyCoord(ra_cen*u.deg, dec_cen*u.deg,frame='icrs') 
	except:	
		c0 = SkyCoord(ra_cen, dec_cen, unit=(u.hourangle, u.deg),frame='icrs')  
	#Read the position (ra, dec) of stars in the file 
	
	data0=Table.read(source,hdu=2)
	flag=data0['FLAGS']
	idx=np.where(flag==0)[0]
	data=data0[idx]
	ra1=data['ALPHA_J2000']
	dec1=data['DELTA_J2000']
	c = SkyCoord(ra=np.array(ra1)*u.degree, dec=np.array(dec1)*u.degree,frame='icrs')
	#Find the file containing 40 arcmin near the center position from the catalog files
	
	catinfo=fits.open(path_cat+'sdssdr17catalog.info') #allwise catalogue
	cramax=catinfo[1].data['ramax']
	cramin=catinfo[1].data['ramin']
	cdecmax=catinfo[1].data['decmax']
	cdecmin=catinfo[1].data['decmin']
	idxca=np.where((c0.ra.degree>=cramin-1/np.cos(np.radians(c0.dec.degree))) & (c0.ra.degree<=cramax+1/np.cos(np.radians(c0.dec.degree))) & (c0.dec.degree>=cdecmin-1) & (c0.dec.degree<=cdecmax+1))[0]

	catalogname=catinfo[1].data['name'][idxca]

	cattable=[Table()]*len(idxca)
	for i in range(len(idxca)):
		
		#cattablei=pd.read_table(path_cat+catalogname[i],header=1000,sep=',')
		#cattablei=Table.from_pandas(cattablei) 
		cattablei=Table.read(path_cat+catalogname[i])
		rai=np.array(cattablei['ra'])
		deci=np.array(cattablei['dec'])      
		clogi=SkyCoord(ra=rai*u.degree, dec=deci*u.degree)
		region=c0.separation(clogi)<40*u.arcminute
		idxdi=np.where(region==True)[0]
		if len(idxdi)==0:
			continue
		else:
			tablei=cattablei[idxdi]
			cattable[i]=tablei
	catalog=cattable[0]
	for j in range(1,len(cattable)):
		catalog=vstack([catalog, cattable[j]])
	# Matching with catalog 
	if len(catalog)==0:
		print('No match catalog data')
		return None
	else:
		ra2=np.array(catalog['ra'])
		dec2=np.array(catalog['dec'])      
		clog=SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
		start2=time.time()
		
		#region=c0.separation(clog)<40*u.arcminute
		#idxd=np.where(region==True)[0]
		 
		idx, d2d, d3d = c.match_to_catalog_sky(clog) 

		index=np.where(d2d.arcsec<match_r)[0]

		tablel=data[index]
		tabler=catalog[idx[index]]
		if len(tablel)==0 or len(tabler)==0 :
		
			print('No match data')
			return None
		number=tablel['NUMBER']
		offset=d2d.arcsec[index]
		tabler.add_column(offset,name='sepdis')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,join_type='left',keys='NUMBER')  
		logger.info(f'Original Num:{len(data)} , Xmatch with SDSS Num:{len(tablel)}')
		#print('Original Num:',len(data),', Xmatch with sdss Num:',len(tablel))	
		return newtable







def xmach_fig(filename,figsavedir,stype='2mass'):
	hl=Table.read(filename)
	ra=hl['ALPHA_J2000'].data
	dec=hl['DELTA_J2000'].data
	# try:
	if(stype=='2mass'):
		ra_2mass=hl['TMASS_RA_J2023'].data
		dec_2mass=hl['TMASS_DEC_J2023'].data
	elif(stype=='ps2'):
		ra_2mass=np.array(hl['raMean_J2023']).reshape(1,-1)[0]
		dec_2mass=np.array(hl['decMean_J2023']).reshape(1,-1)[0]
	elif(stype=='gaia'):
		ra_2mass=np.array(hl['ra_J2023']).reshape(1,-1)[0]
		dec_2mass=np.array(hl['dec_J2023']).reshape(1,-1)[0]
	
	# elif(stype=='gaia'):
	# 	ra_2mass=hl['ALPHA_J2000'].data
	# 	dec_2mass=hl['DELTA_J2000'].data
	else:
		ra_2mass=hl['ALPHA_J2000'].data
		dec_2mass=hl['DELTA_J2000'].data
    
	deltara=(ra-ra_2mass)*np.cos(np.radians(dec))*3600.0
	deltadec=(dec-dec_2mass)*3600.0
	
	fig=plt.figure(figsize=(14,14))
	plt.subplot(2,2,1)
	#plt.plot(ra,dec,'.r')
	plt.scatter(ra_2mass,dec_2mass,s=4,color='r')
	for j in range(len(ra)):
		plt.arrow(ra_2mass[j],dec_2mass[j],deltara[j]/30,deltadec[j]/30,width = 0.001,ec='b',head_width=0.01)
	plt.scatter(np.median(ra),np.max(dec)+0.05,s=4,color='r')
	plt.arrow(np.median(ra),np.max(dec)+0.05,1/30,0.,width=0.001,ec='b',head_width=0.01)
	plt.text(np.median(ra)+9/180,np.max(dec)+0.06,'1 arcsec')
	plt.xlabel('RA (deg)',size=12)
	plt.ylabel('DEC (deg)',size=12)
	plt.gca().invert_xaxis()
	plt.subplot(2,2,2)
	try:
		plt.hist(hl['sepdis'],40)
	except:
		plt.hist(hl['offset'],40)
	plt.xlabel('Distance (arcsec)',size=12)
	plt.ylabel('Number',size=12)
	plt.subplot(2,2,3)
	clipra=sigma_clip(deltara,sigma=3,maxiters=1)
	deltara=deltara[clipra.mask==False]
	medra=round(np.median(deltara),2)
	stdra=round(np.std(deltara),2)
	histra=plt.hist(deltara,40)
	plt.plot([medra,medra],[0,np.max(histra[0])+5],'k')
	plt.plot([medra-stdra,medra-stdra],[0,np.max(histra[0])+5],'--k')
	plt.plot([medra+stdra,medra+stdra],[0,np.max(histra[0])+5],'--k')
	
	plt.xlabel('$\Delta$RA (arcsec)',size=12)
	plt.ylabel('Number',size=12)
	plt.title('$\mu$='+str(medra)+', $\sigma$='+str(stdra),color='r')
	plt.subplot(2,2,4)
	clipdec=sigma_clip(deltadec,sigma=3,maxiters=1)
	deltadec=deltadec[clipdec.mask==False]
	meddec=np.round(np.median(deltadec),2)
	stddec=np.round(np.std(deltadec),2)
	histdec=plt.hist(deltadec,40)
	plt.plot([meddec,meddec],[0,np.max(histdec[0])+5],'k')
	plt.plot([meddec-stddec,meddec-stddec],[0,np.max(histdec[0])+5],'--k')
	plt.plot([meddec+stddec,meddec+stddec],[0,np.max(histdec[0])+5],'--k')
	plt.xlabel('$\Delta$DEC (arcsec)',size=12)
	plt.ylabel('Number',size=12)
	plt.title('$\mu$='+str(meddec)+', $\sigma$='+str(stddec),color='r')
	plt.suptitle(ntpath.basename(filename))
	
	#mini_figsavedir = 
 
	 
	mini_figsavedir = figsavedir+'mini/'
	mkdir(mini_figsavedir)
	logger.info(os.path.join(figsavedir,ntpath.basename(filename)+'.png'))
	logger.info(os.path.join(mini_figsavedir,ntpath.basename(filename)+'.png'))
	fig.savefig(os.path.join(figsavedir,ntpath.basename(filename)+'.png'))
	fig.savefig(os.path.join(mini_figsavedir,ntpath.basename(filename)+'.png'))
	plt.close()
	return len(ra),medra,stdra,meddec,stddec 





def match_pro(date,tcsipath,refpath,figpath ):

    #filesavedir  = '/home/guohl/xmatch/xmatch50cm_%s/'%obsDate
    #if not os.path.exists(filesavedir): os.mkdir(filesavedir)
    #filedir='/data2/workspace/50cm_pro/reception/'+obsDate+'/sci/'
    #objects=glob.glob(filedir+'*')
    source=glob.glob(tcsipath+'*sexcat.fits')
    logger.info(source)
    path_catalog_gaia = refpath + 'GaiaDR3FITS_pm/'
    path_catalog_pans = refpath + 'Panstarrs2New/'
    path_catalog_2mass = refpath + '2massNew/'
    path_catalog_wise = refpath + 'allwiseFITS/'
    path_catalog_sdss = refpath + 'sdssFITS/'
      
    
    if len(source)>0:  
        for j in range(len(source)):
            #try:
            if(1==1):
                filaname = source[j]
                logger.info('%%%%%%%%%%%%%%%%%%%%%%%')
                logger.info(filaname)
                
                if(os.path.exists(path_catalog_gaia)):
                    try:
                        matchtable_gaia=match_gaiadr3_50cm(source=source[j],path_cat=path_catalog_gaia)
                        matchtable_gaia.write(filaname[:-5] +'_gaia.fits',format='fits',overwrite=True)
                        c_match,medra,stdra,meddec,stddec =xmach_fig(filaname[:-5] +'_gaia.fits',figpath,stype='gaia')
                        hdul = fits.open(filaname)
                        iheader1 = hdul[1].header
                        iheader1['ASTPRA'] = stdra
                        iheader1['ASTPDEC'] = stddec
                        iheader1['ASTMRA'] = medra
                        iheader1['ASTMDEC'] = meddec
                        iheader1['ASTCNT'] = c_match
                        fits.writeto(filaname[:-5] +'_gaia.fits',data=matchtable_gaia,header=iheader1,overwrite=True)
                        logger.info(f'{filaname[:-5]}_gaia.fits')
                    except Exception as e:
                        logger.error(traceback.format_exc())
                        #continue

                if(os.path.exists(path_catalog_2mass)):
                    try:
                        matchtable_2mass=match_2mass_50cm(source=source[j],path_cat=path_catalog_2mass)
                        matchtable_2mass.write(filaname[:-5] +'_2mass.fits',format='fits' ,overwrite=True)
                        c_match,medra,stdra,meddec,stddec =xmach_fig(filaname[:-5] +'_2mass.fits',figpath,stype='2mass')
                        hdul = fits.open(filaname)
                        iheader1 = hdul[1].header
                        iheader1['ASTPRA'] = stdra
                        iheader1['ASTPDEC'] = stddec
                        iheader1['ASTMRA'] = medra
                        iheader1['ASTMDEC'] = meddec
                        iheader1['ASTCNT'] = c_match
                        fits.writeto(filaname[:-5] +'_2mass.fits',data=matchtable_2mass,header=iheader1,overwrite=True)
                         
                    except Exception as e:
                        logger.error(traceback.format_exc())
                
                                         
                if(os.path.exists(path_catalog_pans)):
                    try:
                        matchtable_pans=match_ps2_50cm(source=source[j],path_cat=path_catalog_pans)
                        matchtable_pans.write(filaname[:-5] +'_ps2.fits',format='fits' ,overwrite=True)
                        c_match,medra,stdra,meddec,stddec =xmach_fig(filaname[:-5] +'_ps2.fits',figpath,stype='ps2')
                        hdul = fits.open(filaname)
                        iheader1 = hdul[1].header
                        iheader1['ASTPRA'] = stdra
                        iheader1['ASTPDEC'] = stddec
                        iheader1['ASTMRA'] = medra
                        iheader1['ASTMDEC'] = meddec
                        iheader1['ASTCNT'] = c_match
                        fits.writeto(filaname[:-5] +'_ps2.fits',data=matchtable_pans,header=iheader1,overwrite=True)
                     
			
                        
                        
                    except Exception as e:
                        logger.error(traceback.format_exc())
                        continue  
                
                
                if(os.path.exists(path_catalog_wise)):
                    try:
                        matchtable_wise=match_wise_50cm(source=source[j],path_cat=path_catalog_wise)
                        matchtable_wise.write(filaname[:-5] +'_wise.fits',format='fits' ,overwrite=True)
                        c_match,medra,stdra,meddec,stddec =xmach_fig(filaname[:-5] +'_wise.fits',figpath,stype='wise')
                        hdul = fits.open(filaname)
                        iheader1 = hdul[1].header
                        iheader1['ASTPRA'] = stdra
                        iheader1['ASTPDEC'] = stddec
                        iheader1['ASTMRA'] = medra
                        iheader1['ASTMDEC'] = meddec
                        iheader1['ASTCNT'] = c_match
                        
                        fits.writeto(filaname[:-5] +'_wise.fits',data=matchtable_wise,header=iheader1,overwrite=True)
                    except Exception as e:
                        logger.error(traceback.format_exc())
                        continue  


 

def pro(ttfname,date):
    rootdir=get_current_dir()
    try:
        tid,objid,filterid=ttfname.split('_')
    except:
        tid,obj1,obj2,filterid=ttfname.split('_')
        objid = obj1+'_'+obj2
    rawpath=rootdir+'reception/'+str(date)+'/raw/'
    fileguide_raw_path=rawpath+'fileguide/'
    calpath=rootdir+'reception/'+str(date)+'/cal/'
    scipath=rootdir+'reception/'+str(date)+'/sci/'+ttfname+'/'
    
    fileguide_sci_path=scipath+'fileguide/'
    fileguide_raw_path=rawpath+'fileguide/'
    
    figsave_path = '/home/50cm/dugking/CrossPic/'+date+'/'+tid+'/'
    logdir = rootdir+'run_pipe_log/'+str(date)+'/'+ttfname+'/'

    sexcatlist = glob.glob(scipath+'*sexcat.fits')
    sublist = os.path.join(scipath, date+'T*')
 
    subdirs = [d for d in glob.glob(sublist) if os.path.isdir(d)]
    logger.info(subdirs)
    #file = glob.glob(fileguide_raw_path+'sci_'+tid+'_*'+objid+'_'+filterid+'.npy')
    #logger.info(fileguide_raw_path+'sci_'+tid+'_*'+objid+'_'+filterid+'.npy')
    #logger.info(file[0])
    
    if(subdirs):
        for j in range(0,len(subdirs)):
                refpath='/home/50cm/data2/'
                 
                match_pro(date,subdirs[j]+'/',refpath,figsave_path)
                logger.info(subdirs[j])
	
        # # filelist=np.load([0])
        # # logger.info(filelist)
        #     filelist= glob.glob(subdirs+'*sexcat.fits')
        #     nfile=len(filelist)
        #     for i in range(0,nfile):
        #         # dirname,filename=os.path.split(filelist[i])
        #         # subpath = datename(filelist[i]) 
				
        #         # tscidir = scipath +  subpath +'/'
    elif(sexcatlist):
        nfile=len(sexcatlist)
        refpath='/home/50cm/data2/'
                 
        match_pro(date,scipath,refpath,figsave_path)
 
    # if not(os.path.exists(scipath+'fileguide/'+ttfname+'_sciimg.npy')):
    #     logger.info("There is no such observation: %s" % (ttfname)) 
    #     return None
    # filelist=np.load(fileguide_raw_path+'sci_'+ttfname+'.npy')
    # #filelist= glob.glob(rawinpath+tid+'_'+target+'_'+filterid+'_'+'*.fit')
    # nfile=len(filelist)
    # for i in range(0,nfile):
    #     dirname,filename=os.path.split(filelist[i])
    #     subpath = datename(filelist[i]) 
    #     tscidir = scipath +  subpath +'/'
    #     refpath='/home/50cm/data2/'
    #     match_pro(date,tscidir,refpath,)
 

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()  # 创建parser
    # parser.add_argument('--date', type=str, default='20230313', help='输入处理日期')
    # args = parser.parse_args()  # 参数解析
    # loggerloguru = Logger(args.date, '', "_yrecp").get_logger
    import argparse
    import time
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20230207', help='输入处理日期')
    parser.add_argument('--tid_target_filter', default='y50a_H04HZ2_mr',type=str, help='')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "ymatch").get_logger
    starttime=time.time()
    pro(args.tid_target_filter,args.date)
    endtime=time.time()
      
    logger.info("[_ycom] "+"time spending:{}s".format(endtime-starttime))