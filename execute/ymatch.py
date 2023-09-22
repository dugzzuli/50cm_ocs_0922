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

def match_gaiadr3_50cm(source,path_cat,match_r = 3.0):
	if not(os.path.exists(path_cat)):
		print('there is no such directory')
		return None
	name_catalog='gaiadr3'
	#Read the center position from header
	header=fits.getheader(source[:-21]+'.fits')
	ra_cen=header['OBJCTRA']
	dec_cen=header['OBJCTDEC']
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
		#print(i,len(idxca))
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
		tabler.add_column(offset,name='offset')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,keys='NUMBER')  
		deltara=np.array(newtable['ALPHA_J2000'])*3600-np.array(newtable['ra'])*3600
		deltadec=np.array(newtable['DELTA_J2000'])*3600-np.array(newtable['dec'])*3600
		newtable.add_column(deltara,name='deltara')
		newtable.add_column(deltadec,name='deltadec')
		print('Original Num:',len(data),', Xmatch with gaiadr3 Num:',len(tablel))		
		return newtable

def match_ps2_50cm(source,path_cat,match_r = 3.0):
	name_catalog='ps2'
	
	#Read the center position from header
	header=fits.getheader(source[:-21]+'.fits')
	ra_cen=header['OBJCTRA']
	dec_cen=header['OBJCTDEC']
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
		#print(i,len(idxca))
		#cattablei=pd.read_table(path_cat+catalogname[i],header=1000,sep=',')
		#cattablei=Table.from_pandas(cattablei) 
		cattablei=Table.read(path_cat+catalogname[i])
		rai=np.array(cattablei['raMean']).reshape(1,-1)[0]
		deci=np.array(cattablei['decMean']).reshape(1,-1)[0]    
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
		print('No match catalog data')
		return None
	else:
		ra2=np.array(catalog['raMean']).reshape(1,-1)[0]
		dec2=np.array(catalog['decMean']).reshape(1,-1)[0]
		clog=SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
		start2=time.time()
		max_sep = 5.0 
		#region=c0.separation(clog)<40*u.arcminute
		#idxd=np.where(region==True)[0]
		 
		idx, d2d, d3d = c.match_to_catalog_sky(clog) 

		index=np.where(d2d.arcsec<match_r)[0]

		tablel=data[index]
		tabler=catalog[idx[index]]
		if len(tablel)==0 or len(tabler)==0 :
		
			print('No match data')
			return None
		time3=time.time()
		number=tablel['NUMBER']
		offset=d2d.arcsec[index]
		tabler.add_column(offset,name='offset')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,keys='NUMBER')  
		deltara=np.array(newtable['ALPHA_J2000'])*3600-np.array(newtable['raMean']).reshape(1,-1)[0]*3600
		deltadec=np.array(newtable['DELTA_J2000'])*3600-np.array(newtable['decMean']).reshape(1,-1)[0]*3600
		newtable.add_column(deltara,name='deltara')
		newtable.add_column(deltadec,name='deltadec')
		time4=time.time()
		print('Original Num:',len(data),', Xmatch with ps2 Num:',len(tablel))	
		print(time1,time2,time3,time4)	
		return newtable,catalog

def match_2mass_50cm(source,path_cat,match_r = 3.0):
	name_catalog='2mass'
	#path_catalog='/media/huohe/hlserver/'
	#catalog=glob.glob(path_catalog+'gaiasource/'+'*.gz')
	#path_file='/home/50cm/50cm_scripts/50cm_pro/reception/'+date+'/sci/'+source+'/'
	#filename=glob.glob(path_file+'*subbkg_sexcat.fits')
	#Read the center position from header
	header=fits.getheader(source[:-21]+'.fits')
	ra_cen=header['OBJCTRA']
	dec_cen=header['OBJCTDEC']
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
		#print(i,len(idxca))
		#cattablei=pd.read_table(path_cat+catalogname[i],header=1000,sep=',')
		#cattablei=Table.from_pandas(cattablei) 
		cattablei=Table.read(path_cat+catalogname[i])
		rai=np.array(cattablei['TMASS_RA'])
		deci=np.array(cattablei['TMASS_DEC'])       
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
		ra2=np.array(catalog['TMASS_RA'])
		dec2=np.array(catalog['TMASS_DEC'])      
		clog=SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
		start2=time.time()
		max_sep = 5.0 
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
		tabler.add_column(offset,name='offset')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,keys='NUMBER')  
		deltara=np.array(newtable['ALPHA_J2000'])*3600-np.array(newtable['TMASS_RA'])*3600
		deltadec=np.array(newtable['DELTA_J2000'])*3600-np.array(newtable['TMASS_DEC'])*3600
		newtable.add_column(deltara,name='deltara')
		newtable.add_column(deltadec,name='deltadec')
		print('Original Num:',len(data),', Xmatch with 2mass Num:',len(tablel))		
		return newtable

def match_wise_50cm(source,path_cat,match_r = 3.0):
	name_catalog='wise'
	
	#Read the center position from header
	header=fits.getheader(source[:-21]+'.fits')
	ra_cen=header['OBJCTRA']
	dec_cen=header['OBJCTDEC']
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
		#print(i,len(idxca))
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
		tabler.add_column(offset,name='offset')
		tabler.add_column(number,name='NUMBER')
		newtable=join(tablel,tabler,keys='NUMBER')  
		deltara=np.array(newtable['ALPHA_J2000'])*3600-np.array(newtable['ra'])*3600
		deltadec=np.array(newtable['DELTA_J2000'])*3600-np.array(newtable['dec'])*3600
		newtable.add_column(deltara,name='deltara')
		newtable.add_column(deltadec,name='deltadec')
		print('Original Num:',len(data),', Xmatch with wise Num:',len(tablel))		
		return newtable

def match_pro(date,tcsipath,refpath):

    #filesavedir  = '/home/guohl/xmatch/xmatch50cm_%s/'%obsDate
    #if not os.path.exists(filesavedir): os.mkdir(filesavedir)
    #filedir='/data2/workspace/50cm_pro/reception/'+obsDate+'/sci/'
    #objects=glob.glob(filedir+'*')
    source=glob.glob(tcsipath+'*sexcat.fits')
    path_catalog_gaia = refpath + 'GaiaDR3FITS/'
    path_catalog_pans = refpath + 'Panstarrs2FITS/'
    path_catalog_2mass = refpath + '2massFITS/'
    path_catalog_wise = refpath + 'allwiseFITS/'
    if len(source)>0:  
        for j in range(len(source)):
            #try:
            if(1==1):
                filaname = source[j]
                matchtable_gaia=match_gaiadr3_50cm(source=source[j],path_cat=path_catalog_gaia)
                matchtable_gaia.write(tcsipath+filaname[:-5] +'_gaia.fits',format='fits' ,overwrite=True)
                logger.info(f'{filaname[:-5]}_gaia.fits')
                matchtable_pans=match_ps2_50cm(source=source[j],path_cat=path_catalog_pans)
                matchtable_pans.write(tcsipath+filaname[:-5] +'_ps2.fits',format='fits' ,overwrite=True)
                
                matchtable_2mass=match_2mass_50cm(source=source[j],path_cat=path_catalog_2mass)
                matchtable_2mass.write(tcsipath+filaname[:-5] +'_2mass.fits',format='fits' ,overwrite=True)
		
                matchtable_wise=match_wise_50cm(source=source[j],path_cat=path_catalog_wise)
                matchtable_wise.write(tcsipath+filaname[:-5] +'_wise.fits',format='fits' ,overwrite=True)


def pro(ttfname,date):
    rootdir=get_current_dir()
    rawpath=rootdir+'reception/'+str(date)+'/raw/'
    fileguide_raw_path=rawpath+'fileguide/'
    calpath=rootdir+'reception/'+str(date)+'/cal/'
    scipath=rootdir+'reception/'+str(date)+'/sci/'+ttfname+'/'
    
    fileguide_sci_path=scipath+'fileguide/'
     
    logdir = rootdir+'run_pipe_log/'+str(date)+'/'+ttfname+'/'
 
    if not(os.path.exists(scipath+'fileguide/'+ttfname+'_sciimg.npy')):
        logger.info("There is no such observation: %s" % (ttfname)) 
        return None
    filelist=np.load(fileguide_raw_path+'sci_'+ttfname+'.npy')
    #filelist= glob.glob(rawinpath+tid+'_'+target+'_'+filterid+'_'+'*.fit')
    nfile=len(filelist)
    for i in range(0,nfile):
        dirname,filename=os.path.split(filelist[i])
        subpath = datename(filelist[i]) 
        tscidir = scipath +  subpath +'/'
        refpath='/home/50cm/data2/'
        match_pro(date,tscidir,refpath)


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