# -*- coding: UTF-8 -*-
import concurrent
import ccdproc
import logging
from astropy.stats import sigma_clipped_stats
import time as ts
import os
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from loguru import logger as loggerevaluate, logger
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
import config

from concurrent.futures import ThreadPoolExecutor

from astropy.time import Time
from astropy import coordinates as coord, units as u
 
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
# from PyAstronomy import pyas
import glob

import math
import time
import datetime as dt
 
#from astropy.time.core import TimezoneInfo
from astropy.coordinates import EarthLocation
from astropy import units as u

from datetime import datetime,timezone,timedelta
from pytz import timezone

def pk_2_utc(pktime_str:str)->datetime:
     
    # 构造出没有时区的datetime对象
    naive_time = datetime.strptime(pktime_str,'%Y-%m-%dT%H:%M:%S')
    #datetime.strptime('2023-05-23T21:39:45', '%Y-%m-%dT%H:%M:%S').replace(tzinfo=beijing_tz)
    # 将上述对象转化为时区为Asia/Shanghai的datetime对象
    pktime = naive_time.astimezone(timezone('Asia/Shanghai'))
    # 将时区为上海的datetime对象转化为时区为utc的时间对象
    utctime = pktime.astimezone(timezone('UTC'))
    utc_time_format = utctime.strftime("%Y-%m-%dT%H:%M:%S")
    return utc_time_format

def ams_cal(sciimg_name,ira,idec):
    # res=fits.open(filename)
    # data = res[1].data
    
    hdr = fits.getheader(sciimg_name)
    date_obs = hdr['DATE_OBS']
    obstime=pk_2_utc(date_obs)
    expt     = float(hdr['EXPOSURE'])
    amss, jd_mid, tdb_mid = mairmass.cam(obstime, expt, ira, idec)
 
    return amss, jd_mid, tdb_mid


def TransLDAC_time(sciimg_name):
    hdr=fits.getheader(sciimg_name) 
    date_obs = hdr['DATE']
    obstime=pk_2_utc(date_obs)
    exptime = float(hdr['EXPOSURE'])
    obstime=datetime.strptime(obstime,'%Y-%m-%dT%H:%M:%S')
  
    # 计算偏移量
    offset = timedelta(seconds=float(exptime)/2)
    # 获取修改后的时间并格式化
    timeplus = obstime + offset
    timeplus=str(timeplus).replace(' ','T')
    t=Time(timeplus)
    tjd=Time(np.array(str(t.jd)),format='jd')
    tb=tjd.tdb
    return str(obstime),str(timeplus),str(tb)

def getindex_location(tra,tdec,ralist,declist):
    tra = float(tra)
    tdec = float(tdec)
    #print(tra,tdec)
    declist=declist.astype('float64')
    ralist=ralist.astype('float64')
    for i in range(len(ralist)):
       
        try:
            d=np.sqrt((ralist[i]-tra)**2+(declist[i]-tdec)**2)
            
        except Exception as e:
            logger.info(e)
            continue
 
        if(d<8/3600):
            return i,round(ralist[i],4),round(declist[i],4)
    return 0

def mag_reg_target(filename,tra,tdec,targetname,regname):
     
    regListn   = regname #ildac_b[:-5]+ "_select.reg"
    regList    = open(regListn,"w")
    regList.write('# Region file for DS9'+"\n")
    # Region file format: DS9 version 4.1
 
    #circle(14:03:38.564, +54:18:42.02,6.000") # width=4 font="helvetica 8 bold roman" text={ngc5457_sn}
    regList.write('global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n")
    regList.write('fk5'+"\n")
    sline='circle('+str(tra)+','+ str(tdec)+', 10.00)'  + '# width=4 font="helvetica 8 bold roman" text={ngc5457_sn}#+'
    #  font="helvetica 8 normal" text={' + 'imag:'+str(round(mag[i],2))+'}'
        #logger.info(sline)
    regList.write(sline+"\n")
    regList.close()
    logger.info('the reg is ok')


#def mag_targetlist(listpoint,regname)

def mag_target(sciimg_name,tra,tdec):
    tra= float(tra)
    tdec = float(tdec)
    hdr_sci=fits.getheader(sciimg_name) 
    date_obs = hdr_sci['DATE']
    obstime =pk_2_utc(date_obs)
    exptime = float(hdr_sci['EXPOSURE'])
    tra = float(tra)
    tdec = float(tdec)
    #obstime=datetime.datetime.strptime(date_obs,'%Y-%m-%dT%H:%M:%S')
    ldac_name= sciimg_name.split('.fits')[0]+'_sexcat.fits'
    ldac_aper_name = sciimg_name.split('.fits')[0]+'_sexcat_apcorr.fits'
    
    mag_auto =0.0
    magerr_auto =0.0
    mag_psf =0.0
    magerr_psf =0.0
    mag_aper =0.0
    magerr_aper =0.0
    amss =0.0
    jd_mid =0
    tdb_mid =0
    ximage =0
    yimage =0
    #index,ira,idec=getindex_location(tra,tdec,ra, dec)
    amss, jd_mid, tdb_mid = ams_cal(sciimg_name,tra,tdec)
    obstime,utctime,tdbtime=TransLDAC_time(sciimg_name)
    filename = os.path.basename(sciimg_name)
    fileid = filename.split('.fits')[0] 
    iflag = 0
    if(os.path.exists(ldac_aper_name)):
        try:
            logger.info(ldac_aper_name)
            ildac = fits.getdata(ldac_aper_name,ext=2)
            #mag_aper,magerr_aper=ildac["MAG_FAPER"],ildac["MAGERR_FAPER"]
            ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
            lflag = ildac["FLAGS"]
            lmag_auto, lmagerr_auto=ildac["MAG_AUTO_S"],ildac["MAGERR_AUTO_S"]
            lmag_aper,lmagerr_aper=ildac["MAG_FAPER_S"],ildac["MAGERR_FAPER_S"]
            lmag_psf,lmagerr_psf=ildac["MAG_PSF_S"],ildac["MAGERR_PSF_S"]
            index,ira,idec=getindex_location(tra,tdec,ra, dec)
            mag_auto=lmag_auto[index]
            magerr_auto=lmagerr_auto[index]
            mag_psf=lmag_psf[index]
            magerr_psf=lmagerr_psf[index]
            mag_aper=lmag_aper[index]
            magerr_aper=lmagerr_aper[index]
            ximage,yimage= ildac["X_IMAGE"][index],ildac["Y_IMAGE"][index]
            iflag =lflag[index]
        except:
            return 0
    elif (os.path.exists(ldac_name)):
        try:
            logger.info('there is no aper_cor mag')
            logger.info(ldac_name)
            ildac = fits.getdata(ldac_name,ext=2)
             
            #mag_aper,magerr_aper=ildac["MAG_FAPER"],ildac["MAGERR_FAPER"]
            ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
            lflag = ildac["FLAGS"]
            lmag_auto, lmagerr_auto=ildac["MAG_AUTO_S"],ildac["MAGERR_AUTO_S"]
            lmag_psf,lmagerr_psf=ildac["MAG_PSF_S"],ildac["MAGERR_PSF_S"]
             
            index,ira,idec=getindex_location(tra,tdec,ra, dec)
            mag_auto=lmag_auto[index]
            magerr_auto=lmagerr_auto[index]
            mag_psf=lmag_psf[index]
            magerr_psf=lmagerr_psf[index]
            ximage,yimage= ildac["X_IMAGE"][index],ildac["Y_IMAGE"][index]
            iflag =lflag[index]
        except:
            return 0
    else: 
            logger.info('there is no any mag')
    mag_dict={'fileid':fileid,'mag_auto':mag_auto,'magerr_auto':magerr_auto,'mag_psf':mag_psf,'magerr_psf':magerr_psf,'mag_aper':mag_aper,'magerr_aper':magerr_aper,'amss':amss,'jd':jd_mid,'tdb':tdbtime,'utctime':utctime,'exptime':exptime,'ximage':ximage,'yimage':yimage,'ra':ira,'dec':idec,'iflag':iflag}
    return mag_dict

#写到这里，要加本地时间转换成ut时间的程序
def lightcurve_txt(savepath,listpath,tra,tdec):
    filelen = len(listpath)
    data_array= []
    for i in range(0,filelen):
        sciimg_name =listpath[i]
        
        mdict = mag_target(sciimg_name,tra,tdec)
        if (mdict ==0):
            continue
        data_array.append([mdict['fileid'],mdict['utctime'],mdict['mag_auto'],mdict['magerr_auto'],mdict['mag_psf'],mdict['magerr_psf'],mdict['mag_aper'],
                           mdict['magerr_aper'],mdict['amss'],mdict['jd'],mdict['tdb'],mdict['exptime'],mdict['ximage'],mdict['yimage'],mdict['ra'],mdict['dec'],mdict['iflag']])
        #data_array.append(np.array(mdict))
     
    txtheader = 'fileid,utctime,mag_auto,magerr_auto,mag_psf,magerr_psf,mag_aper,magerr_aper,airmass,JD,TDB,exptime,ximage,yimage,ra,dec,flag'
    
    t_light=Table(rows=data_array,names=['fileid','utctime','mag_auto','magerr_auto','mag_psf','magerr_psf','mag_aper','magerr_aper','airmass','JD','TDB','exptime','ximage','yimage','ra','dec','flag'])
    t_light.sort('utctime')
    
    #txtname=savepath+date+'_'+ttfname+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_l.txt'
    txtname=savepath+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_l.txt'
    tablename  = savepath+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_l.csv'
    np.savetxt(txtname, np.array(data_array), fmt="%s", delimiter=' ',header = txtheader)
    t_light.write(tablename,overwrite=True)
    return data_array,t_light
        

def delta_date():
    start_date = dt.date(2023, 5, 20)
    end_date =dt.date(2023, 7, 2)
    delta = dt.timedelta(days=1)

    dates = []
    while start_date < end_date:
        dates.append(start_date.strftime('%Y%m%d'))
        start_date += delta
    logger.info(dates)
    return dates

def getfilelist(datelist,ttfname):
    sci_filelist=[]
    rootdir = config["dirPath"]["Origin_Raw_Path"]
    for i in range(0,len(datelist)):
        date=datelist[i]
        tid,objid,filterid=ttfname_split(ttfname)
  
        scipath=rootdir+'reception/'+str(date)+'/sci/'+ttfname+'/'
          
        sublist = os.path.join(scipath, date+'T*')
        subdirs = [d for d in glob.glob(sublist) if os.path.isdir(d)]
        #logger.info(subdirs)
        if(subdirs):
            for j in range(0,len(subdirs)):
                source=glob.glob(subdirs[j]+'/'+'*sexcat.fits')
                if source:
                    sciname = source[0].split('_sexcat.fits')[0]+'.fits'
                    sci_filelist.append(sciname)
                   
        logger.info(sci_filelist) 
        #logger.info(sci_filelist.sorted())
    return sci_filelist
    
    

def pro(ttfname,tra,tdec):
    #dates=delta_date()
    datelist = ['20230527']
    
    getfilelist(datelist,ttfname)
    mag_target(sciimg_name,tra,tdec)
    
     
#def target_ref(targetlist,reflist):
     

if __name__ == "__main__":
    # sciimg_name = '/home/mpilot/1m6/reception/20230527/sci/mb_ngc5457sn_v/20230527T191942/mb_sc_tngc5457sn_v_20230527191941_230_sciimg.fits'
    # #RA,DEC: 210.9106833, +54.3116718
    
    # dict = mag_target(sciimg_name,'210.9106833','+54.3116718')
    # print(dict)
    dates=delta_date()
    #sci_list = getfilelist(dates,'my_ngc5457sn_r')
    # for i in range(len(sci_list)):
    #     dict = mag_target(sci_list[i],210.9106833, +54.3116718)
    #     print(dict)
    # circle(14:03:38.5640,+54:18:42.020,6.000") # color=blue width=4 font="helvetica 8 bold roman" text={ngc5457_sn}
    # circle(14:03:45.1544,+54:16:15.341,8.666") # text={ref}
    # circle(14:04:02.8278,+54:18:18.226,9.839") # text={check}
    
    
    sci_list1 = getfilelist(dates,'mb_ngc5457sn_u')
    sci_list2 = getfilelist(dates,'mb_ngc5457sn_v')
    sci_list3 = getfilelist(dates,'my_ngc5457sn_g')
    sci_list4 = getfilelist(dates,'my_ngc5457sn_r')
    
    
    savepath1 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_u_check'
    savepath2 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_v_check'
    savepath3 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_g_check'
    savepath4 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_r_check'
    
    # savepath1 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_u_ref'
    # savepath2 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_v_ref'
    # savepath3 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_g_ref'
    # savepath4 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_r_ref'
    
    
    # savepath1 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_u'
    # savepath2 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_v'
    # savepath3 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_g'
    # savepath4 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_r'
    #mkdir('/home/mpilot/1m6/reception/lightcurve/')
    tra0,tdec0=210.9106833, +54.3116718
    tra1,tdec1 = d2hms('14:03:45.1544', '+54:16:15.341', conv=1)
    tra2,tdec2 = d2hms('14:04:02.8278','+54:18:18.226', conv=1)
    #print(tra1,tdec1)
    
    array1,table1 =lightcurve_txt(savepath1,sci_list1,tra2,tdec2)
    array2,table2 =lightcurve_txt(savepath2,sci_list2,tra2,tdec2)
    array3,table3 =lightcurve_txt(savepath3,sci_list3,tra2,tdec2)
    array4,table4 =lightcurve_txt(savepath4,sci_list4,tra2,tdec2)
    
    #print(array1)
    print(table1)
    #print(table1.colnames)
