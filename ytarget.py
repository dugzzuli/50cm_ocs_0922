# -*- coding: UTF-8 -*-
import concurrent
import logging
from astropy.stats import sigma_clipped_stats,sigma_clip
import time as ts
import os
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from loguru import logger as loggerevaluate, logger
from config import config
from lib.phot.ybias import get_current_dir
from concurrent.futures import ThreadPoolExecutor

from astropy.time import Time
from astropy import coordinates as coord, units as u
 
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
# from PyAstronomy import pyas
import glob
import math
import time
import datetime as dt
 
from execute.flatdiff import mkdir 
from  lib.phot import cal_airmass  
from execute.yphotutils import  datename,d2hms,crossmatch,ttfname_header,ttfname_split

 
#from astropy.time.core import TimezoneInfo
from astropy.coordinates import EarthLocation
from astropy import units as u

from datetime import datetime,timezone,timedelta
from pytz import timezone
from tables.atom import Float32Atom
import csv
 

def ams_cal(sciimg_name,ira,idec):
    # res=fits.open(filename)
    # data = res[1].data
    
    hdr = fits.getheader(sciimg_name)
    date_obs = hdr['DATE-OBS']
    #DATE-OBS= '2023-05-27T15:21:27' /YYYY-MM-DDThh:mm:ss observation start, UT   
    
    obstime = date_obs
    # print('obstime,ira,idec')
    # print(obstime,ira,idec)
    expt     = float(hdr['EXPTIME'])
    amss, jd_mid, tdb_mid = cal_airmass.cam(obstime, expt, ira, idec)
 
    return amss, jd_mid, tdb_mid


def TransLDAC_time(sciimg_name):
    hdr=fits.getheader(sciimg_name) 
    date_obs = hdr['DATE-OBS']
    #obstime=pk_2_utc(date_obs)
    obstime = date_obs
    exptime = float(hdr['EXPTIME'])
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
    date_obs = hdr_sci['DATE-OBS']
    obstime = date_obs#pk_2_utc(date_obs)
    exptime = float(hdr_sci['EXPTIME'])
     
    #obstime=datetime.datetime.strptime(date_obs,'%Y-%m-%dT%H:%M:%S')
    ldac_name= sciimg_name.split('_sciimg.fits')[0]+'_sciimg_sexcat.fits'
    print(ldac_name)
    ldac_aper_name = sciimg_name.split('_sciimg.fits')[0]+'_sexcat_apcorr.fits'
    
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
    # print('amss, jd_mid, tdb_mid')
    # print(amss)
    obstime,utctime,tdbtime=TransLDAC_time(sciimg_name)
    filename = os.path.basename(sciimg_name)
    fileid = filename.split('.fits')[0] 
    iflag = 0
    if(os.path.exists(ldac_name)):
        try:
            logger.info('there is no aper_cor mag')
            logger.info(ldac_name)
            ildac = fits.getdata(ldac_name,ext=2)
            lfwhm=ildac['FWHM_IMAGE']
             
            #mag_aper,magerr_aper=ildac["MAG_FAPER"],ildac["MAGERR_FAPER"]
            ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
            lflag = ildac["FLAGS"]
            lmag_auto, lmagerr_auto=ildac["MAG_AUTO_S"],ildac["FLUXERR_AUTO"]/ildac["FLUX_AUTO"]#ildac["MAGERR_AUTO"]/1.0857*exptime
            
            
            lmag_psf,lmagerr_psf=ildac["MAG_PSF_S"],ildac["FLUXERR_PSF"]/ildac["FLUX_PSF"]
            #lmagerr_psf =  1.0857*(np.array(ildac["FLUXERR_PSF"],dtype=float)/np.array(ildac["FLUX_PSF"],dtype=float) )
            lk_radius=ildac["KRON_RADIUS"]
            llbkg = ildac['BACKGROUND']
            lsnr = ildac["SNR_WIN"]
            index,ira,idec=getindex_location(tra,tdec,ra, dec)
            mag_auto=lmag_auto[index]
            magerr_auto=lmagerr_auto[index]
            mag_psf=lmag_psf[index]
            magerr_psf=lmagerr_psf[index]
            ximage,yimage= ildac["X_IMAGE"][index],ildac["Y_IMAGE"][index]
            iflag =lflag[index]
            ifwhm =lfwhm[index]
            ibkg=llbkg[index]
            isnr = lsnr[index]
            ikradius = lk_radius[index]
            print(ifwhm )

        except:
            return 0
    elif (os.path.exists(ldac_aper_name)):
        try:
            logger.info(ldac_name)
            ildac = fits.getdata(ldac_aper_name,ext=2)
            lfwhm=ildac['FWHM_IMAGE']
             
            #mag_aper,magerr_aper=ildac["MAG_FAPER"],ildac["MAGERR_FAPER"]
            ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
            lflag = ildac["FLAGS"]
            lmag_auto, lmagerr_auto=ildac["MAG_AUTO_S"],ildac["FLUXERR_AUTO"]/ildac["FLUX_AUTO"]*1.0857 #ildac["MAGERR_AUTO_S"]/1.0857*exptime
          
            #lmagerr_auto = 1.0857*(np.array(ildac["FLUXERR_AUTO"],dtype=float)/np.array(ildac["FLUX_AUTO"],dtype=float) )
            lmag_aper,lmagerr_aper=ildac["MAG_FAPER_S"],ildac["FLUXERR_FAPER"]/ildac["FLUX_FAPER"] #ildac["MAGERR_FAPER"]/1.0857*exptime
            #lmagerr_aper =1.0857*(ildac["FLUXERR_FAPER"]/ildac["FLUX_FAPER"]) 
            lmag_psf,lmagerr_psf=ildac["MAG_PSF_S"],ildac["FLUXERR_PSF"]/ildac["FLUX_PSF"]*1.0857#ildac["MAGERR_PSF"]/1.0857*exptime
            #lmagerr_psf =  1.0857*(np.array(ildac["FLUXERR_PSF"],dtype=float)/np.array(ildac["FLUX_PSF"],dtype=float) )
            lk_radius=ildac["KRON_RADIUS"]
            llbkg = ildac['BACKGROUND']
            lsnr = ildac["SNR_WIN"]
            index,ira,idec=getindex_location(tra,tdec,ra, dec)
            mag_auto=lmag_auto[index]
            magerr_auto=lmagerr_auto[index]
            mag_psf=lmag_psf[index]
            magerr_psf=lmagerr_psf[index]
            mag_aper=lmag_aper[index]
            magerr_aper=lmagerr_aper[index]
            ximage,yimage= ildac["X_IMAGE"][index],ildac["Y_IMAGE"][index]
            iflag =lflag[index]
            ifwhm =lfwhm[index]
            ibkg=llbkg[index]
            isnr = lsnr[index]
            ikradius = lk_radius[index]
            
            print(ibkg)
            

        except:
            return 0
    else: 
            logger.info('there is no any mag')
    mag_dict = 0
    if(ifwhm<9):
        mag_dict={'fileid':fileid,'mag_auto':mag_auto,'magerr_auto':magerr_auto,'mag_psf':mag_psf,'magerr_psf':magerr_psf,'mag_aper':mag_aper,'magerr_aper':magerr_aper,'amss':amss,'jd':jd_mid,'tdb':tdbtime,'utctime':utctime,'exptime':exptime,'ximage':ximage,'yimage':yimage,'ra':ira,'dec':idec,'iflag':iflag,'fwhm':ifwhm,'snr':isnr,'kron_radius':ikradius,'lbkg':ibkg}
    return mag_dict


def lightcurve_txt(savepath,listpath,tra,tdec):
    filelen = len(listpath)
    data_array= []
    for i in range(0,filelen):
        sciimg_name =listpath[i]
        
        mdict = mag_target(sciimg_name,tra,tdec)
        if (mdict ==0):
            continue
        print(mdict)
        data_array.append([mdict['fileid'],mdict['utctime'],mdict['mag_auto'],mdict['magerr_auto'],mdict['mag_psf'],mdict['magerr_psf'],mdict['mag_aper'],
                           mdict['magerr_aper'],mdict['amss'],mdict['jd'],mdict['tdb'],mdict['exptime'],mdict['ximage'],mdict['yimage'],mdict['ra'],mdict['dec'],mdict['iflag'],mdict['fwhm'],mdict['snr'],mdict['kron_radius'],mdict['lbkg']])
        #data_array.append(np.array(mdict))

    txtheader = 'fileid,utctime,mag_auto,magerr_auto,mag_psf,magerr_psf,mag_aper,magerr_aper,airmass,JD,TDB,exptime,ximage,yimage,ra,dec,flag,fwhm,snr,kron_radius,lbkg'
    
    t_light=Table(rows=data_array,names=['fileid','utctime','mag_auto','magerr_auto','mag_psf','magerr_psf','mag_aper','magerr_aper','airmass','JD','TDB','exptime','ximage','yimage','ra','dec','flag','fwhm','snr','kron_radius','lbkg'])
     
    t_light.sort('utctime')
    
    #txtname=savepath+date+'_'+ttfname+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_l.txt'
    #txtname=savepath+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_l.txt'
    tablename  = savepath+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_l.csv'
    #np.savetxt(txtname, np.array(data_array), fmt="%s", delimiter=' ',header = txtheader)
    t_light.write(tablename,overwrite=True)
    return data_array,t_light
            

def delta_date():
    start_date = dt.date(2023, 5, 20)
    end_date =dt.date(2023, 7,1 )
    delta = dt.timedelta(days=1)

    dates = []
    while start_date < end_date:
        dates.append(start_date.strftime('%Y%m%d'))
        start_date += delta
    logger.info(dates)
    return dates

def getfilelist(datelist,ttfname):
    sci_filelist=[]
    #rootdir = config["dirPath"]["Origin_Raw_Path"]
    rootdir = get_current_dir()
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


def flux_alignment(kdict,filelist,ref_file):
    # tra= float(tra)
    # tdec = float(tdec)
    hdr_sci=fits.getheader(ref_file) 
    date_obs = hdr_sci['DATE-OBS']
    obstime = date_obs 
    exptime_ref = float(hdr_sci['EXPTIME'])
    #tid,target,filterid = ttfname.split('_')

    tid,target,filterid = ttfname_header(ref_file)
    k1 = kdict[filterid]
    fsList = {}
    
    for i in range (0,len(filelist)):
         
        iimg    = filelist[i]
        iimgKey = filelist[i] #iimg.split("/")[-1]
        logger.info(iimgKey)
        ihdr  = fits.getheader(iimg)
        ifluxCorr = 10**(0.4 * float(ihdr["AIRMASS"]) * k1)  # "u": 0.7639, "v": 0.5411, "g": 0.2482, "r": 0.1846, "i": 0.0087, "z": 0.0056
        fsList[iimgKey] = ifluxCorr

    logger.info(fsList)
    refImg = ref_file[:-5] + "_sexcat.fits"
    refImg = ref_file.split('_sciimg.fits')[0] + "_sciimg_sexcat.fits"
    logger.info(refImg)
    
    refLdac  = fits.getdata(refImg, ext=2)
    flagRef  = refLdac["FLAGS"]
    snrRef   = refLdac["SNR_WIN"]
    #fluxRef  = refLdac["FLUX_PSF"]
    fluxRef  = refLdac["FLUX_AUTO_S"]
     
    gid      = (flagRef==0) & (snrRef>50) & (fluxRef>0.0)
    fluxRef  = fluxRef[gid] * fsList[ref_file]
    
    raRef    = refLdac["ALPHA_J2000"][gid]
    decRef   = refLdac["DELTA_J2000"][gid]
    
    for iimgKey in fsList.keys():
        if iimgKey == refImg: continue
    #if iimgKey != "ngc157.u.exp5.sub.fits": continue
        ifs    = fsList[iimgKey]
        ildacn = iimgKey[:-5] + "_sexcat.fits"
        logger.info(ildacn)
        ildac  = fits.getdata(ildacn, ext=2)
        iflag  = ildac["FLAGS"]
        isnr   = ildac["SNR_WIN"]
        #iflux  = ildac["FLUX_PSF"]
        iflux  = ildac["FLUX_AUTO_S"]
        gid    = (iflag==0) & (isnr>50) & (iflux>0.0)
        iflux  = iflux[gid] * ifs
        ira    = ildac["ALPHA_J2000"][gid]
        idec   = ildac["DELTA_J2000"][gid]
        logger.info("    Total %d stars are selected in Image %s"%(len(ira), iimgKey))
        
        # match the catalgs
        idRef, iid = crossmatch(raRef, decRef, ira, idec, aperture=2.0)
        irflux     = fluxRef[idRef]/iflux[iid]
        irfluxNew  = sigma_clip(irflux, sigma=3, maxiters=5, masked=False)
        imrflux    = np.nanmedian(irfluxNew)
        istd       = np.nanstd(irfluxNew)

        # show the plots
        # plt.plot([-9.4, -5.2], [imrflux, imrflux], "r-", linewidth=1.5, label="$s=%.3f \pm %.3f$"%(imrflux,istd))
        # plt.plot([-9.4, -5.2], [imrflux+istd, imrflux+istd], "r--", linewidth=1.5)
        # plt.plot([-9.4, -5.2], [imrflux-istd, imrflux-istd], "r--", linewidth=1.5)
        # plt.plot(-2.5*np.log10(fluxRef[idRef]), irflux, "k.", ms=8)
        # plt.xlabel("AUTO MAG", fontsize=10)
        # plt.ylabel("FLUX SCALE", fontsize=10)
        # plt.legend(loc="upper left",fontsize=12)
        # plt.savefig(figdir + "photScale" + iimgKey.split("/")[-1][:-4]+"pdf")
        # plt.clf()
        # plt.close()

        fsList[iimgKey] = ifs * imrflux
    
    return fsList 


def lightcurve_fluxalign(savepath,listpath,tra,tdec,fluxs_dict):
    filelen = len(listpath)
    data_array= []
    for i in range(0,filelen):
        sciimg_name =listpath[i]
        flux_a = fluxs_dict[sciimg_name]
        mdict = mag_target(sciimg_name,tra,tdec)
        if (mdict ==0):
            continue
        print(mdict)
        data_array.append([mdict['fileid'],mdict['utctime'],mdict['mag_auto'],mdict['magerr_auto'],mdict['mag_psf'],mdict['magerr_psf'],mdict['mag_aper'],
                           mdict['magerr_aper'],mdict['amss'],mdict['jd'],mdict['tdb'],mdict['exptime'],mdict['ximage'],mdict['yimage'],mdict['ra'],mdict['dec'],mdict['iflag'],
                           mdict['fwhm'],mdict['snr'],mdict['kron_radius'],mdict['lbkg'],mdict['mag_auto']-(2.5*math.log10(flux_a)),mdict['mag_psf']-(2.5*math.log10(flux_a)),mdict['mag_aper']-(2.5*math.log10(flux_a)),flux_a])
                           
         
    txtheader = 'fileid,utctime,mag_auto,magerr_auto,mag_psf,magerr_psf,mag_aper,magerr_aper,airmass,JD,TDB,exptime,ximage,yimage,ra,dec,flag,fwhm,snr,k_radius,lbkg,mag_auto_a,mag_psf_a,mag_aper_a,flux_scale'
    
    t_light=Table(rows=data_array,names=['fileid','utctime','mag_auto','magerr_auto','mag_psf','magerr_psf','mag_aper','magerr_aper','airmass','JD','TDB','exptime','ximage','yimage','ra','dec','flag','fwhm','snr','k_radius','lbkg','mag_auto_a','mag_psf_a','mag_aper_a','flux_scale'])
     
    t_light.sort('utctime')
    
    #txtname=savepath+date+'_'+ttfname+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_l.txt'
    txtname=savepath+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_al.txt'
    tablename  = savepath+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_al.csv'
    np.savetxt(txtname, np.array(data_array), fmt="%s", delimiter=' ',header = txtheader)
    t_light.write(tablename,overwrite=True)
    return data_array,t_light
    
    

if __name__ == "__main__":
    # sciimg_name = '/home/mpilot/1m6/reception/20230527/sci/mb_ngc5457sn_v/20230527T191942/mb_sc_tngc5457sn_v_20230527191941_230_sciimg.fits'
    # #RA,DEC: 210.9106833, +54.3116718
    
    # dict = mag_target(sciimg_name,'210.9106833','+54.3116718')
    # print(dict)
    #dates=delta_date()
    dates= ['20230520','20230521','20230522','20230523','20230524','20230525','20230526','20230527','20230529','20230530','20230531','20230601','20230602','20230604','20230606','20230701']
    #dates = ['20230529']
    #sci_list = getfilelist(dates,'my_ngc5457sn_r')
    # for i in range(len(sci_list)):
    #     dict = mag_target(sci_list[i],210.9106833, +54.3116718)
    #     print(dict)
    # circle(14:03:38.5640,+54:18:42.020,6.000") # color=blue width=4 font="helvetica 8 bold roman" text={ngc5457_sn}
    # circle(14:03:45.1544,+54:16:15.341,8.666") # text={ref}
    # circle(14:04:02.8278,+54:18:18.226,9.839") # text={check}
    
    #y50a_H10HZ21_mi
    #y50b_H10HZ21_mz
    #y50a_H10HZ21_mr
    #y50b_H10HZ21_mg
    
    # 'y50b_NGC5457SN_mg'
    # 'y50b_NGC5457SN_mz'
    # 'y50a_NGC5457SN_mi'
    # 'y50a_NGC5457SN_mr'

    
    
    # savepath21 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_u_check'
    # savepath22 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_v_check'
    # savepath23 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_g_check'
    # savepath24 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_r_check'
    
    # savepath11 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_u_ref'
    # savepath12 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_v_ref'
    # savepath13 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_g_ref'
    # savepath14 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_r_ref'
    
    
    #savepath01 = '/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+'y50b_NGC5457SN_mg'
    # savepath02 = '/home/mpilot/1m6/reception/lightcurve/'+'mb_ngc5457sn_v'
    # savepath03 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_g'
    # savepath04 = '/home/mpilot/1m6/reception/lightcurve/'+'my_ngc5457sn_r'
    
    #mkdir('/home/mpilot/1m6/reception/lightcurve/')

    sci_list1 = getfilelist(dates,'y50b_NGC5457SN_mg')
    sci_list2 = getfilelist(dates,'y50b_NGC5457SN_mz')
    sci_list3 = getfilelist(dates,'y50a_NGC5457SN_mi')
    sci_list4 = getfilelist(dates, 'y50a_NGC5457SN_mr')
    
    tra0,tdec0=210.9106833, +54.3116718
    tra1,tdec1 = d2hms('14:03:45.1544', '+54:16:15.341', conv=1)
    tra2,tdec2 = d2hms('14:04:02.8278','+54:18:18.226', conv=1)

    savepath01 = '/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+'y50b_NGC5457SN_mg'
    savepath02 = '/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+'y50b_NGC5457SN_mz'
    savepath03 ='/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+ 'y50a_NGC5457SN_mi'
    savepath04 = '/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+'y50a_NGC5457SN_mr'

    #y50a，r波段：y50a_NGC5457SN_mr_N0034_20230521155137.fits_sciimg.fits
    sci_ref_g ='/home/50cm/50cm_scripts/50cm_pro/reception/20230521/sci/y50b_NGC5457SN_mg/20230521T155427/y50b_NGC5457SN_mg_N0034_20230521155427.fits_sciimg.fits'
    #sci_ref_r='/home/50cm/50cm_scripts/50cm_pro/reception/20230521/sci/y50a_NGC5457SN_mr/20230521T155137/y50a_NGC5457SN_mr_N0034_20230521155137.fits_sciimg.fits'
    sci_ref_r ='/home/50cm/50cm_scripts/50cm_pro/reception/20230521/sci/y50a_NGC5457SN_mr/20230522T183047/y50a_NGC5457SN_mr_N0071_20230522183047.fits_sciimg.fits'
    #sci_ref_i='/home/50cm/50cm_scripts/50cm_pro/reception/20230521/sci/y50a_NGC5457SN_mi/20230521T154235/y50a_NGC5457SN_mi_N0033_20230521154235_sciimg.fits'
    sci_ref_i='/home/50cm/50cm_scripts/50cm_pro/reception/20230521/sci/y50a_NGC5457SN_mi/20230521T154235/y50a_NGC5457SN_mi_N0033_20230521154235_sciimg.fits'
    sci_ref_z='/home/50cm/50cm_scripts/50cm_pro/reception/20230521/sci/y50b_NGC5457SN_mz/20230521T154656/y50b_NGC5457SN_mz_N0033_20230521154656.fits_sciimg.fits'
    #kdict = {"u": 0.7639, "v": 0.5411, "g": 0.2482, "r": 0.1846, "i": 0.0087, "z": 0.0056}
    kdict = {"u": 0.576, "v": 0.373, "mg": 0.254, "mr": 0.184, "mi": 0.074, "mz": 0.116}

    #kdict_g = flux_alignment(kdict,sci_list1,sci_ref_g)
    #kdict_z = flux_alignment(kdict,sci_list2,sci_ref_z)
    kdict_i = flux_alignment(kdict,sci_list3,sci_ref_i)
    #kdict_r = flux_alignment(kdict,sci_list4,sci_ref_r)

    #array01a,table01a =lightcurve_fluxalign(savepath01,sci_list1,tra0,tdec0,kdict_g)
    #array02a,table02a =lightcurve_fluxalign(savepath02,sci_list2,tra0,tdec0,kdict_z)
    array03a,table03a =lightcurve_fluxalign(savepath03,sci_list3,tra0,tdec0,kdict_i)
    #array04a,table04a =lightcurve_fluxalign(savepath04,sci_list4,tra0,tdec0,kdict_r)




    # array21,table21 =lightcurve_txt(savepath21,sci_list1,tra2,tdec2)
    # array22,table22 =lightcurve_txt(savepath22,sci_list2,tra2,tdec2)
    # array23,table23 =lightcurve_txt(savepath23,sci_list3,tra2,tdec2)
    # array24,table24 =lightcurve_txt(savepath24,sci_list4,tra2,tdec2)


    #print(tra1,tdec1)
    
    # array21,table21 =lightcurve_txt(savepath21,sci_list1,tra2,tdec2)
    # array22,table22 =lightcurve_txt(savepath22,sci_list2,tra2,tdec2)
    # array23,table23 =lightcurve_txt(savepath23,sci_list3,tra2,tdec2)
    # array24,table24 =lightcurve_txt(savepath24,sci_list4,tra2,tdec2)
    
    
    # array11,table11 =lightcurve_txt(savepath11,sci_list1,tra1,tdec1)
    # array12,table12 =lightcurve_txt(savepath12,sci_list2,tra1,tdec1)
    # array13,table13 =lightcurve_txt(savepath13,sci_list3,tra1,tdec1)
    # array14,table14 =lightcurve_txt(savepath14,sci_list4,tra1,tdec1)
    
    # array01,table01 =lightcurve_txt(savepath01,sci_list1,tra0,tdec0)
    # array02,table02 =lightcurve_txt(savepath02,sci_list2,tra0,tdec0)
    # array03,table03 =lightcurve_txt(savepath03,sci_list3,tra0,tdec0)
    # array04,table04 =lightcurve_txt(savepath04,sci_list4,tra0,tdec0)
    
    
    # sci_list1 = getfilelist(dates,'y50b_NGC5457SN_mg')
    # array01a,table01a =lightcurve_txt(savepath01,sci_list1,tra0,tdec0)



    # sci_ref_r = '/home/mpilot/1m6/reception/20230521/sci/my_ngc5457sn_r/20230521T153127/my_sc_tngc5457sn_r_20230521153125_106_sciimg.fits'
    # sci_ref_v='/home/mpilot/1m6/reception/20230521/sci/mb_ngc5457sn_v/20230521T153127/mb_sc_tngc5457sn_v_20230521153125_120_sciimg.fits'
    # sci_ref_u='/home/mpilot/1m6/reception/20230521/sci/mb_ngc5457sn_u/20230521T152836/mb_sc_tngc5457sn_u_20230521152835_117_sciimg.fits'
    # sci_ref_g = '/home/mpilot/1m6/reception/20230521/sci/my_ngc5457sn_g/20230521T152836/my_sc_tngc5457sn_g_20230521152835_103_sciimg.fits'
    # #kdict = {"u": 0.7639, "v": 0.5411, "g": 0.2482, "r": 0.1846, "i": 0.0087, "z": 0.0056}
    # kdict = {"u": 0.5943, "v": 0.4100, "g": 0.2482, "r": 0.1846, "i": 0.0087, "z": 0.0056}
    # kdict_u = flux_alignment(kdict,sci_list1,sci_ref_u)
    # kdict_v = flux_alignment(kdict,sci_list2,sci_ref_v)
    # kdict_g = flux_alignment(kdict,sci_list3,sci_ref_g)
    # kdict_r = flux_alignment(kdict,sci_list4,sci_ref_r)
     
    
    # array01a,table01a =lightcurve_fluxalign(savepath01,sci_list1,tra0,tdec0,kdict_u)
    # array02a,table02a =lightcurve_fluxalign(savepath02,sci_list2,tra0,tdec0,kdict_v)
    # array03a,table03a =lightcurve_fluxalign(savepath03,sci_list3,tra0,tdec0,kdict_g)
    # array04a,table04a =lightcurve_fluxalign(savepath04,sci_list4,tra0,tdec0,kdict_r)

#h10hz21: ra :183.4838990,dec:+32.9421783158

#y50a_H10HZ21_mi
    #y50b_H10HZ21_mz
    #y50a_H10HZ21_mr
    #y50b_H10HZ21_mg
    # h10_ra=183.4838990
    # h10_dec=32.9421783158
    
    # sci_list1 = getfilelist(dates,'y50a_H10HZ21_mi')
    # sci_list2 = getfilelist(dates,'y50b_H10HZ21_mz')
    # sci_list3 = getfilelist(dates,'y50a_H10HZ21_mr')
    # sci_list4 = getfilelist(dates,'y50b_H10HZ21_mg')
    
    # savepath01 = '/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+'y50a_H10HZ21_mi'
    # savepath02 = '/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+'y50b_H10HZ21_mz'
    # savepath03 ='/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+ 'y50a_H10HZ21_mr'
    # savepath04 = '/home/50cm/50cm_scripts/50cm_pro/reception/lightcurve/'+'y50b_H10HZ21_mg'

     
    # array01a,table01a =lightcurve_txt(savepath01,sci_list1,h10_ra,h10_dec)
    # array02a,table02a =lightcurve_txt(savepath02,sci_list2,h10_ra,h10_dec)
    # array03a,table03a =lightcurve_txt(savepath03,sci_list3,h10_ra,h10_dec)
    # array04a,table04a =lightcurve_txt(savepath04,sci_list4,h10_ra,h10_dec)
