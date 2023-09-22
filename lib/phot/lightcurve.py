import os
import sys
from loguru import logger
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy import coordinates as coord, units as u
import datetime
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
# from PyAstronomy import pyas
import glob
import lib.phot.cal_airmass as cal_airmass
import math
from loguru import logger


def ams_cal(filename,ira,idec):
    res=fits.open(filename)
    data = res[1].data

    date_obs = data.field(0)[0,5][11:30]

    expt     = float(data.field(0)[0,6][12:].split(' ')[0])
    amss,jd_mid,tbd_mid = cal_airmass.cam(date_obs,expt,ira,idec)  
    return amss
    
from lib.phot.ybias import get_current_dir

def plotshow(x,y,ye,ylabel,title,filename):
    plt.switch_backend('agg')
    plt.errorbar(x, y, yerr=ye,fmt='o',ecolor='r',color='b',elinewidth=2,capsize=4)
    #plt.xlim(0, 0.7)
    plt.xlabel('MJD-2459334.3(days)')
    plt.ylabel(ylabel)
    plt.title(title)
    #plt.clf()
    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(filename)

def getindex_location(tra,tdec,ralist,declist):
    for i in range(len(ralist)):
        d=np.sqrt((ralist[i]-tra)**2+(declist[i]-tdec)**2)
        #logger.info(ra[i],',',dec[i],',d=',d)
        if(d<5/3600):
            return i,ralist[i],declist[i]

def TransLDAC_time(ldacfilename):
    data=fits.getdata(ldacfilename,ext=1)
    obstime_0=np.array(data[0][0])[5]
    exptime=np.array(data[0][0])[7]
    obstime=obstime_0.split('/')[0].split('=')[1].replace('\'',"").strip().replace('T',' ')#.strip('\'')
    exptime=exptime.split('/')[0].split('=')[1]
    obstime=datetime.datetime.strptime(obstime,'%Y-%m-%d %H:%M:%S')
    # 计算偏移量
    offset = datetime.timedelta(seconds=float(exptime)/2)
    # 获取修改后的时间并格式化
    timeplus = obstime + offset
    timeplus=str(timeplus).replace(' ','T')
    t=Time(timeplus)
    tjd=Time(np.array(str(t.jd)),format='jd')
    tb=tjd.tdb
    return str(obstime),str(timeplus),str(tb)

def exposur(ldacfilename):
    data=fits.getdata(ldacfilename,ext=1)
    exptime=np.array(data[0][0])[7]
    exptime=exptime.split('/')[0].split('=')[1]
    return float(exptime)

def reg(rootpath,target,date,tra,tdec):
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'
    ttf_g='y50a_'+target+'_mg'
    ttf_r='y50b_'+target+'_mr'
    tscidir_g = upath+'sci/'+ttf_g+'/'
    tscidir_r = upath+'sci/'+ttf_r+'/'
    
    flist_g = glob.glob(tscidir_g +'*_sexcat.fits')
    flist_r = glob.glob(tscidir_r +'*_sexcat.fits')
    flist_g =np.sort(flist_g)
    flist_r =np.sort(flist_r)
    
    data_g=[]
    for i in range(len(flist_g)):
        ildac = fits.getdata(flist_g[i],ext=2)
        logger.info(flist_g[i][:-21]+'.fits') 
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        #mag_auto, magerr_auto=ildac["MAG_AUTO"],ildac["MAGERR_AUTO"]
        #mag_aper,magerr_aper=ildac["MAG_APER"],ildac["MAGERR_APER"]
        #mag_psf,magerr_psf=ildac["MAG_PSF"],ildac["MAGERR_PSF"]
        img_x,imgy=ildac["X_IMAGE"],ildac["Y_IMAGE"]
        index,ira,idec=getindex_location(tra,tdec,ra, dec)
        data_g.append([img_x[index],imgy[index]])
        

    regListn   = rootpath + ttf_g+".reg"
    regList    = open(regListn,"w")
    regList.write('# Region file for DS9'+"\n")
    regList.write('global color=green font=\'helvetica 10 normal\' select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'+"\n")
    regList.write('image'+"\n")
    for i in range(0,len(data_g)):
            sline='circle('+str(data_g[i][0])+','+ str(data_g[i][1])+', 8.00)' +'# color=cyan font="helvetica 8 normal" text={' + 'id='+str(i+1)+'}'
            regList.write(sline+"\n")
    regList.close()
    logger.info('the reg_g is ok')
    

    data_r=[]
    for i in range(len(flist_r)):
        ildac = fits.getdata(flist_r[i],ext=2)
        logger.info(flist_r[i][:-21]+'.fits') 
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        img_x,imgy=ildac["X_IMAGE"],ildac["Y_IMAGE"]
        index,ira,idec=getindex_location(tra,tdec,ra, dec)
        data_r.append([img_x[index],imgy[index]])
        

    regListn   = rootpath + ttf_r+".reg"
    regList    = open(regListn,"w")
    regList.write('# Region file for DS9'+"\n")
    regList.write('global color=green font=\'helvetica 10 normal\' select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'+"\n")
    regList.write('image'+"\n")
    for i in range(0,len(data_r)):
            sline='circle('+str(data_r[i][0])+','+ str(data_r[i][1])+', 8.00)' +'# color=cyan font="helvetica 8 normal" text={' + 'id='+str(i+1)+'}'
            regList.write(sline+"\n")
    regList.close()
    logger.info('the reg_r is ok')
def lc_txt_aper(rootpath,date,ttfname,tra,tdec):
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'+ttfname+'/'
    falist = glob.glob(tscidir +'*sexcat.fits')
    #fblist = glob.glob(tscidir +'*sexcat_apcorr.fits')
    #falist=np.sort(falist)
    ldacList=np.sort(falist)
    #fblist=np.sort(fblist)
    # if (len(falist)>len(fblist)):
    #     ldacList=np.sort(fblist)
    #     logger.info('the aperture corrected sexcat files may not be complete! It may lead to some records are lack of Aper_cor Magnitude!')
    # else:
    #     ldacList=np.sort(fblist)
 
    filtert=ttfname[-1]
    #logger.info(ldacList)
    logger.info(ldacList)
    data_array=[]
    logger.info(str(len(ldacList))+'################################## mags of target in this file!')
    for i in range(len(ldacList)):
        ildac1name=ldacList[i]
        aper_cor_ildac=ldacList[i][:-5]+'_apcorr.fits'
        apermag =0.0
        apererr =0.0
        #apermag=np.nan
        if(os.path.exists(aper_cor_ildac)):
            ildac = fits.getdata(aper_cor_ildac,ext=1)
            mag_aper,magerr_aper=ildac["MAG_FAPER"],ildac["MAGERR_FAPER"]
        else: 
            ildac=fits.getdata(ldacList[i],ext=1)
            mag_aper=[]
            logger.info('there is no aper_cor mag')

        exp=exposur(ildac1name)
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        mag_auto, magerr_auto=ildac["MAG_AUTO"],ildac["MAGERR_AUTO"]
        #mag_aper,magerr_aper=ildac["MAG_FAPER"],ildac["MAGERR_FAPER"]
        mag_psf,magerr_psf=ildac["MAG_PSF"],ildac["MAGERR_PSF"]
        
        index,ira,idec=getindex_location(tra,tdec,ra, dec)
         
        ams=ams_cal(ildac1name,ira,idec)
        psfmag=mag_psf[index]
        #apermag=mag_aper[index]#+math.log10(exp)*2.5
        psfmag=mag_psf[index]#+math.log10(exp)*2.5#+ztmag
        if len(mag_aper)>0:
            apermag=mag_aper[index]+math.log10(exp)*2.5
            apererr=magerr_aper[index]
        #logger.info('index',index)
        #logger.info('psfmag',mag_psf[index])
        psferr=magerr_psf[index]
        ximage,yimage= ildac["X_IMAGE"][index],ildac["Y_IMAGE"][index]

  
        fileid=str(aper_cor_ildac).split('-')[1][:4]
        obstime,utctime,tdbtime=TransLDAC_time(ildac1name)
        data_array.append([utctime,tdbtime,ximage,yimage,ira,idec,ams,psfmag,psferr,apermag,apererr,ximage,yimage])
        logger.info('fileid,utctime,ira,idec,ximage,yimage,ams,psfmag,aper_mag=',fileid,utctime,ira,idec,ximage,yimage,ams,psfmag,apermag)


    txtname=rootpath+date+'_'+ttfname+'_'+str(round(tra,4))+'+'+str(round(tdec,4))+'_l.txt'
    np.savetxt(txtname, np.array(data_array), fmt="%s", delimiter=' ')
def pro(ttfname,date,tra,tdec):
    #tra=322.7362855
    #tdec=44.346236
    #J043610.27+195226.1
    #tra,tdec=69.0428297,19.87394052
    #g193-74
    #tra,tdec=118.3636,52.4920
    #hz21
    #tra,tdec=183.48443,32.942044
    #g191-b2b
    #tra,tdec=76.3775,52.830
   #refcheck #
    #tra,tdec=118.41283983333334,52.51145416666667
    #refcheck#tra,tdec= 76.43471333333333 ,52.82121611111111
    
    rootpath=get_current_dir()
    lc_txt_aper(rootpath,date,ttfname,tra,tdec)
    #reg(rootpath,ttfname,date,tra,tdec)
    logger.info('the lightcurve process is done')
    

# 国王注释
# ttfname = sys.argv[1]
# date = sys.argv[2]
# tra=sys.argv[3]
# tdec=sys.argv[4]
# tra=float(tra)
# tdec=float(tdec)
# pro(ttfname,date,tra,tdec)
#python lightcurve.py y50a_HZ21_mg 20220104 183.48443 32.942044
