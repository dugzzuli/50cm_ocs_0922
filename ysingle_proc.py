import glob
import traceback
import os
import argparse
import gc
from multiprocessing import JoinableQueue, Process
import time
from config import config

from lib.LogInstance import Logger
import traceback
import matplotlib
import glob
from tqdm import tqdm
from lib.phot.ybias import  get_current_dir
matplotlib.use('Agg')
 
# from execute import yprecom1 as ypre
# from execute import yphostcom1 as ypho
# from execute import yrecpcom1 as yr
from execute import yprecom0_new as ypre
from execute import yphostcom0 as ypho
from execute import yrecpcom1_new as yr
 
# -*- coding: UTF-8 -*-
import argparse
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
import datetime
import os
import time
import traceback
from loguru import logger

from astropy.time import Time
from astropy import coordinates as coord, units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
# from PyAstronomy import pyas
import glob
import math
import time
import datetime as dt
import numpy as np
from execute.flatdiff import mkdir 
from  lib.phot import cal_airmass  
 

from execute.yphotutils import make_get_request,datename,ttfname_header,d2hms,mad,wds9reg,read_param,crossmatch,reg,ref_expt_select,readfits,hdr_edit_fits,re_hdr_edit_fits,read_list,timeplus,HJDTrans,get_scamp_head,overlapRect,pointRect,_spherical_to_cartesian,_great_circle_distance,mkdir

 
#from astropy.time.core import TimezoneInfo
from astropy.coordinates import EarthLocation
from astropy import units as u

from datetime import datetime,timezone,timedelta
from pytz import timezone
from tables.atom import Float32Atom

 
from execute import yprecom0_new as ypre
 
from execute import yrecpcom1_new as yr
from execute.yphost0 import single_cmd_deal
from yapercor import single_apercor_pro
from ymatch import match_pro
# from mphot_parallel import 
 
  

 
def executerYrecp_parallel(date):
    logger.info("starting Mrecp_parallel....")
    yr.pro(date)
    logger.info("ending Mrecp_parallel....")
     

 
def photometry(date,i_obj_file):
    logger.info('begin to photometry')
    rootpath=get_current_dir()
    
    #################20230911 112900  杜国王######################
    try:
        tid,target,fid = ttfname_header(i_obj_file)
    except:
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": i_obj_file, "reason": f"preprocess is failed----{traceback.format_exc()}","stage": 0}
            response = make_get_request(url, params)
            return None,None    ###202309121658 zhangjh



        except:
            pass
    #######################################
    
    ttfname = tid+'_'+target+'_'+fid
    
    subpath = datename(i_obj_file)
    if(subpath is None):
        logger.error('the file can not be open successfully')
        return None
    confdir = rootpath + "config/"
    ancdir  = rootpath + "ancData/"
    #figdir  = rootpath + "figures/"+str(date)+'/'+str(date)+'_'+ttfname+"/"
    logdir = rootpath+'run_pipe_log/'+str(date)+'/'+ttfname+'/'
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir=rootpath+'reception/'+str(date)+'/sci/'+ttfname+'/'+subpath+'/'
    #figsave_path = '/home/50cm/dugking/CrossPic/'+date+'/'+tid+'/'
    figdir= '/home/50cm/dugking/CrossPic/'+date+'/'+tid+'/'
     
    #scipath = upath+'sci/'+ttfname+'/'
    mkdir(upath)
    mkdir(tscidir)
    mkdir(logdir)
    mkdir(figdir)
    
    
    # 3S configuration
    sexParam1  = confdir + "default.param"
    sexParam2  = confdir + "default.paramNew"
    sexConf    = confdir + "default.sex"
    swarpConf  = confdir + "default.swarp"
    scampConf  = confdir + "default.scamp"
    scampConf2  = confdir + "second.scamp"
    psfexConf  = confdir + "default.psfex"

    sexComd1   = "sex %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s"
    sexComd2   = " -DETECT_MINAREA 3 -DETECT_THRESH 1.0 -ANALYSIS_THRESH 1.5"
    sexComd2_astro   = " -DETECT_MINAREA 5 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5"
    sexComd3   = " -SATUR_LEVEL %s -CHECKIMAGE_TYPE NONE"
    sexComd4   = " -PSF_NAME %s"
    sexComd    = sexComd1 + sexComd2 + sexComd3
    sexComdAstro =sexComd1 + sexComd2_astro + sexComd3
    sexComdPSF = sexComd + sexComd4
    scampComd1 = "scamp %s -c %s -ASTREFCAT_NAME %s -MOSAIC_TYPE UNCHANGED" # LOOSE "
    scampComd2 = " -FWHM_THRESHOLDS 2,20 -SN_THRESHOLDS 20,1000"
    scampComd  = scampComd1 + scampComd2
    
    psfexComd  = "psfex %s -c %s -PSF_SIZE %s,%s -PSF_DIR %s"
    
        
    libpath = rootpath + 'lib/cal3/'
    
    rawpath = rootpath + 'reception/' + str(date) + '/raw/'
    fileguide_raw_path = rawpath + 'fileguide/'
    calpath = rootpath + 'reception/' + str(date) + '/cal/'
    logger.info('the calpath is:' + calpath)
    scipath = rootpath + 'reception/' + str(date) + '/sci/' + ttfname + '/'
    mkdir(scipath)
    fileguide_sci_path = scipath + 'fileguide/'
    mkdir(fileguide_sci_path)
    logdir = rootpath + 'run_pipe_log/' + str(date) + '/' + ttfname + '/'
    mkdir(logdir)
    # if not (os.path.exists(fileguide_raw_path)):
    #     executerYrecp_parallel(date)
    # if not (os.path.exists(calpath + tid + '_master_bias.fits')):
    #     executerYrecp_parallel(date)
    #     ypre.calibration1(rawpath, calpath)
    # ycal.calibration_gen(rawpath, calpath)
    # dirname, filename = os.path.split(i_obj_file)
    if not (os.path.exists(fileguide_raw_path)):
        executerYrecp_parallel(date)
    if not (os.path.exists(calpath + tid + '_master_bias.fits')):
        ypre.calibration1(rawpath, calpath)
         
    if(os.path.exists(i_obj_file)):
        dirname, filename = os.path.split(i_obj_file)
        fileresult_list=ypre.combine_process(filename, dirname + '/', calpath, tscidir, libpath, date,repro=0)
        if not(fileresult_list is None):
            try:
                url = 'http://12.12.12.251:8888/process_request'
                params = {"eid": i_obj_file, "reason": "preprocess ok","stage": 0}
                response = make_get_request(url, params)
            except:
                pass
        else:
            try:
                url = 'http://12.12.12.251:8888/process_request'
                params = {"eid": i_obj_file, "reason": f"the preprocess is failed!  :{i_obj_file}","stage": 0}
                response = make_get_request(url, params)
            except:
                pass
            return None
 

    else:
        logger.info(f'there is no such file :{i_obj_file}')
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": i_obj_file, "reason": f"there is no such file  :{i_obj_file}","stage": 0}
            response = make_get_request(url, params)
        except:
            pass
        return None
            
            # return
    logger.info( tscidir+filename[:-5]+'_sciimg.fits')
    logger.info( tscidir+filename.split('.fit')[0]+'_sciimg.fits')
    #i_obj_sci =tscidir+filename[:-5]+'_sciimg.fits'
    #single_pho=single_cmd_deal(i_obj_sci,i_obj_file,ttfname,date,sexComd,sexComdAstro,sexConf,sexParam1,psfexComd,psfexConf,scampComd, scampConf, scampConf2,sexComdPSF,sexParam2,rootpath,figdir,ancdir,tscidir,confdir)
    filename1=filename.split('.fit')[0]
    if(filename1[-1]=='_'):
        filename1 = filename1[:-1]
    if os.path.exists( tscidir+filename1+'_sciimg.fits'):
         
        i_obj_sci =tscidir+filename1+'_sciimg.fits'
        single_pho=single_cmd_deal(i_obj_sci,i_obj_file,ttfname,date,sexComd,sexComdAstro,sexConf,sexParam1,psfexComd,psfexConf,scampComd, scampConf, scampConf2,sexComdPSF,sexParam2,rootpath,figdir,ancdir,tscidir,confdir)
        i_obj_sexcat =tscidir+filename1+'_sciimg_sexcat.fits'
        logger.info('++++++++++++++++++++++++++single_pho=')
        logger.info(single_pho)
        logger.info(single_pho[0])
        if not single_pho[0] is None:
            logger.info('++++++++++++++++++++++++++single_pho is none')
            logger.info(single_pho[0])
            if (os.path.exists(i_obj_sexcat)):
                catdata=fits.getdata(i_obj_sexcat,ext=2)
                try:
                    #logger.info(np.median(catdata['MAG_PSF']))
                    if not (np.median(catdata['MAG_PSF'])==0.0):
                        try:
                            url = 'http://12.12.12.251:8888/process_request'
                            params = {"eid": i_obj_file, "reason": "photometry is ok","stage": 2}
                            response = make_get_request(url, params)
                        except:
                            pass
            
                except Exception as e:
                    
                    logger.error(str(traceback.format_exc()))
                    
                    try:
                        url = 'http://12.12.12.251:8888/process_request'
                        params = {"eid": i_obj_file, "reason": f"{str(traceback.format_exc())}","stage": 2}
                        response = make_get_request(url, params)
                    except:
                        pass
    if os.path.exists( tscidir+filename1+'_sciimg_sexcat.fits'):
        i_obj_sexcat = filename1+'_sciimg_sexcat.fits'

    
        if not single_pho[0] is None:
            single_apercor_pro(date, ttfname, tscidir,figdir, i_obj_sexcat)
            i_obj_sexcat_aper = tscidir+filename1+'_sciimg_sexcat_apcorr.fits'
            logger.info(i_obj_sexcat_aper)
            if (os.path.exists(i_obj_sexcat_aper )):
                try:
                    url = 'http://12.12.12.251:8888/process_request'
                    params = {"eid": i_obj_file, "reason": "aperture correction is ok","stage": 3}
                    response = make_get_request(url, params)
                except:
                    pass
    
    if os.path.exists( tscidir+filename1+'_sciimg_sexcat.fits'):
        i_obj_sexcat =filename1+'_sciimg_sexcat.fits'
        refpath='/home/50cm/data2/'
        if not single_pho[0] is None:
            match_pro(date,tscidir,refpath,figdir)
            i_obj_sexcat_gaia = tscidir+filename1+'_sciimg_sexcat_gaia.fits'
        
            if (os.path.exists(i_obj_sexcat_gaia)):
                try:
                    url = 'http://12.12.12.251:8888/process_request'
                    params = {"eid": i_obj_file, "reason": "data process is ok","stage": 4}
                    response = make_get_request(url, params)
                except:
                    pass
        #executerMapercor_parallel(date, ttfname, tscidir,figdir, i_obj_sexcat)
    
                 
                

if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default="20221029", help='输入处理日期')
    #parser.add_argument('--filename', type=str, default="/home/50cm/data2/raw/20230316_OCS/50B/y50b_TIC20212631_mg_N0107_20230316165050.fits", help='输入处理日期')
    #parser.add_argument('--date', type=str, default="20221102", help='输入处理日期')
    # parser.add_argument('--filename', type=str, default="/home/50cm/data2/raw/20230302/y50b_NGC5915_mg_0018_20230302213927.fit", help='输入处理日期')
    parser.add_argument('--filename', type=str, default="/home/50cm/data2/raw/20221029/y50b_H02 GD50_-0003mr_.fit", help='输入处理日期')
    
     
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_ysingle").get_logger 
     
    photometry(args.date,args.filename)
    logger.info('the process is down')



