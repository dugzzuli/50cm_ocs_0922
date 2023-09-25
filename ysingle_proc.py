# -*- coding:utf-8 -*-
import traceback
import os
import argparse
from lib.LogInstance import Logger
import matplotlib
from lib.phot.ybias import  get_current_dir
matplotlib.use('Agg')
from loguru import logger
from astropy.io import fits
import numpy as np
from execute.yphotutils import make_get_request,datename,mkdir
from execute import yphotutils
from execute import yprecom0_new as ypre
from execute import yrecpcom1_new as yr
from execute.yphost0 import single_cmd_deal
from yapercor import single_apercor_pro
from ymatch import match_pro

def executerYrecp_parallel(date):
    logger.info("starting Mrecp_parallel....")
    yr.pro(date)
    logger.info("ending Mrecp_parallel....")
def photometry(date,i_obj_file):
    logger.info('begin to photometry')
    rootpath=get_current_dir()
    
    
    logger.warning(f"i_obj_file:{i_obj_file}")
    try:
        tid,target,fid = yphotutils.ttfname_header(i_obj_file)
    except:
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": i_obj_file, "reason": f"preprocess is failed----{traceback.format_exc()}","stage": 0}
            response = make_get_request(url, params)
            return None,None    ###202309121658 zhangjh
        except:
            pass
    ######################avoid to obtain two "y50a" or y50b#################
    if target.startswith("y50"):
        ttfname=target+'_'+fid
    else:
        ttfname = tid+'_'+target+'_'+fid
    ttfname=ttfname.replace("__", "_")
    
    subpath = datename(i_obj_file) # obtain date_sub path
    
    if(subpath is None):
        logger.error('the file can not be open successfully')
        return None
    
    confdir = rootpath + "config/"
    ancdir  = rootpath + "ancData/"
    
    logdir = rootpath+'run_pipe_log/'+str(date)+'/'+ttfname+'/'
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir=rootpath+'reception/'+str(date)+'/sci/'+ttfname+'/'+subpath+'/'
    
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
    logger.info( tscidir+filename[:-5]+'_sciimg.fits')
    logger.info( tscidir+filename.split('.fit')[0]+'_sciimg.fits')
    
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
           
                

if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default="20221103")
    
    # parser.add_argument('--filename', type=str, default="/home/50cm/data2/raw/20230914_OCS/50A/y50a_p52251-Kepler_mi_N0036_20230914155818.fits")
    
    parser.add_argument('--filename', type=str, default="/home/50cm/data2/raw/20221103/y50b_ZTFJ2130+4420_-0012mg_.fit")
    args = parser.parse_args() 
    loggerloguru = Logger(args.date, '', "_ysingle").get_logger 
     
    photometry(args.date,args.filename)
    logger.info('the process is down')



