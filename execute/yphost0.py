import datetime
import gc
import logging
import math
import os
import sys
import time
import traceback
from lib.phot.ybias import get_current_dir
import ccdproc
import numpy as np
#import matplotlib
#import pylab as pl
from astropy import coordinates as coord
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clip,sigma_clipped_stats
from astropy.table import Table,hstack
from astropy.time import Time
from loguru import logger
from scipy.optimize import leastsq
from scipy.spatial import cKDTree as ckdt
import ntpath
import lib.phot.PlotFun as PlotFun
 
from lib.LogInstance import Logger

from lib.phot.yrefsed import refSed_gaia,refSed_gaia_filter
import subprocess
import config

from execute.yphotutils import make_get_request,datename,d2hms,mad,wds9reg,read_param,crossmatch,reg,ref_expt_select,readfits,hdr_edit_fits,re_hdr_edit_fits,read_list,timeplus,HJDTrans,get_scamp_head,overlapRect,pointRect,_spherical_to_cartesian,_great_circle_distance,mkdir


def astronet_comd(astrocomd):
     
    astroComdlist=[]
    astroComdlist=astrocomd.split(" ")
    #try:
    if(1==1):
        # 执行命令并设置超时时间为300秒
        completed_process_astrocomd = subprocess.run(astroComdlist, timeout=300, shell=False, capture_output=True)
        #completed_process_astrocomd = subprocess.call(astroComdlist, timeout=300, shell=False, capture_output=True)
        #logger.info(iscampComdlist)           
        # 处理命令执行结果
        if completed_process_astrocomd.returncode == 0:
           
            logger.info(completed_process_astrocomd.stdout)
        else:
            logger.info(completed_process_astrocomd.stderr)
            logger.info('XXXXXXXXXXXXXXXastro failedXXXXXXXXXXXXXXXX')
    
    # except subprocess.TimeoutExpired:
    #     logger.info("Command scamp_comd execution timed out.")
    #     return None, None


def iwcs(iobjFrame_o0,config_cfg):
    wcsfile=iobjFrame_o0[:-5]+'.wcs'
    if not os.path.exists(wcsfile):
        wcs(iobjFrame_o0,config_cfg)
    if os.path.exists(iobjFrame_o0[:-5]+'.head'):
        ihead = iobjFrame_o0[:-5]+'.head'
        try:
            os.system("rm %s" % ihead)
        except:
            logger.info(ihead +' delete failed')
    if os.path.exists(iobjFrame_o0[:-5]+'.ahead'):
        ahead = iobjFrame_o0[:-5]+'.ahead'
        try:
            os.system("rm %s" % ahead)
        except:
            logger.info(ahead +' delete failed')
    iimgMat,hdr_sciimg = fits.getdata(iobjFrame_o0,header=True)
    if "PV1_1" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_1"]
    if "PV1_2" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_2"]
    if "PV1_0" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_0"]
    if "PV1_4" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_4"]
    if "PV1_5" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_5"]
    if "PV1_6" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_6"]
    if "PV1_7" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_7"]
    if "PV1_8" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_8"]
    if "PV1_9" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_9"]
    if "PV1_10" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_10"]
    if "PV2_1" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_1"]
    if "PV2_2" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_2"]
    if "PV2_0" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_0"]
    if "PV2_4" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_4"]
    if "PV2_5" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_5"]
    if "PV2_6" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_6"]
    if "PV2_7" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_7"]
    if "PV2_8" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_8"]
    if "PV2_9" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_9"]
    if "PV2_10" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_10"]
        
        
    if(os.path.exists(wcsfile)):
        try:
            logger.info(wcsfile)
            
            hdr_wcs = fits.getheader(wcsfile)
            logger.info(hdr_wcs["CD1_1"])
            logger.info(hdr_wcs["CD1_2"])
            logger.info(hdr_wcs["CD2_1"])
            logger.info(hdr_wcs["CD2_2"])

            hdr_sciimg["CTYPE1"] = "RA---TPV"
            hdr_sciimg["CTYPE2"] = "DEC--TPV"
            hdr_sciimg["RADESYS"] = "ICRS"
            hdr_sciimg["CUNIT1"] = "deg"
            hdr_sciimg["CUNIT2"] = "deg"
            
            hdr_sciimg["CRVAL1"]  = hdr_wcs["CRVAL1"]
            hdr_sciimg["CRVAL2"]  = hdr_wcs["CRVAL2"]
            
            hdr_sciimg["CRPIX1"]  = hdr_wcs["CRPIX1"]
            hdr_sciimg["CRPIX2"]  = hdr_wcs["CRPIX2"]
            

            hdr_sciimg["CD1_1"]   = hdr_wcs["CD1_1"]
            hdr_sciimg["CD1_2"]   = hdr_wcs["CD1_2"]
            hdr_sciimg["CD2_1"]   = hdr_wcs["CD2_1"]
            hdr_sciimg["CD2_2"]   = hdr_wcs["CD2_2"]
            
            fits.writeto(iobjFrame_o0, iimgMat, hdr_sciimg, overwrite=True)

        except:
            logger.info('these are no wcs of in header of this fits')
            logger.error(str(traceback.format_exc()))
 

def wcs(iobjFrame_o0,config_cfg):
    ihdr = fits.getheader(iobjFrame_o0)
    if os.path.exists(iobjFrame_o0[:-5]+'.head'):
        ihead = iobjFrame_o0[:-5]+'.head'
        try:
            os.system("rm %s" % ihead)
        except:
            logger.info(ihead +' delete failed')
    if os.path.exists(iobjFrame_o0[:-5]+'.ahead'):
        ahead = iobjFrame_o0[:-5]+'.ahead'
        try:
            os.system("rm %s" % ahead)
        except:
            logger.info(ahead +' delete failed')
    
    if("OBJCTRA" in ihdr.keys()):
        ira  = ":".join(ihdr["OBJCTRA"].split())
        idec = ":".join(ihdr["OBJCTDEC"].split())
        tra, tdec = d2hms(ira, idec, conv=1)
    elif("RA" in ihdr.keys()):
        ira  = ":".join(ihdr["RA"].split())
        idec = ":".join(ihdr["DEC"].split())
        tra, tdec = d2hms(ira, idec, conv=1)


    tra= str(tra) 
    tdec=str(tdec)  
    #astroComd = "time solve-field %s --scale-low 0.4 --scale-high 0.7  --scale-units app --config %s  --use-sextractor --no-plots --overwrite"    
    # isastroComd = astroComd%(iobjFrame_o0,config_cfg)#,str(tra),str(tdec))
    
    ##
    binfact=int(ihdr['XBINNING'])
    pscale  = 0.3 * binfact # pixel scale in unit of arcsec
    
    scale_low= pscale-0.1
    scale_high=pscale+0.1
    astroComd = "time solve-field %s --scale-low %s --scale-high %s --scale-units app --config %s --ra %s --dec %s --radius 3 --use-sextractor --no-plots --overwrite"    
    isastroComd = astroComd%(iobjFrame_o0,str(scale_low),str(scale_high),config_cfg,tra,tdec)
     
    ##
    logger.info(isastroComd)
    astronet_comd(isastroComd)
    
     
    try:
        os.system("rm %s" % iobjFrame_o0[:-5]+'.match')
        os.system("rm %s" % iobjFrame_o0[:-5]+'.new')
        os.system("rm %s" % iobjFrame_o0[:-5]+'.rdls')
        os.system("rm %s" % iobjFrame_o0[:-5]+'.axy')
        os.system("rm %s" % iobjFrame_o0[:-5]+'.corr')
        os.system("rm %s" % iobjFrame_o0[:-5]+'-indx.xyls')
        os.system("rm %s" % iobjFrame_o0[:-5]+'.solved')
            
    except:
        logger.error('there is no any pdf file')
        
    wcsfile=iobjFrame_o0[:-5]+'.wcs'
    iimgMat,hdr_sciimg = fits.getdata(iobjFrame_o0,header=True)
    #logger.info(hdr_sciimg)
    if(os.path.exists(wcsfile)):
        try:
            
            logger.info(wcsfile)
            
            hdr_wcs = fits.getheader(wcsfile)
            logger.info(hdr_wcs["CD1_1"])
            logger.info(hdr_wcs["CD1_2"])
            logger.info(hdr_wcs["CD2_1"])
            logger.info(hdr_wcs["CD2_2"])
            logger.info(f'sci_CRPIX1='+str(hdr_sciimg["CRPIX1"]))
            logger.info(f'sci_CRPIX2='+str(hdr_sciimg["CRPIX2"]))
            logger.info(f'wcs_CRPIX1='+str(hdr_wcs["CRPIX1"]))
            logger.info(f'wcs_CRPIX2='+str(hdr_wcs["CRPIX2"]))
            hdr_sciimg['EQUINOX']=2000.0
            
            
            
            hdr_sciimg["CUNIT1"]  = hdr_wcs["CUNIT1"]
            hdr_sciimg["CUNIT2"]  = hdr_wcs["CUNIT2"]
            
            hdr_sciimg["CRVAL1"]  = hdr_wcs["CRVAL1"]
            hdr_sciimg["CRVAL2"]  = hdr_wcs["CRVAL2"]
            
            hdr_sciimg["CRPIX1"]  = hdr_wcs["CRPIX1"]
            hdr_sciimg["CRPIX2"]  = hdr_wcs["CRPIX2"]
            
            hdr_sciimg["CTYPE1"] = "RA---TPV"
            hdr_sciimg["CTYPE2"] = "DEC--TPV"
            hdr_sciimg["RADESYS"] = "ICRS"
             
            

            hdr_sciimg["CD1_1"]   = hdr_wcs["CD1_1"]
            hdr_sciimg["CD1_2"]   = hdr_wcs["CD1_2"]
            hdr_sciimg["CD2_1"]   = hdr_wcs["CD2_1"]
            hdr_sciimg["CD2_2"]   = hdr_wcs["CD2_2"]
            
            fits.writeto(iobjFrame_o0, iimgMat, hdr_sciimg, overwrite=True)

        except:
            logger.info('these are no wcs of in header of this fits')
            logger.error(str(traceback.format_exc()))


  
def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)

 

def single_cmd_deal(i_obj_file,source_file,ttfname,date,sexComd,sexComd_Astro,sexConf,sexParam1,psfexComd,psfexConf,scampComd, scampConf, scampConf2,sexComdPSF,sexParam2,rootpath,figdir,ancdir,tscidir,configdir):
    loggerloguru = Logger(date, '', "_yphost").get_logger
    logger.info('####### the processed fits is: %s' % (i_obj_file))
    objfilename = ntpath.basename(i_obj_file)
    objfilename.split('_sciimg')[0]+'.fits'

    ildacn,iobjFrame_o0=None,None
    
    #202309081631 增加try catch 杜国王  目的捕获其中某个异常，先能查看日志
    try:
    # if(1==1):
        # 对每一个frame run sextractor
         
        iobjFrameNew = i_obj_file #.split('sciimg')[0] + 'fcimg.fits'
         
        iobjFrame_o0 = i_obj_file
        hdr = fits.getheader(iobjFrameNew)
        
        ########202309081631 增加radec判断 杜国王##########
        
        if("OBJCTRA" in hdr.keys()):
            pass
        elif("RA" in hdr.keys()):
            pass
        else:
            
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": source_file, "reason": "RA or DEC doesn't exist","stage": 0}
            response = make_get_request(url, params)
            
            raise("XXXXXXXXXXXXXXXXX the fits file is broken: %s"%(iobjFrameNew))
            
            
            
    
        ##################
        logger.info(iobjFrameNew)
        logger.info(iobjFrame_o0)
        ################################################
        ipsfFrameNew = iobjFrame_o0[:-4] + "psf"
         
        
        
        
        ihdr,iimgMat= readfits(iobjFrame_o0)
        r,c=np.shape(iimgMat)

 
        binfact=int(ihdr['XBINNING'])
        pscale  = 0.297 * binfact # pixel scale in unit of arcsec
        ira='0'
        idec='0'
        try:
            if "OBJCTRA" in ihdr.keys():
                ira  = ":".join(ihdr["OBJCTRA"].split())
                idec = ":".join(ihdr["OBJCTDEC"].split())
                ira, idec = d2hms(ira, idec, conv=1)
                ira =str(ira)
                idec=str(idec)
            elif "RA" in ihdr.keys(): 
                #logger.info(ihdr["RA"],ihdr["DEC"])
                ira  = ":".join(ihdr["RA"].split())
                idec = ":".join(ihdr["DEC"].split())
                # logger.info('split --- idec')
                # logger.info(idec)
                ira =str(ira)
                idec=str(idec)
                ira, idec = d2hms(ira, idec, conv=1)
            logger.info(f'ira={ira}, idec={idec}')

 
        except:
            logger.error(str(traceback.format_exc()))
            
            logger.info('these are no parameters of object ra and dec in header of this fits')

            try:
                url = 'http://12.12.12.251:8888/process_request'
                params = {"eid": source_file, "reason": "these are no parameters of object ra and dec in header of this fits: -------"+str(traceback.format_exc()),"stage": 2}
                response = make_get_request(url, params)
                return None, None
            except: 
                pass
            return None, None
        
             
        
 
    
        #ira,idec = hdr_edit_fits(i_obj_file)
        
        if(c>8000 and binfact==1) or (c> 4000 and binfact ==2) :

            iimgMat = iimgMat[:,int(1000/binfact):c-int(1000/binfact)]
            r,c=np.shape(iimgMat)
            logger.info(f"#################imgcut################## :{np.shape(iimgMat)}")
           
        fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True)
        ira,idec =hdr_edit_fits(iobjFrameNew)

        
        #######astrometry##############
        
        iimgMat,ihdr=fits.getdata(iobjFrame_o0,header=True)
        # logger.info(ihdr)          ###202309121737  zhangjh
        exptime=ihdr['EXPTIME']
        ildacn = iobjFrame_o0[:-4] + "ldac"
        
        saturate_value=ihdr['SATURATE']*0.9
        exptlast = ihdr['EXPTIME']
        logger.error(ttfname)
        
        objID = hdr['OBJECT'].strip().replace(" ", "")
        if str(objID).__contains__("y50a") or str(objID).__contains__("y50b"):
            objID=objID[5:]
        filterid = 'm'+hdr['FILTER'].strip()
        ############################################
    
        refCatn  = ancdir + "GaiaStar_%s_%s_%s.ldac"%(objID,filterid,date)
        refCatn=refCatn.replace("__", "_")
        refCat2n = ancdir + "GaiaStar_%s_%s_%s.cat"%(objID,filterid,date)
        refCat2n=refCat2n.replace("__", "_")
        refCat2=None
    
        logger.info(refCat2n)
        if not (os.path.exists(refCat2n)):#(1==1):
            logger.info(f'ira={ira}, idec={idec}')
            maxmaglim=ref_expt_select(configdir,filterid,binfact,exptlast)
            refSed_gaia_filter(rootpath,ttfname,date,ira,idec,10,maxmaglim)
            #refSed_gaia(rootpath,ttfname,date,ira,idec,10,maxmaglim)
        else:
            try:
                refCat2  = Table.read(refCat2n, format="ascii")
            except:
                maxmaglim=ref_expt_select(configdir,filterid,binfact,exptlast)
                refSed_gaia_filter(rootpath,ttfname,date,ira,idec,10,maxmaglim)
            finally:
                
                logger.info(f'ira={ira}, idec={idec}')
                
                refCat2  = Table.read(refCat2n, format="ascii")
                
        if not refCat2:
            raise(f'''Table.read({refCat2n}, format="ascii") is failed''')
        
        logger.info(refCat2n)
        raRef    = refCat2["ra"]
        decRef   = refCat2["dec"]
        magRef  = refCat2["mag"]
        
        
        wcsfile=source_file.split('.fit')[0]+'.wcs'
        logger.info(wcsfile)
        if not (os.path.exists(wcsfile)):
            config_cfg = configdir + 'fy_astro.cfg'
                 
            wcs(iobjFrame_o0,config_cfg)
            
        try:
            
            isexComd = sexComd_Astro % (iobjFrame_o0, sexConf, sexParam1, ildacn, str(saturate_value))
            logger.info(isexComd)
            
            s_return=sex_comd(isexComd)
            if(s_return==0):
                try:
                    url = 'http://12.12.12.251:8888/process_request'
                    params = {"eid": source_file, "reason": f"source detection by sextractor is failed","stage": 1}
                    response = make_get_request(url, params)
                    return None,None
                except:
                    pass
                logger.info(f"source detection by sextractor is failed")
                return None,None
                
            
            ildac=fits.getdata(ildacn,ext=2)
         
            try:
                ixcenloc=np.where((ildac["X_IMAGE"]>c/2+10) | (ildac["X_IMAGE"]<c/2-10))
                ildac =ildac[ixcenloc]
                hdul=fits.open(ildacn)###
                hdul[1].data=ildac
                hdul.close()
            except: 
                logger.info('XXXXXXXXXXXXXXXXXXXX  ldac cut failed ')
            logger.info('astrometry begin')
            astrometry(ildacn,refCatn,iobjFrame_o0,scampComd, scampConf, scampConf2 ,figdir)

            ildac = fits.getdata(ildacn, ext=2)
            ira = ildac["ALPHA_J2000"]
            idec = ildac["DELTA_J2000"] 
            idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)     
            logger.info(f"crossmatch(raRef, decRef, ira, idec, aperture=1.5),len(idRef)={len(idRef)},len(ids)={len(ids)}")
            #he=fits.getheader(iobjFrame_o0)
            #logger.info(he)
            if len(idRef)<20:
                logger.info(f'len(idRef)={len(idRef)}')
                re_hdr_edit_fits(iobjFrame_o0)
                config_cfg = configdir + 'fy_astro.cfg'
                iwcs(iobjFrame_o0,config_cfg)
                isexComd = sexComd_Astro % (iobjFrame_o0, sexConf, sexParam1, ildacn, str(saturate_value))
                sex_comd(isexComd)
                try:
                    ixcenloc=np.where((ildac["X_IMAGE"]>c/2+10) | (ildac["X_IMAGE"]<c/2-10))
                    #ixyloc=np.where((iximg<c-(1000/binfact))&(iximg>1000/binfact))#&(iyimg<6100)&(iyimg>30))

                    ildac =ildac[ixcenloc]
                    hdul=fits.open(ildacn)###
                    hdul[1].data=ildac
                    hdul.close()
                except: 
                    logger.info('XXXXXXXXXXXXXXXXXXXX  ldac cut failed ')
                astro=astrometry(ildacn,refCatn,iobjFrame_o0,scampComd, scampConf, scampConf2,configdir)

                if astro is None:
                    #logger.error(f"os.path.exists(ihead)==False: {os.path.exists(ihead)==False}")
                    try:
                        url = 'http://12.12.12.251:8888/process_request'
                        params = {"eid": source_file, "reason": f"astrometry is failed,crossmatch len(idRef)={len(idRef)},len(ids)={len(ids)}","stage": 1}
                        #params = {"eid": source_file, "reason": "astrometry is failed","stage": 1}
                        response = make_get_request(url, params)
                        return None,None
                    except:
                        pass
                    logger.info(f"the astrometry is failed, ending this loop")
                    return None,None
                
            ildac = fits.getdata(ildacn, ext=2)
            ira = ildac["ALPHA_J2000"]
            idec = ildac["DELTA_J2000"] 
             
            idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)     
            
            if len(idRef)<20 :
                    
                logger.info(f'len(idRef)={len(idRef)}')
                logger.info(f"astrometry is failed!{iobjFrame_o0}")
                logger.error(f"astrometry is failed! crossmatch(raRef, decRef, ira, idec, aperture=1.5),len(idRef)={len(idRef)},len(ids)={len(ids)}")
                logger.info(f"astrometry is failed! crossmatch(raRef, decRef, ira, idec, aperture=1.5),len(idRef)={len(idRef)},len(ids)={len(ids)}")             
                try:
                    url = 'http://12.12.12.251:8888/process_request'
                    params = {"eid": source_file, "reason": f"astrometry is failed,crossmatch len(idRef)={len(idRef)},len(ids)={len(ids)}","stage": 1}
                    response = make_get_request(url, params)
                    return None,None
                except:
                    pass
                logger.info(f"the astrometry is failed, ending this loop")
                return None,None
            
            else:
                isexComd = sexComd_Astro % (iobjFrame_o0, sexConf, sexParam1, ildacn, str(saturate_value))
                logger.info(isexComd)
                sex_comd(isexComd)
                
                ildac = fits.getdata(ildacn, ext=2)
                logger.info(f'#########the star count of detection for astrometry={len(ildac)}')
                ira = ildac["ALPHA_J2000"]
                idec = ildac["DELTA_J2000"] 

                idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)     
                logger.info(f"crossmatch(raRef, decRef, ira, idec, aperture=1.5),len(idRef)={len(idRef)},len(ids)={len(ids)}")
                idra = 3600.0 * (ira[ids] - raRef[idRef])*np.cos(np.radians(decRef[idRef]))
                iddec = 3600.0*(idec[ids]-decRef[idRef])
                iraSig, idecSig = mad(idra), mad(iddec)
                iraMu,idecMu= sigma_clipped_stats(idra)[1],sigma_clipped_stats(iddec)[1]
                #'内部精度'

                ifigname1 = figdir + "astrometryIntra_%s.png" % (ildacn.split("/")[-1][:-5])
                ifigname = figdir + "astrometry_%s.png" % (ildacn.split("/")[-1][:-5])
                mkdir(figdir+'mini/')
                mini_ifigname=figdir+"mini/astrometry_%s.png" % (ildacn.split("/")[-1][:-5])
                cfigname = figdir + "astroCross_%s.png" % (ildacn.split("/")[-1][:-5])
                #logger.info(f'plot len(ra)={len(ira[ids])},len(raRef)={len(raRef[idRef])},len(idra)={len(idra)},len(idRef)={len(idRef)},len(ids)={len(ids)}')
                try:
                    PlotFun.astrometryPlot(idra, iddec, sigma=(iraSig, idecSig), figout=ifigname1)
                    PlotFun.astrometryPlot_new(ira[ids], idec[ids], raRef[idRef], decRef[idRef],idra,iddec, sigma=(iraSig, idecSig),mu=(iraMu,idecMu),figout=ifigname,mini_figout=mini_ifigname)
                    logger.info("[_mphot]:"+'astrometry plotok')
                    try:
                        url = 'http://12.12.12.251:8888/process_request'
                        params = {"eid": source_file, "reason": 'the astrometry process is successful! ',"stage": 1}
                        response = make_get_request(url, params)
                    except:
                        pass
                except Exception as e:
                    logger.error(str(traceback.format_exc()))
                    logger.info(' Astrometric plot generating failed, Astrometry process is failed!')
                    try:
                        url = 'http://12.12.12.251:8888/process_request'
                        params = {"eid": source_file, "reason": 'Astrometric plot generating failed, Astrometry process is failed!',"stage": 1}
                        response = make_get_request(url, params)
                    except:
                        pass
                    return None,None
                
                if(iraSig>0 and idecSig>0 and len(idRef)>10):
                    iimgMat,ihdr=fits.getdata(iobjFrame_o0,header=True) 
                    ihdr['AST_RA'] = iraSig
                    ihdr['AST_DEC'] = idecSig
                    
                    ihdr['AST_CNT'] =len(idRef)
                    fits.writeto(iobjFrame_o0,iimgMat,header=ihdr,overwrite=True)

                logger.info("    Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))
                logger.info("Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))
        except Exception as e:
            logger.error(str(traceback.format_exc()))
            try:
                url = 'http://12.12.12.251:8888/process_request'
                params = {"eid": source_file, "reason": str(traceback.format_exc()),"stage": 1}
                response = make_get_request(url, params)
                return None,None  
            except:
                pass
            return None,None  
            
        logger.info(f"[_mphot]: 2) determine the astrometric accuracy")
        try:
            ################psf phot#####################
            ipsfFrameNew = iobjFrame_o0[:-4] + "psf"
            psfexComd  = "psfex %s -c %s -PSF_SIZE %s,%s -PSF_DIR %s"
            hdul=fits.open(ildacn)###
            hdul[1].data=ildac
            ixcenloc=np.where((ildac["FLAGS"]==0)) 
            ildac1 =ildac[ixcenloc]
             
            hdul[1].data=ildac1
            hdul.close()
             
            fwhm_value=np.nanmedian(ildac1["FWHM_IMAGE"])
            kernel_size= int(fwhm_value)*7
            logger.info(f'##########the psf kernal is 7 size of fwhm={fwhm_value},kernel_size={kernel_size}')
            ipsfComd = psfexComd % (ildacn, psfexConf, str(kernel_size), str(kernel_size), tscidir)
            logger.info(ipsfComd) #增加日志输出，202309081945
            psf_comd(ipsfComd)
            logger.info("[_yphot]:"+"    Generate final SExtractor catalog")
            isexComd = sexComdPSF % (iobjFrame_o0, sexConf, sexParam2, ildacn, str(saturate_value), ipsfFrameNew)
            logger.info(f"sexComdPSF:{isexComd}")
            sex_comd(isexComd)

            ihdr = fits.getheader(iobjFrame_o0)
            if('ASTPRA' in hdr):
               
                hdul = fits.open(ildacn)
                ildac = hdul[2].data
                iheader1 = hdul[1].header
                iheader1['ASTPRA'] = ihdr['ASTPRA']
                iheader1['ASTPRA'] = ihdr['ASTPDEC']
                iheader1['ASTMRA'] = ihdr['ASTMRA']
                iheader1['ASTMDEC'] = ihdr['ASTMDEC']
                iheader1['ASTCNT'] = ihdr['ASTCNT']
                hdul.writeto(iobjFrame_o0[:-5] + "_sexcat.fits", overwrite=True)
                #fits.writeto(iobjFrame_o0[:-5] + "_sexcat.fits",data=hdul, header= overwrite=True)
                
                hdul.close()
            
             
        
        ## some diagnosis charts and files
            ildac = fits.getdata(ildacn, ext=2)
             
            iflag = ildac["FLAGS"]
            isnr = ildac["SNR_WIN"]
            iximg = ildac["X_IMAGE"]
            iyimg = ildac["Y_IMAGE"]
            ira = ildac["ALPHA_J2000"]
            idec = ildac["DELTA_J2000"]

            imagAUTO = ildac["MAG_AUTO"]  
            imagPSF = ildac["MAG_PSF"] 
            imagAPER =ildac["MAG_APER"] 
            imagErrAUTO  = ildac["MAGERR_AUTO"] 
            imagErrPSF  = ildac["MAGERR_PSF"]
            imagErrAPER  = ildac["MAGERR_APER"]


            ifluxAUTO = ildac["FLUX_AUTO"]
            ifluxErrAUTO = ildac["FLUXERR_AUTO"]

            ifluxPSF = ildac["FLUX_PSF"]
            ifluxErrPSF = ildac["FLUXERR_PSF"]

            ifluxAPER =ildac["FLUX_APER"]
            ifluxErrAPER = ildac["FLUXERR_APER"]

            ifwhm = ildac["FWHM_IMAGE"]
            itheta = ildac["THETAWIN_IMAGE"]
            iflux_radius = ildac["FLUX_RADIUS"]

 

        #logger.info(f"[_mphot]:{np.array(ildac["MAG_AUTO"])[1]}")
        except Exception as e:
            logger.error(str(traceback.format_exc()))
            try:
                url = 'http://12.12.12.251:8888/process_request'
                params = {"eid": source_file, "reason": "photometry is failed-------"+str(traceback.format_exc()),"stage": 2}
                response = make_get_request(url, params)
                return None,None
            except:
                pass
    
            return None,None    
        try:
            hdul = fits.open(ildacn)
            ildac = hdul[2].data 
            
            ######fy20230329
            ifluxAUTO_s = ifluxAUTO/ float(exptime)
            ifluxPSF_s = ifluxPSF/ float(exptime)
            ifluxAPER_s = ifluxAPER/ float(exptime) 

            ifluxErrAUTO_s = ifluxErrAUTO/float(exptime)
            ifluxErrPSF_s = ifluxErrPSF/float(exptime)
            ifluxErrAPER_s = ifluxErrAPER/float(exptime)

            imagAUTO_s = imagAUTO + math.log10(exptime) * 2.5
            imagPSF_s = imagPSF + math.log10(exptime) * 2.5 
            imagAPER_s = imagAPER + math.log10(exptime) * 2.5
            
            imagErrAUTO_s = ildac["MAGERR_AUTO"] 
            imagErrPSF_s = ildac["MAGERR_PSF"]
            imagErrAPER_s = ildac["MAGERR_APER"]
            
             
             
            cat = Table(ildac)
            dat = Table({ 'MAG_AUTO_S':imagAUTO_s, 'MAGERR_AUTO_S':imagErrAUTO_s, 'MAG_PSF_S':imagPSF_s, 'MAGERR_PSF_S':imagErrPSF_s,'MAG_APER_S':imagAPER_s, 'MAGERR_APER_S':imagErrAPER_s,
                         'FLUX_AUTO_S':ifluxAUTO_s ,'FLUXERR_AUTO_S':ifluxErrAUTO_s ,'FLUX_PSF_S':ifluxPSF_s,'FLUXERR_PSF_S':ifluxErrPSF_s, 'FLUX_APER_S':ifluxAPER_s,'FLUXERR_APER_S':ifluxErrAPER_s })
            zat = hstack([cat, dat])
            hdul[2].data =np.array(zat)
            
            #hdul.writeto(ildacn, overwrite=True)
            #hdul.writeto(iobjFrame_o0[:-5] + "_sexcat.fits", overwrite=True)
            hdul.writeto(iobjFrame_o0[:-5] + "_sexcat.fits", overwrite=True)
            hdul.close()
            

            #ildac = fits.getdata(ildacn, ext=2)
         
            
        except Exception as e:
            logger.error(str(traceback.format_exc()))
             
            try:
                url = 'http://12.12.12.251:8888/process_request'
                params = {"eid": source_file, "reason": "photometry is failed-------"+str(traceback.format_exc()),"stage": 2}
                response = make_get_request(url, params)
                return None,None
            except:
                pass
            return None,None  
             
         
        ##fy 20230217#####

        logger.info(f"[_mphot]: 1) region file")
         
    

        # 3) determine the zeropoint
        ildac = fits.getdata(ildacn, ext=2)
        ifluxAUTO = ildac["FLUX_AUTO"]
        ifluxErrAUTO = ildac["FLUXERR_AUTO"]
        ifluxPSF = ildac["FLUX_PSF"]
        ifluxErrPSF = ildac["FLUXERR_PSF"]
        ifwhm = ildac["FWHM_IMAGE"]
        iflux_radius = ildac["FLUX_RADIUS"]        
        imagAUTO = ildac["MAG_AUTO"]  
        imagPSF = ildac["MAG_PSF"] 
        imagAPER =ildac["MAG_APER"] 
        imagErrAUTO  = ildac["MAGERR_AUTO"] 
        imagErrPSF  = ildac["MAGERR_PSF"]
        imagErrAPER  = ildac["MAGERR_APER"]         

        iximg = ildac["X_IMAGE"]
        iyimg = ildac["Y_IMAGE"]
        iflag = ildac["FLAGS"]
        isnr = ildac["SNR_WIN"]

        ira = ildac["ALPHA_J2000"]
        idec = ildac["DELTA_J2000"] 
    
        ifwhm = ildac["FWHM_IMAGE"]
        itheta = ildac["THETAWIN_IMAGE"]
         
        psf_id = iflux_radius > 1.5
        # select high-snr stars
        gid = (iflag == 0) & (isnr > 30.0)
        
        ira = ira[gid]
        idec = idec[gid]
        imagPSF = imagPSF[gid]

        ifwhm = ifwhm[gid]
        iximg = iximg[gid]
        iyimg = iyimg[gid]
        itheta = itheta[gid]
        idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)
        ipsffig = figdir + "psfPhot_%s.pdf" % (ildacn.split("/")[-1][:-5])
        izpArr = magRef[idRef] - imagPSF[ids]
        imask = sigma_clip(izpArr, sigma=3.0, maxiters=3, masked=True)
        izpArr = izpArr[~imask.mask]
        #izp = np.median(izpArr)
        izp=0
        try:
            if (len(izpArr) < 20):
                nameonj = iobjFrame_o0.rsplit('/', 1)[1]
                logger.error("[_mphot]:"+'XXXXXXXXX' + nameonj + ':  the astrometric measurement may not success, please check all print information and source fits!'+ f"len(izpArr):{len(izpArr)}")

            bkgNoise = np.sqrt(np.pi) * np.median(ifwhm) * (np.min([ihdr["SKYRMS1"], ihdr["SKYRMS2"]]))
            PlotFun.mag2Err(ifluxPSF[psf_id], ifluxErrPSF[psf_id], ifluxAUTO[psf_id], ifluxErrAUTO[psf_id],
                            bkgNoise=bkgNoise, zpoint=izp, gain=ihdr["GAIN"], figout=ipsffig)
            loggerloguru.info(iobjFrame_o0+' has been successfully processed!')
        except Exception as e:
                logger.error(str(traceback.format_exc()))
                try:
                    url = 'http://12.12.12.251:8888/process_request'
                    params = {"eid": source_file, "reason": "psf photometry is failed-------"+str(traceback.format_exc()),"stage": 2}
                    response = make_get_request(url, params)
                except:
                    pass
                return None,None
                logger.info(' psfphot generating failed!')

        
        
            
        # 4) show the Whisker plot
        # try:
        #     iwskfig     = figdir + "ewhiskerEll_%s.pdf"%(ildacn.split("/")[-1][:-5])
        #     iwskfwhmfig = figdir + "ewhiskerFWHM_%s.pdf"%(ildacn.split("/")[-1][:-5])
        #     PlotFun.ewhisker(iximg, iyimg, iaimg, ibimg, itheta, figout=iwskfig)
        #     #plotFun.ewhiskerEllipse(iximg, iyimg, iaimg, ibimg, itheta, scale=30.0, figout=iwskfig)
        #     PlotFun.ewhiskerFWHM(iximg,iyimg,ifwhm,itheta,scale=20,figout=iwskfwhmfig)
        # except Exception as e:
        #     logger.info(str(e))
        #     logger.info('the Whisker figure can not be generated!') 
        #     logger.info('XXXXXXXXXXXXXXXXX the Whisker figure can not be generated!')      
        #     continue    

        # 5) update image header
        # ihdr["ZPMAG"] = izp
        # ihdr["ZPERR"] = izpErr
 
 
            
        # except Exception as e:
        #     loggerloguru.info(str(e))
        #     loggerloguru.info('the zpPlot figure can not be generated!')         
 
          
    except Exception as e:
            
            logger.info(traceback.format_exc())
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": source_file, "reason": str(traceback.format_exc()),"stage": 2}
            response = make_get_request(url, params)
            return None, None
          
    
    
    return ildacn,iobjFrame_o0

def scamp_comd(scampcomd):
     
    iscampComdlist=[]
    iscampComdlist=scampcomd.split(" ")
    try:
        # 执行命令并设置超时时间为10秒
        completed_process_scampcomd = subprocess.run(iscampComdlist, timeout=300, shell=False, capture_output=True)
        #logger.info(iscampComdlist)           
        # 处理命令执行结果
        if completed_process_scampcomd.returncode == 0:
           
            logger.info(completed_process_scampcomd.stdout)
        else:
            logger.info(completed_process_scampcomd.stderr)
            logger.info('XXXXXXXXXXXXXXXscamp failedXXXXXXXXXXXXXXXX')
    
    except subprocess.TimeoutExpired:
        logger.info("Command scamp_comd execution timed out.")
        return None, None
    
def sex_comd(isexComd):    
    sexComdlist=[]
    sexComdlist=isexComd.split(" ")
    #logger.info(sexComdlist)
    try:
    # if(1==1):
    #     os.system(isexComd)
        # 执行命令并设置超时时间为10秒
        
        completed_process_isexcomd = subprocess.run(sexComdlist, timeout=900, shell=False, capture_output=True)
                         
                         
    # 处理命令执行结果
        if completed_process_isexcomd.returncode == 0:
           
            logger.info(completed_process_isexcomd.stdout)
        else:
            logger.info('XXXXXXXXXXXXXXXXXSEXTRACTOR')
            logger.info(completed_process_isexcomd.stderr)
            return 0

    except subprocess.TimeoutExpired:
        logger.info("Command sex_comd execution timed out.")  
        return None, None
    
def psf_comd(ipsfComd):
    ipsfComdlist=[]
    #os.system(ipsfComd)
    ipsfComdlist=ipsfComd.split(" ")
    objfilename = ipsfComdlist[1] 
    #objfilename = ntpath.basename(objfilename)
    objfilename =objfilename.split('_sciimg')[0]+'.fits'
    logger.info('*****************psfcomd************'+ objfilename)
    #logger.info(ipsfComdlist)
    try:
        # 执行命令并设置超时时间为10秒
        completed_process_ipsfcomd = subprocess.run(ipsfComdlist, timeout=900, shell=False, capture_output=True)
                        
        # 处理命令执行结果
        if completed_process_ipsfcomd.returncode == 0:
            logger.info(completed_process_ipsfcomd .stdout)
        else:
            logger.info(completed_process_ipsfcomd.stderr)
            try:
                url = 'http://12.12.12.251:8888/process_request'
                params = {"eid": objfilename, "reason": "psf photometry is failed"+str(completed_process_ipsfcomd.stderr),"stage": 2}
                response = make_get_request(url, params)
            except:
                pass
    
    except subprocess.TimeoutExpired:
        logger.info("Command psf_comd execution timed out.")
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": objfilename, "reason": "photometry is failed for time out","stage": 2}
            response = make_get_request(url, params)
        except:
            pass
        return None, None
    
def astrometry(ildacn,refCatn,iobjFrame_o0,scampComd, scampConf, scampConf2,figdir):
    try:
        iscampComd = scampComd%(ildacn,scampConf,refCatn)
        #os.system(iscampComd)
        logger.info(iscampComd)
        scamp_comd(iscampComd)


        ihead   = ildacn[:-4] + "head"
        ihead1 = ildacn[:-4] + "ahead"


        os.system("mv %s %s"%(ihead,ihead1))
        iscampComd2 = scampComd%(ildacn,scampConf2,refCatn)
         
        logger.info(iscampComd2)
        scamp_comd(iscampComd2)
 
    except Exception as e:
        logger.error('scamp is calling failed!')
         
    ihead = ildacn[:-4] + "head"
    if(os.path.exists(ihead)):
        
        iimgMat, ihdr = fits.getdata(iobjFrame_o0, header=True)

        with open(ihead, "r") as iheadList_fits:
            
            iheadList=iheadList_fits.read().splitlines()[3:-1]
            nhead = len(iheadList)
            #logger.info(iheadList)
            for ih in range(nhead):
                iheadx = iheadList[ih].split("/")[0].split("=")
                try:
                    ihdr[iheadx[0]] = float(iheadx[1])
                except:
                    ihdr[iheadx[0]] = str(iheadx[1])
        ihdr["CTYPE1"] = "RA---TPV"
        ihdr["CTYPE2"] = "DEC--TPV"
        ihdr["RADESYS"] = "ICRS"
        ihdr["CUNIT1"] = "deg"
        ihdr["CUNIT2"] = "deg"
        fits.writeto(iobjFrame_o0, iimgMat, ihdr, overwrite=True)
        logger.info(ihdr["CTYPE1"])
        logger.info(ihdr["CTYPE2"])
        try:
            if (os.path.exists(ihead)):
                os.system("rm %s" % ihead)
            if (os.path.exists(ihead1)):
                os.system("rm %s" % ihead1)
            os.system("mv *pdf %s" % (figdir)) #这个存在报错mv: cannot stat ‘*pdf’: No such file or directory
        except:
            logger.error('there is no any pdf file')
    else:
        logger.error(f"os.path.exists(ihead)==False: {os.path.exists(ihead)==False}")
        return None,None


def photometry(rootpath,ttfname,date):
    
    upath=rootpath+'reception/'+str(date)+'/'
    scipath = upath+'sci/'+ttfname+'/'
    logger.info(scipath+'fileguide/'+ttfname+'_sciimg.npy')
    if not(os.path.exists(scipath+'fileguide/'+ttfname+'_sciimg.npy')):
        logger.info("There is no such observation: %s" % (ttfname)) 
        return None
    lpathsciimg=np.load(scipath+'fileguide/'+ttfname+'_sciimg.npy')
        
    logger.info(lpathsciimg)
    lpathsciimg=np.sort(lpathsciimg)
    #basic setting
    #dirname, filename = os.path.split(os.path.abspath(__file__))
    imgdir  = rootpath + "images/"
    confdir = rootpath + "config/"
    ancdir  = rootpath + "ancData/"
    #figdir  = rootpath + "figures/"+str(date)+'/'+str(date)+'_'+ttfname+"/"
    pdfdir = rootpath+'50cmpy/'
    logdir = rootpath+'run_pipe_log/'+str(date)+'/'+ttfname+'/'
    ttf=ttfname.split('_')
    tid=ttf[0]
    targetid=ttf[1]
    figdir= '/home/50cm/dugking/CrossPic/'+date+'/'+tid+'/'
    upath=rootpath+'reception/'+str(date)+'/'
    mkdir(logdir)
    mkdir(figdir)   
    mkdir(upath)

 
    

    if len(lpathsciimg)==0:
        logger.error (f"lpathsciimg:{lpathsciimg}")
        return None
    
    ira='0'
    idec='0'
    
    for i in range(0,len(lpathsciimg)):
        # subpath = datename(lpathsciimg)
    
        # if(subpath is None):
        #     logger.error('the file can not be open successfully')
        #     return None
        # tscidir=rootpath+'reception/'+str(date)+'/sci/'+ttfname+'/'+subpath+'/'
        # mkdir(tscidir)
        try:
            hrs,imgs=readfits(lpathsciimg[i])
            #logger.info(hrs)
            if "OBJCTRA" in hrs.keys():
                ira  = ":".join(hrs["OBJCTRA"].split())
                idec = ":".join(hrs["OBJCTDEC"].split())
                ira, idec = d2hms(ira, idec, conv=1)
                ira =str(ira)
                idec=str(idec)
            elif "RA" in hrs.keys(): 
                #logger.info(hrs["RA"],hrs["DEC"])
                ira  = ":".join(hrs["RA"].split())
                idec = ":".join(hrs["DEC"].split())
                # logger.info('split --- idec')
                # logger.info(idec)
                ira =str(ira)
                idec=str(idec)
                ira, idec = d2hms(ira, idec, conv=1)
            logger.info(f'ira={ira}, idec={idec}')
        except:
            logger.error(str(traceback.format_exc()))
            continue
             
    
    if(ira=='0' or idec=='0'):
        logger.info('these are no parameters of object ra and dec in header of this fits')
        return None
    # else:
    #     ira, idec = d2hms(ira, idec, conv=1)
 
     
    r,c=np.shape(imgs)
    filtID  = ttf[2][-1]
    objID   = ttf[1] 

    binfact = int(hrs['XBINNING'])
    pscale  = 0.297 * binfact # pixel scale in unit of arcsec
    
    # saturate_value=hrs['SATURATE']*0.9
    # exptlast = hrs['EXPTIME']
    # refCatn  = ancdir + "GaiaStar_%s_%s.ldac"%(objID,date)
    # refCat2n = ancdir + "GaiaStar_%s_%s.cat"%(objID,date)
    #if not (os.path.exists(refCatn)):
    # maxmaglim=ref_expt_select(confdir,filtID,binfact,exptlast)
    # logger.info(f'ira={ira}, idec={idec}')
    # refSed_gaia(rootpath,ttfname,date,ira,idec,10,maxmaglim)
    # refCat2  = Table.read(refCat2n, format="ascii")
    # raRef    = refCat2["ra"]
    # decRef   = refCat2["dec"]
    # try:
    #     # magRef   = refCat2["mag%s"%filtID.upper()]
    #      magRef  = refCat2["mag"]
    # except:
    #      logger.info('the filteris is '+filtID.upper()+' and use mag_nofilter only') 
    #      magRef  = refCat2["mag"]    
    # nrefStar = len(raRef)

    scipath=rootpath+'reception/'+str(date)+'/sci/'+ttfname+'/' 

     # 3S configuration
    sexParam1  = confdir + "default.param"
    sexParam2  = confdir + "default.paramNew"
    sexConf    = confdir + "default.sex"
    swarpConf  = confdir + "default.swarp"
    scampConf  = confdir + "default.scamp"
    scampConf2  = confdir + "second.scamp"
    psfexConf  = confdir + "default.psfex"
    
    
    miniarea = str(6/binfact)
    sexComd1   = "sex %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s"
    sexComd2_astro   = " -DETECT_MINAREA "  +miniarea+ " -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5"
    sexComd2   = " -DETECT_MINAREA "  +miniarea+ " -DETECT_THRESH 1.0 -ANALYSIS_THRESH 1.5"
    sexComd3   = " -SATUR_LEVEL %s -CHECKIMAGE_TYPE NONE"
    sexComd4   = " -PSF_NAME %s"
    sexComd    = sexComd1 + sexComd2 + sexComd3
    sexComd_Astro    = sexComd1 + sexComd2_astro + sexComd3
    sexComdPSF = sexComd + sexComd4
    scampComd1 = "scamp %s -c %s -ASTREFCAT_NAME %s -MOSAIC_TYPE LOOSE" #UNCHANGED" # LOOSE "
    scampComd2 = " -FWHM_THRESHOLDS 2,20 -SN_THRESHOLDS 20,1000"
    scampComd  = scampComd1 + scampComd2
    

    psfexComd  = "psfex %s -c %s -PSF_SIZE %d,%d -PSF_DIR %s"

    # logger.info('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    # iscampComd_test = scampComd%('ildacn',scampConf2,refCatn)
    # logger.info(iscampComd_test)
    # logger.info('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    # set logging
    skey  = str(date)+'_'+ttfname+'_phos'
    loggn = logdir + skey + ".log"
    if os.path.exists(loggn): os.system("rm %s"%loggn)
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logFile = logging.FileHandler(loggn)
    logFile.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger(skey).addHandler(logFile)
    loggerpho = logging.getLogger(skey)
    loggerpho.info("^_^ Telescope_filter_target: %s"%(ttfname))
    

    for i,file_name in enumerate(lpathsciimg):
        filei=lpathsciimg[i].rsplit('/',1)[1]
        #logger.info('####### the precessed fits is: '+ filei)
        loggerpho.info('####### the processed fits is: %s'%(filei))
        #对每一个frame run sextractor
        subpath = datename(lpathsciimg[i]) 
        tscidir = scipath +  subpath +'/'
        mkdir(tscidir)
        single_cmd_deal(lpathsciimg[i],ttfname,date,sexComd,sexConf,sexParam1,psfexComd,psfexConf,scampComd, scampConf, scampConf2,sexComdPSF,sexParam2,rootpath,figdir,ancdir,tscidir,confdir)
        #single_cmd_deal(i_obj_sci,i_obj_file,ttfname,date,sexComd,sexComdAstro,sexConf,sexParam1,psfexComd,psfexConf,scampComd, scampConf, scampConf2,sexComdPSF,sexParam2,rootpath,figdir,ancdir,tscidir,confdir)
        #i_obj_sexcat =tscidir+ filename[:-5]+'_sciimg_sexcat.fits'
    loggerpho.info("^_^ ALL astrometry and photometry processes are ok!")


def pro(tid_target_filter,date):
    rootpath=get_current_dir()
    photometry(rootpath,tid_target_filter,date)
    logger.info('the astrometric process is done')

 

if __name__ == "__main__":
    import argparse
    import time
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20221215', help='输入处理日期')
    parser.add_argument('--tid_target_filter', default='y50b_M33_mr',type=str, help='')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_ycom").get_logger
    starttime=time.time()
    pro(args.tid_target_filter,args.date)
    endtime=time.time()
    loggerloguru.info("[_yphostcom1] "+"time spending:{}s".format(endtime-starttime))
