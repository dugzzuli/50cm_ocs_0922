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
from astropy.stats import sigma_clip
from astropy.table import Table,hstack
from astropy.time import Time
from loguru import logger
from scipy.optimize import leastsq
from scipy.spatial import cKDTree as ckdt
 
import lib.phot.PlotFun as PlotFun
from execute.flatdiff import mkdir
from lib.LogInstance import Logger

from lib.phot.yrefsed import refSed_gaia
import subprocess
import config

from execute.yphotutils import datename,d2hms,mad,wds9reg,read_param,crossmatch,reg,ref_expt_select,readfits,hdr_edit_fits,re_hdr_edit_fits,read_list,timeplus,HJDTrans,get_scamp_head,overlapRect,pointRect,_spherical_to_cartesian,_great_circle_distance



  
def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)
  



def single_cmd_deal(i_obj_file,ttfname,date,sexComd,sexConf,sexParam1,psfexComd,psfexConf,scampComd, scampConf, scampConf2,sexComdPSF,sexParam2,rootpath,figdir,ancdir,tscidir,configdir):
    loggerloguru = Logger(date, '', "_yphost").get_logger
    logger.info('####### the processed fits is: %s' % (i_obj_file))
    ildacn,iobjFrame_o0=None,None
    #try:
    if(1==1):
        # 对每一个frame run sextractor
         
        iobjFrameNew = i_obj_file.split('sciimg')[0] + 'fcimg.fits'
        iobjFrame_o0 = i_obj_file

        hdr = fits.getheader(iobjFrameNew )
         
        logger.info(iobjFrameNew)
        # if os.path.exists( iobjFrame_o0[:-5] + "_sexcat.fits") and config["log"]["record"]:
        #     #print(iobjFrame_o0[:-5] + '_sexcat.fits'+' was been phot processd! ')
        #     ildacn = iobjFrame_o0[:-4] + "ldac"
        #     return ildacn,iobjFrame_o0

        logger.info(iobjFrame_o0)
        ################################################
        ipsfFrameNew = iobjFrame_o0[:-4] + "psf"
         
        
        
        try:
            ihdr,iimgMat= readfits(iobjFrame_o0)
            r,c=np.shape(iimgMat)
        except:
            logger.info(f'the fits file is broken:{iobjFrameNew}')
            logger.info("XXXXXXXXXXXXXXXXX the fits file is broken: %s"%(iobjFrameNew))
            return None, None
 
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
        exptime=ihdr['EXPTIME']
        ildacn = iobjFrame_o0[:-4] + "ldac"
         
        saturate_value=ihdr['SATURATE']*0.9
        exptlast = ihdr['EXPTIME']
        tid, objID, filterid=ttfname.split('_')
        
        isexComd = sexComd % (iobjFrameNew, sexConf, sexParam1, ildacn, str(saturate_value))
        sex_comd(isexComd)
        
        refCatn  = ancdir + "GaiaStar_%s_%s.ldac"%(objID,date)
        refCat2n = ancdir + "GaiaStar_%s_%s.cat"%(objID,date)

        if not (os.path.exists(refCatn)):#(1==1):
 
            maxmaglim=ref_expt_select(configdir,filterid,binfact,exptlast)
            logger.info(f'ira={ira}, idec={idec}')
            refSed_gaia(rootpath,ttfname,date,ira,idec,10,maxmaglim)
            
        
        # try:
        #     # magRef   = refCat2["mag%s"%filtID.upper()]
        #     refCat2  = Table.read(refCat2n, format="ascii")
        #     raRef    = refCat2["ra"]
        #     decRef   = refCat2["dec"]
        #     magRef  = refCat2["mag"]
        # except:
        #     logger.info('the filteris is '+filterid.upper()+' and use mag_nofilter only')
        #     magRef  = refCat2["mag"]    
        #nrefStar = len(raRef)
        
        
        refCat2  = Table.read(refCat2n, format="ascii")
        raRef    = refCat2["ra"]
        decRef   = refCat2["dec"]
        magRef  = refCat2["mag"]
        
        ildac=fits.getdata(ildacn,ext=2)
            
        #if(filterid =='v' or filterid =='u'):
        try:
            ixcenloc=np.where((ildac["X_IMAGE"]>c/2+10) | (ildac["X_IMAGE"]<c/2-10))
            #ixyloc=np.where((iximg<c-(1000/binfact))&(iximg>1000/binfact))#&(iyimg<6100)&(iyimg>30))

            ildac =ildac[ixcenloc]
            hdul=fits.open(ildacn)###
            hdul[1].data=ildac
            hdul.close()
        except: 
            logger.info('XXXXXXXXXXXXXXXXXXXX  ldac cut failed ')
        logger.info('astrometry begin')
        astrometry(ildacn,refCatn,iobjFrame_o0,scampComd, scampConf, scampConf2 ,figdir)
        
        
        try:
            isexComd = sexComd % (iobjFrame_o0, sexConf, sexParam1, ildacn, str(saturate_value))
            sex_comd(isexComd)

            ildac = fits.getdata(ildacn, ext=2)
            ira = ildac["ALPHA_J2000"]
            idec = ildac["DELTA_J2000"] 
            idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)     
            if len(idRef)==0:
                logger.info(f"astrometry is failed!{iobjFrame_o0}")
                logger.error(f"astrometry is failed! crossmatch(raRef, decRef, ira, idec, aperture=1.5),len(idRef)={len(idRef)},len(ids)={len(ids)}")
                logger.info(f"astrometry is failed! crossmatch(raRef, decRef, ira, idec, aperture=1.5),len(idRef)={len(idRef)},len(ids)={len(ids)}")
                re_hdr_edit_fits(iobjFrame_o0)
                isexComd = sexComd % (iobjFrame_o0, sexConf, sexParam1, ildacn, str(saturate_value))
                sex_comd(isexComd)
                astro=astrometry(ildacn,refCatn,iobjFrame_o0,scampComd, scampConf, scampConf2,configdir)
                if astro is None:
                    #logger.error(f"os.path.exists(ihead)==False: {os.path.exists(ihead)==False}")
                    logger.info(f"ending this loop")
                    return None,None
        except Exception as e:
            logger.error(str(traceback.format_exc()))
            return None,None        
        logger.info(f"[_mphot]: 2) determine the astrometric accuracy")
     
        try:
            ################psf phot#####################
            ipsfFrameNew = iobjFrame_o0[:-4] + "psf"
            psfexComd  = "psfex %s -c %s -PSF_SIZE %s,%s -PSF_DIR %s"
            # isexComd = sexComd % (iobjFrame_o0, sexConf, sexParam1, ildacn, saturate_value)
            # os.system(isexComd)
            ipsfComd = psfexComd % (ildacn, psfexConf, str(21), str(21), tscidir)
            psf_comd(ipsfComd)
            
            logger.info("[_mphot]:"+"    Generate final SExtractor catalog")
            isexComd = sexComdPSF % (iobjFrame_o0, sexConf, sexParam2, ildacn, str(saturate_value), ipsfFrameNew)
            sex_comd(isexComd)
            
            hdul = fits.open(ildacn)
            ildac = hdul[2].data 
            hdul.writeto(iobjFrame_o0[:-5] + "_sexcat.fits", overwrite=True)
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
            imagErrPSF = ildac["MAGERR_PSF"]
            ifluxAUTO = ildac["FLUX_AUTO"]
            ifluxErrAUTO = ildac["FLUXERR_AUTO"]
            ifluxPSF = ildac["FLUX_PSF"]
            ifluxErrPSF = ildac["FLUXERR_PSF"]
            ifwhm = ildac["FWHM_IMAGE"]
            itheta = ildac["THETAWIN_IMAGE"]
            iflux_radius = ildac["FLUX_RADIUS"]
        #logger.info(f"[_mphot]:{np.array(ildac["MAG_AUTO"])[1]}")
        except Exception as e:
            logger.error(str(traceback.format_exc()))
            return None,None    
        try:
            hdul = fits.open(ildacn)
            ildac = hdul[2].data 
            
            ######fy20230329
            imagAUTO_s = ildac["MAG_AUTO"] + math.log10(exptime) * 2.5
            imagPSF_s = ildac["MAG_PSF"] + math.log10(exptime) * 2.5 
            imagAPER_s = ildac["MAG_APER"] + math.log10(exptime) * 2.5
            
            imagErrAUTO_s = ildac["MAGERR_AUTO"] 
            imagErrPSF_s = ildac["MAGERR_PSF"]
            imagErrAPER_s = ildac["MAGERR_APER"]
            
            len_ldac = len(ildac["FLUXERR_AUTO"])
            
            for i in range(0,len_ldac):
                if(float(ildac["FLUXERR_AUTO"][i])==0.0):
                    imagErrAUTO_s[i] =  ildac["MAGERR_AUTO"][i]
                else:
                    imagErrAUTO_s[i] = 1.0857 * (float(ildac["FLUXERR_AUTO"][i])/float(ildac["FLUX_AUTO"][i]))
                    
                if(float(ildac["FLUXERR_APER"][i][0])==0.0): 
                    imagErrAPER_s[i] =ildac["MAGERR_APER"][i]
                else:
                    imagErrAPER_s[i] = 1.0857 * (ildac["FLUXERR_APER"][i]/ildac["FLUX_APER"][i])
                try:
                    if(float(ildac["FLUXERR_PSF"][i])==0.0):
                        imagErrPSF_s[i]=  ildac["MAGERR_PSF"][i]
                    else:
                        imagErrPSF_s[i]= 1.0857 * (float(ildac["FLUXERR_PSF"][i])/float(ildac["FLUX_PSF"][i]))
                except Exception as e:
                        logger.error(str(traceback.format_exc()))
 
            
          
            ifluxAUTO_s = ifluxAUTO/ float(exptime)
            ifluxPSF_s = ifluxPSF/ float(exptime)
            ifluxAPER_s = ildac["FLUX_APER"]/ float(exptime) 
            
                
            cat = Table(ildac)
            dat = Table({ 'MAG_AUTO_S':imagAUTO_s, 'MAGERR_AUTO_S':imagErrAUTO_s, 'MAG_PSF_S':imagPSF_s, 'MAGERR_PSF_S':imagErrPSF_s,'MAG_APER_S':imagAPER_s, 'MAGERR_APER_S':imagErrAPER_s,
                         'FLUX_AUTO_S':ifluxAUTO_s ,'FLUX_PSF_S':ifluxPSF_s, 'FLUX_APER_S':ifluxAPER_s })
            zat = hstack([cat, dat])
            hdul[2].data =np.array(zat)
             
            #hdul.writeto(ildacn, overwrite=True)
            hdul.writeto(iobjFrame_o0[:-5] + "_sexcat.fits", overwrite=True)
            hdul.close()
            
            ildac = fits.getdata(ildacn, ext=2)
         
            imagAUTO = ildac["MAG_AUTO"]  
            imagPSF = ildac["MAG_PSF"] 
            ifluxAUTO = ildac["FLUX_AUTO"]
            ifluxErrAUTO = ildac["FLUXERR_AUTO"]
            ifluxPSF = ildac["FLUX_PSF"]
            ifluxErrPSF = ildac["FLUXERR_PSF"]
        except Exception as e:
            logger.error(str(traceback.format_exc()))
            return None,None    
         
        ##fy 20230217#####

        logger.info(f"[_mphot]: 1) region file")
         
    # 1) region file
        iregn  = iobjFrame_o0[:-4] + "reg"
        ildac_select_reg = iobjFrame_o0[:-5] + "_mag.reg"
        wds9reg(ira,idec,flag=None,radius=8.0,unit="arcsec",color="green",outfile=iregn)
        reg(ildacn, ildac_select_reg)
        psf_id    = iflux_radius>1.5
        # select high-snr stars
        gid        = (iflag==0) & (isnr>30.0)
        ira        = ira[gid]
        idec       = idec[gid]
        imagAUTO   = imagAUTO[gid]
        imagPSF    = imagPSF[gid]
        imagErrPSF = imagErrPSF[gid]
        ifwhm      = ifwhm[gid]
        iximg      = iximg[gid]
        iyimg      = iyimg[gid]
        # iaimg      = iaimg[gid]
        # ibimg      = ibimg[gid]
        itheta     = itheta[gid]
 
        idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)

        # 2) determine the astrometric accuracy
         
        idra = 3600.0 * (ira[ids] - raRef[idRef])*np.cos(np.radians(decRef[idRef]))
        iddec = 3600.0*(idec[ids]-decRef[idRef])
        iraSig, idecSig = mad(idra), mad(iddec)
        #'内部精度'

        ifigname1 = figdir + "astrometryIntra_%s.png" % (ildacn.split("/")[-1][:-5])
        ifigname = figdir + "astrometry_%s.png" % (ildacn.split("/")[-1][:-5])
        mkdir(figdir+'mini/')
        mini_ifigname=figdir+"mini/astrometry_%s.png" % (ildacn.split("/")[-1][:-5])
        cfigname = figdir + "astroCross_%s.png" % (ildacn.split("/")[-1][:-5])
        #logger.info(f'plot len(ra)={len(ira[ids])},len(raRef)={len(raRef[idRef])},len(idra)={len(idra)},len(idRef)={len(idRef)},len(ids)={len(ids)}')
        try:
            PlotFun.astrometryPlot(idra, iddec, sigma=(iraSig, idecSig), figout=ifigname1)
            PlotFun.astrometryPlot_new(ira[ids], idec[ids], raRef[idRef], decRef[idRef],idra,iddec, sigma=(iraSig, idecSig), figout=ifigname,mini_figout=mini_ifigname)
            logger.info("[_mphot]:"+'astrometry plotok')
        except Exception as e:
                logger.error(str(traceback.format_exc()))
                logger.info(' Astrometric scatter generating failed!')
         
        if(iraSig>0 and idecSig>0 and len(idRef)>0):
            iimgMat,ihdr=fits.getdata(iobjFrame_o0,header=True)
            ihdr['AST_RA'] = iraSig
            ihdr['AST_DEC'] = idecSig
            ihdr['AST_CNT'] =len(idRef)
            fits.writeto(iobjFrame_o0,iimgMat,header=ihdr,overwrite=True)

        logger.info("    Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))
        logger.info("Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))

        # 3) determine the zeropoint
        ildac = fits.getdata(ildacn, ext=2)
            
        # iflag = ildac["FLAGS"]
        # isnr = ildac["SNR_WIN"]
        # iximg = ildac["X_IMAGE"]
        # iyimg = ildac["Y_IMAGE"]
        # ira = ildac["ALPHA_J2000"]
        # idec = ildac["DELTA_J2000"]
        # imagAUTO = ildac["MAG_AUTO"] 
        # imagPSF = ildac["MAG_PSF"] 
        # ifluxAUTO = ildac["FLUX_AUTO"]
        # ifluxErrAUTO = ildac["FLUXERR_AUTO"]
        # ifluxPSF = ildac["FLUX_PSF"]
        # ifluxErrPSF = ildac["FLUXERR_PSF"]
        # ifwhm = ildac["FWHM_IMAGE"]
        # itheta = ildac["THETAWIN_IMAGE"]
        # iflux_radius = ildac["FLUX_RADIUS"]
        
        # psf_id = iflux_radius > 1.5
        # # select high-snr stars
        # gid = (iflag == 0) & (isnr > 30.0)
        
        # ira = ira[gid]
        # idec = idec[gid]
        # imagPSF = imagPSF[gid]

        # ifwhm = ifwhm[gid]
        # iximg = iximg[gid]
        # iyimg = iyimg[gid]
        # itheta = itheta[gid]
        #idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)
        ipsffig = figdir + "psfPhot_%s.pdf" % (ildacn.split("/")[-1][:-5])
        izpArr = magRef[idRef] - imagPSF[ids]
        imask = sigma_clip(izpArr, sigma=3.0, maxiters=3, masked=True)
        izpArr = izpArr[~imask.mask]
        #izp = np.median(izpArr)
        izp=0
        try:
            if (len(izpArr) < 20):
                nameonj = iobjFrame_o0.rsplit('/', 1)[1]
                logger.info("[_mphot]:"+'XXXXXXXXX' + nameonj + ':  the astrometric measurement may not success, please check all print information and source fits!')

            bkgNoise = np.sqrt(np.pi) * np.median(ifwhm) * (np.min([ihdr["SKYRMS1"], ihdr["SKYRMS2"]]))
            PlotFun.mag2Err(ifluxPSF[psf_id], ifluxErrPSF[psf_id], ifluxAUTO[psf_id], ifluxErrAUTO[psf_id],
                            bkgNoise=bkgNoise, zpoint=izp, gain=ihdr["GAIN"], figout=ipsffig)
            loggerloguru.info(iobjFrame_o0+' has been successfully processed!')
        except Exception as e:
                logger.error(str(traceback.format_exc()))
                logger.info(' psfphot generating failed!')
        
         
            
        # 4) show the Whisker plot
        # try:
        #     iwskfig     = figdir + "ewhiskerEll_%s.pdf"%(ildacn.split("/")[-1][:-5])
        #     iwskfwhmfig = figdir + "ewhiskerFWHM_%s.pdf"%(ildacn.split("/")[-1][:-5])
        #     plotFun.ewhisker(iximg, iyimg, iaimg, ibimg, itheta, figout=iwskfig)
        #     #plotFun.ewhiskerEllipse(iximg, iyimg, iaimg, ibimg, itheta, scale=30.0, figout=iwskfig)
        #     plotFun.ewhiskerFWHM(iximg,iyimg,ifwhm,itheta,scale=20,figout=iwskfwhmfig)
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
 
          
    # except Exception as e:
    #         logger.info(str(e))
    #         logger.info('the phot process is failed !')      
    # finally:
    #     for k in list(locals().keys()):
    #             # if locals[k] is np.nan:
    #             try:
    #                 del locals[k]
    #             except:
    #                 continue
    #     gc.collect()
          
    
    
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

    except subprocess.TimeoutExpired:
        logger.info("Command sex_comd execution timed out.")  
        return None, None
    
def psf_comd(ipsfComd):
    ipsfComdlist=[]
    #os.system(ipsfComd)
    ipsfComdlist=ipsfComd.split(" ")
    #logger.info(ipsfComdlist)
    try:
        # 执行命令并设置超时时间为10秒
        completed_process_ipsfcomd = subprocess.run(ipsfComdlist, timeout=900, shell=False, capture_output=True)
                        
        # 处理命令执行结果
        if completed_process_ipsfcomd.returncode == 0:
            logger.info(completed_process_ipsfcomd .stdout)
        else:
            logger.info(completed_process_ipsfcomd.stderr)
    
    except subprocess.TimeoutExpired:
        logger.info("Command psf_comd execution timed out.")
        return None, None
    
def astrometry(ildacn,refCatn,iobjFrame_o0,scampComd, scampConf, scampConf2,figdir):
    try:
        iscampComd = scampComd%(ildacn,scampConf,refCatn)
        #os.system(iscampComd)
        
        scamp_comd(iscampComd)


        ihead   = ildacn[:-4] + "head"
        ihead1 = ildacn[:-4] + "ahead"


        os.system("mv %s %s"%(ihead,ihead1))
        iscampComd2 = scampComd%(ildacn,scampConf2,refCatn)
         
        #os.system(iscampComd2)
        scamp_comd(iscampComd2)

        

        ihead   = ildacn[:-4] + "head"
        ihead1 = ildacn[:-4] + "ahead"

        os.system("mv %s %s"%(ihead,ihead1))
        iscampComd2 = scampComd%(ildacn,scampConf2,refCatn)
        #os.system(iscampComd2)
        scamp_comd(iscampComd2)
    except Exception as e:
        logger.error('scamp is calling failed!')
         
    
        # if(os.path.exists(ihead)==False):
        #     logger.error(f"os.path.exists(ihead)==False: {os.path.exists(ihead)==False}")
        #     logger.info(f"ending this loop")

        #     return None,None
    ihead = ildacn[:-4] + "head"
    if(os.path.exists(ihead)):
        
        iimgMat, ihdr = fits.getdata(iobjFrame_o0, header=True)

        with open(ihead, "r") as iheadList_fits:
            
            iheadList=iheadList_fits.read().splitlines()[3:-1]
            nhead = len(iheadList)

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
            os.system("rm %s" % ihead)
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
    figdir  = rootpath + "figures/"+str(date)+'/'+str(date)+'_'+ttfname+"/"
    pdfdir = rootpath+'50cmpy/'
    logdir = rootpath+'run_pipe_log/'+str(date)+'/'+ttfname+'/'
    
    mkdir(logdir)
    mkdir(figdir)   
    ttf=ttfname.split('_')
    tid=ttf[0]
    targetid=ttf[1]
    if len(lpathsciimg)==0:
        logger.error (f"lpathsciimg:{lpathsciimg}")
        return None
    
    ira='0'
    idec='0'
    
    for i in range(0,len(lpathsciimg)):
        
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
    sexComd2   = " -DETECT_MINAREA "  +miniarea+ " -DETECT_THRESH 1.0 -ANALYSIS_THRESH 1.5"
    sexComd3   = " -SATUR_LEVEL %s -CHECKIMAGE_TYPE NONE"
    sexComd4   = " -PSF_NAME %s"
    sexComd    = sexComd1 + sexComd2 + sexComd3
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
