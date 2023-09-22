from astropy import stats
from astropy.io import fits
import numpy as np
import glob
import os
import pandas as pd

from loguru import logger
from lib.phot.ybias import save_fit

def combine_dark(dfilelist,expt,nsigma):
    aimga = []
    aimgb = []
    nfile=len(dfilelist)
    for i in range(0,nfile):
        img = fits.open(dfilelist[i])

        r,c=np.shape(np.array(img[0].data))
        aimga.append(img[0].data[:,0:int(c/2)])
        aimgb.append(img[0].data[:,int(c/2):int(c)])
                #aimgb.append(fits.open(filelist[i][0].data))[4088:8175,:]

    mdarka = stats.sigma_clipped_stats(aimga,axis=0,sigma=nsigma)[1]
    mdarkb = stats.sigma_clipped_stats(aimgb,axis=0,sigma=nsigma)[1]

    return mdarka,mdarkb
# from numba import jit
# @jit #(nopython=True)

def masterdark_old(ipath,opath,tid,bins,nsigma=3):
     
     dfilelist1=glob.glob(ipath+tid+'_'+'*dark*'+'.fit*')
     dfilelist=glob.glob(ipath+tid+'_'+'*DARK*'+'.fit*')
     dfilelist=dfilelist.extend(dfilelist1)
     #logger.info('*******************************')
     logger.info(dfilelist)
     nfile = len(dfilelist)
     if(nfile==1):
         return 0
     expt_dict={}
     expt_set={}
     groups = {}
        
     for i in range(0,nfile):
        hdr=fits.getheader(dfilelist[i])
        expt=int(hdr["EXPTIME"])
        expt_dict[expt]=dfilelist[i]
             
     logger.info(expt_dict)
     # for (key, value) in expt_dict:
     #    groups.setdefault(key, []).append(value)
     for key, value in expt_dict.items():
        groups.setdefault(key, []).append(value)
        logger.info('*******************************')
        logger.info(value)
        
      
     for key in groups:
         
        files=groups[key]
        if(len(files)>1):
            mdarka,mdarkb=combine_dark(files,key,nsigma)
            hdr=fits.getheader(files[0])
            hdr['EXPTIME']=key
            try:
                  save_fit(hdr,opath,tid+'_'+'master_dark_'+str(key)+'_'+bins+'_o1.fits',mdarka)
                  save_fit(hdr,opath,tid+'_'+'master_dark_'+str(key)+'_'+bins+'_o2.fits',mdarkb)
            except Exception as e:
                  logger.info('can not generate ydark')

            logger.info('*****************the dark '+tid+'_'+'master_dark_'+str(key)+'_'+bins+'_o1.fits  is ok')

def masterdark(ipath,opath,dfilelist,tid,bins,nsigma=3):
     
    #  dfilelist1=glob.glob(ipath+tid+'_'+'*dark*'+'.fit*')
    #  dfilelist=glob.glob(ipath+tid+'_'+'*DARK*'+'.fit*')
    #  dfilelist=dfilelist.extend(dfilelist1)
    #  #logger.info('*******************************')
     #logger.info(dfilelist)
     nfile = len(dfilelist)
     if(nfile==1):
         return 0
     expt_dict={}
     expt_set={}
     groups = {}
     dlistbin =[]
     for i in range(0,nfile):
        hdr=fits.getheader(dfilelist[i])
        expt=int(hdr["EXPTIME"])
        expt_dict[expt]=dfilelist[i]
        dlistbin.append(list([dfilelist[i], expt]))
     dlistbin=np.array(dlistbin)
     darkdf=pd.DataFrame({'filename':list(dlistbin[:,0]),'exp':list(dlistbin[:,1])})
     
     darkgroup = darkdf.groupby(['exp'])
     #logger.info(expt_dict)

     for key,values in darkgroup:
        files=list(values['filename'])
        # logger.info('#########***************#######****************')
        # logger.info(key)
        if(len(files)>1):
            mdarka,mdarkb=combine_dark(files,key,nsigma)
            hdr=fits.getheader(files[0])
            hdr['EXPTIME']=key
            try:
                  save_fit(hdr,opath,tid+'_'+'master_dark_'+bins+'_'+str(key)+'_o1.fits',mdarka)
                  save_fit(hdr,opath,tid+'_'+'master_dark_'+bins+'_'+str(key)+'_o2.fits',mdarkb)
            except Exception as e:
                  logger.info('can not generate ydark')

            logger.info('*****************the dark '+tid+'_'+'master_dark_'+bins+'_'+str(key)+'_o1.fits  is ok')

    
     # for (key, value) in expt_dict:
     #    groups.setdefault(key, []).append(value)
    #  for key, value in expt_dict.items():
    #     groups.setdefault(key, []).append(value)
        
    #     # logger.info('*******************************')
    #     # logger.info(value)
        
    #  for key in groups:
         
    #     files=groups[key]
    #     logger.info('#########***************#######****************')
    #     logger.info(key)
    #     if(len(files)>1):
    #         mdarka,mdarkb=combine_dark(files,key,nsigma)
    #         hdr=fits.getheader(files[0])
    #         hdr['EXPTIME']=key
    #         try:
    #               save_fit(hdr,opath,tid+'_'+'master_dark_'+bins+'_'+str(key)+'_o1.fits',mdarka)
    #               save_fit(hdr,opath,tid+'_'+'master_dark_'+bins+'_'+str(key)+'_o2.fits',mdarkb)
    #         except Exception as e:
    #               logger.info('can not generate ydark')

    #         logger.info('*****************the dark '+tid+'_'+'master_dark_'+bins+'_'+str(key)+'_o1.fits  is ok')
 
def masterbias(ipath,opath,bfilelist,tid,bins,nsigma=3):
        #:logger.info('tid=',tid)
        #bfilelist = glob.glob(ipath+tid+'_bias_'+'*.fit')

        nfile = len(bfilelist)
        aimga = []
        aimgb = []

        for i in range(0,nfile):
                img = fits.open(bfilelist[i])
                r,c=np.shape(np.array(img[0].data))
                
                #logger.info('r,c=',r,',',c)
                #aimga.append(img[0].data[:,500:int(c/2)])
                #aimgb.append(img[0].data[:,int(c/2):int(c)-500])
                aimga.append(img[0].data[:,0:int(c/2)])
                aimgb.append(img[0].data[:,int(c/2):int(c)])
                #aimgb.append(fits.open(filelist[i][0].data))[4088:8175,:]

        mbiasa = stats.sigma_clipped_stats(aimga,axis=0,sigma=nsigma)[1]
        mbiasb = stats.sigma_clipped_stats(aimgb,axis=0,sigma=nsigma)[1]

        hdr = fits.Header()
        #outpath = opath+date+'/'
        save_fit(hdr,opath,tid+'_'+'master_bias_'+bins+'_o1.fits',mbiasa)
        save_fit(hdr,opath,tid+'_'+'master_bias_'+bins+'_o2.fits',mbiasb)
        logger.info('the bias is ok')










 



    