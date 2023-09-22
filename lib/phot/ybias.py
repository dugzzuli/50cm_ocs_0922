from astropy import stats
from astropy.io import fits
import numpy as np
 

import glob
import os

from loguru import logger
from config import config
def get_current_dir():
    filepath=config["data"]["rootpath"]
    return filepath


def creat_hdu(from_header):
	empty_primary = fits.PrimaryHDU(header=from_header.copy())
	hduf = fits.HDUList([empty_primary])
	return hduf

@logger.catch
def save_fit(from_header,outpath,oname,img):
        try:
                 
                files = os.path.exists(outpath+oname)	
                # if files:
                #         os.remove(outpath+oname)

                hdu_list = creat_hdu(from_header)
                hdu = fits.ImageHDU(img)
                hdu_list.append(hdu)
                hdu_list.writeto(outpath+oname,overwrite=True)
        except:
                logger.info('wrong fits save')
        

def masterbias_nobins(ipath,opath,tid,nsigma=3):
	logger.info('tid={}'.format(tid))
	filelist = glob.glob(ipath+tid+'_bias_'+'*.fit')
	nfile = len(filelist)
	aimga = []
	aimgb = []

	for i in range(0,nfile):
		img = fits.open(filelist[i])
		aimga.append(img[0].data[:,0:2044])
		aimgb.append(img[0].data[:,2044:4088])
		#aimgb.append(fits.open(filelist[i][0].data))[4088:8175,:]

	mbiasa = stats.sigma_clipped_stats(aimga,axis=0,sigma=nsigma)[1]
	mbiasb = stats.sigma_clipped_stats(aimgb,axis=0,sigma=nsigma)[1]

	hdr = fits.Header()
	#outpath = opath+date+'/'
	save_fit(hdr,opath,tid+'_'+'master_bias'+'_o1.fits',mbiasa)
	save_fit(hdr,opath,tid+'_'+'master_bias'+'_o2.fits',mbiasb)
	logger.info('the bias is ok')

def masterbias_nocut(ipath,opath,bfilelist,tid,bins,nsigma=3):
        #:logger.info('tid=',tid)
        #bfilelist = glob.glob(ipath+tid+'_bias_'+'*.fit')
        
        nfile = len(bfilelist)
        aimga = []
        aimgb = []
        
        for i in range(0,nfile):
                img = fits.open(bfilelist[i])
                r,c=np.shape(np.array(img[0].data))
                #logger.info('r,c=',r,',',c)
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

def masterbias(ipath,opath,bfilelist,tid,bins,nsigma=3):
        #:logger.info('tid=',tid)
        #bfilelist = glob.glob(ipath+tid+'_bias_'+'*.fit')

        nfile = len(bfilelist)
        aimga = []
        aimgb = []
        if(nfile>20):
            nfile=20
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




