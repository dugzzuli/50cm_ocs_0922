from astropy import stats
from astropy.io import fits
import numpy as np
import glob
import os

from scipy import ndimage, misc
import math
import time as ts
from astropy.stats import SigmaClip,sigma_clipped_stats

from loguru import logger
from lib.phot.ybias import creat_hdu, save_fit


def weight_iter(fimglist,wlist,masklist):
	nfile=len(fimglist)
	wimg=[]
	fmean=[]
	fstd=[]
	for i in  range(0,nfile):
		wimg.append(wlist[i]*masklist[i])#
	wimg_sum=np.array(wimg).sum(axis=0)
	for i in range (0,nfile):
		fmean.append(fimglist[i]*wlist[i]*masklist[i]/wimg_sum)
	fmean_sum=np.array(fmean).sum(axis=0)
	for i in range (0,nfile):
		fstd.append(wlist[i]*masklist[i]/wimg_sum*np.square(fimglist[i]-fmean_sum))	
	fstd_sum=np.sqrt(np.array(fstd,dtype=float).sum(axis=0))
	return fmean_sum,fstd_sum


def iter_mask(fnimg,fmean_sum,fstd_sum,nfile,sigma=2.5):
	masklist=[]
	for i in range(0,nfile):
		mask=np.ones_like(fmean_sum)
		diff=abs(fnimg[i]-fmean_sum)-sigma*np.array(fstd_sum,dtype=float)
		mask[np.where(diff>0)]=0
		masklist.append(mask)
	return masklist

def median_smooth(img,size):
	data = ndimage.median_filter(np.array(img), size)
	return data

def sigma_iter(fnimglist,wimglist,sigma=2.5,times=3):
	nfile=len(fnimglist)
	fmean_sum,fstd_sum=weight_iter(fnimglist,wimglist,np.ones(nfile))
	for i in range(0,times-1):
		maskalist_new=iter_mask(fnimglist,fmean_sum,fstd_sum,nfile,sigma)
		fmean_sum,fstd_sum=weight_iter(fnimglist,wimglist,maskalist_new)
	return fmean_sum,fstd_sum
	 

def quality(flatimg):
	if(np.median(flatimg)<15000 or np.median(flatimg)>35000):
		return False
	return True

def master_lflat_nobins(ipath,opath,tid,filterid):
	
	mbiasa=fits.open(opath+tid+'_'+'master_bias_o1.fits')[1].data
	mbiasb=fits.open(opath+tid+'_'+'master_bias_o2.fits')[1].data
	rawfilelist = glob.glob(ipath+tid+'_Tflat_'+filterid+'*.fit*')#'y50a_Tflat_jr*.fit'
	logger.info('the len of flat file is:',len(rawfilelist))
	filelist=[]
	nfile = len(filelist)
	fimga = []
	fimgb = []
	fnimga=[]
	fnimgb=[]
	# fmimga=[]
	# fmimgb=[]

    #check quality of flat
	for i in range(0,len(rawfilelist)):
		img = fits.open(rawfilelist[i])
		if(quality(img[0].data)):
			filelist.append(img[0].data)
		else:

			logger.info('the flat'+str(filelist[i]+'is a bad flat'))
			logger.info(np.median(img[0].data)) 
	nfile = len(filelist)

	if(nfile<5):
		logger.info('there many bad flat and the master flat can not be built')
		return 0
	if(nfile>20):
		nfile=20
	for i in range(0,nfile):#len(rawfilelist)):
		#logger.info('===========================================')
		#logger.info('the raw flat is:',filelist[i] )
		#logger.info('the bias is','mbiasa_std=',np.std(mbiasa),'mbiasb_std=',np.std(mbiasb))
		r,c=np.shape(np.array(img[0].data))
		 
		imga=img[0].data[:,0:int(c/2)]+(-1)*mbiasa  
		imgb=img[0].data[:,int(c/2):c]+(-1)*mbiasb
		fimga.append(imga)
		fimgb.append(imgb)
 
		fileda=imga[1484:1583,int(c/2)-144:int(c/2)-44]
		filedb=imgb[1484:1583,44:144]
        
		fnimga.append(np.array(imga,dtype=float)/np.mean(fileda))
		fnimgb.append(np.array(imgb,dtype=float)/np.mean(filedb))
	 
	fmeana_sum,fstda_sum=sigma_iter(fnimga,fimga,sigma=2.5,times=3)
	fmeanb_sum,fstdb_sum=sigma_iter(fnimgb,fimgb,sigma=2.5,times=3)

	hdr = fits.Header()#
	 
	save_fit(hdr,opath,tid+'_master_lflat_'+filterid+'_o1.fits',fmeana_sum)
	save_fit(hdr,opath,tid+'_master_lflat_'+filterid+'_o2.fits',fmeanb_sum)
	logger.info('the master flat is ok')
 

def master_lflat(ipath,opath,flatfilelist,tid,filterid,bins):
        
    mbiasa=fits.open(opath+tid+'_'+'master_bias_'+bins+'_o1.fits')[1].data
    mbiasb=fits.open(opath+tid+'_'+'master_bias_'+bins+'_o2.fits')[1].data
    #rawfilelist = glob.glob(ipath+tid+'_Tflat_'+filterid+'*.fit')#'y50a_Tflat_jr*.fit'
    logger.info('the len of flat file is:',len(flatfilelist))
    filelist=[]
    nfile = len(filelist)
    fimga = []
    fimgb = []
    fnimga=[]
    fnimgb=[]
    # fmimga=[]
    # fmimgb=[]

    #check quality of flat
    for i in range(0,len(flatfilelist)):
        img = fits.open(flatfilelist[i])
        if(quality(img[0].data)):
                filelist.append(np.array(img[0].data))
        else:
                logger.info('the flat'+str(filelist[i]+'is a bad flat'))
                logger.info(np.median(img[0].data)) 
    nfile = len(filelist)

    if(nfile<5):
        logger.info('there many bad flat and the master flat can not be built')
        return 0
    if(nfile>20):
        nfile=20
    for i in range(0,nfile):
        #logger.info('===========================================')
        #logger.info('the raw flat is:',filelist[i] )
        #logger.info('the bias is','mbiasa_std=',np.std(mbiasa),'mbiasb_std=',np.std(mbiasb))
        r,c=np.shape(np.array(filelist[i].data))
        binning=int(bins[-1])
        imga=filelist[i][:,0:int(c/2)]+(-1)*mbiasa
        imgb=filelist[i][:,int(c/2):int(c)]+(-1)*mbiasb
        fimga.append(imga)
        fimgb.append(imgb)

        #fileda=imga[1484:1583,int(c/2)-144:int(c/2)-44]
        #filedb=imgb[1484:1583,44:144]
        fileda=imga[int(r/2)-int(100/binning):int(r/2)+int(100/binning),int(c/2)-int(200/binning)-50:int(c/2)-50]
        filedb=imgb[int(r/2)-int(100/binning):int(r/2)+int(100/binning),50:int(200/binning)+50]

        fnimga.append(np.array(imga,dtype=float)/np.mean(fileda))
        fnimgb.append(np.array(imgb,dtype=float)/np.mean(filedb))

    fmeana_sum,fstda_sum=sigma_iter(fnimga,fimga,sigma=2.5,times=3)
    fmeanb_sum,fstdb_sum=sigma_iter(fnimgb,fimgb,sigma=2.5,times=3)

    hdr = fits.Header()#
    logger.info('the normalized filed gatea:['+str(int(r/2)-(100/binning))+':'+str(int(r/2)+(100/binning))+','+str(int(c/2)-(200/binning+100/binning))+':'+str(int(c/2)-(100/binning))+']')
    save_fit(hdr,opath,tid+'_master_lflat_'+filterid+'_'+bins+'_o1.fits',fmeana_sum)
    save_fit(hdr,opath,tid+'_master_lflat_'+filterid+'_'+bins+'_o2.fits',fmeanb_sum)
    logger.info('the master flat is ok')



def master_flat_nocut(ipath,opath,flatfilelist,tid,filterid,bins):

        mbiasa=fits.open(opath+tid+'_'+'master_bias_'+bins+'_o1.fits')[1].data
        mbiasb=fits.open(opath+tid+'_'+'master_bias_'+bins+'_o2.fits')[1].data
        #rawfilelist = glob.glob(ipath+tid+'_Tflat_'+filterid+'*.fit')#'y50a_Tflat_jr*.fit'
        logger.info('the len of flat file is:',len(flatfilelist))
        filelist=[]
        nfile = len(filelist)
        fimga = []
        fimgb = []
        fnimga=[]
        fnimgb=[]
        # fmimga=[]
        # fmimgb=[]

    #check quality of flat
        for i in range(0,len(flatfilelist)):
                img = fits.open(flatfilelist[i])
                if(quality(img[0].data)):
                        filelist.append(np.array(img[0].data))
                else:
                        logger.info('the flat'+str(flatfilelist[i]+'is a bad flat'))
        nfile = len(filelist)

        if(nfile<5):
            logger.info('there many bad flat and the master flat can not be built')
            return 0
        for i in range(0,len(filelist)):
            #logger.info('===========================================')
            #logger.info('the raw flat is:',filelist[i] )
            #logger.info('the bias is','mbiasa_std=',np.std(mbiasa),'mbiasb_std=',np.std(mbiasb))
            r,c=np.shape(np.array(filelist[i].data))
            binning=int(bins[-1])
            imga=filelist[i][:,0:int(c/2)]+(-1)*mbiasa
            imgb=filelist[i][:,int(c/2):int(c)]+(-1)*mbiasb
            fimga.append(imga)
            fimgb.append(imgb)

            #fileda=imga[1484:1583,int(c/2)-144:int(c/2)-44]
            #filedb=imgb[1484:1583,44:144]
            fileda=imga[int(r/2)-int(200/binning):int(r/2)+int(200/binning),int(c/2)-int(200/binning)-50:int(c/2)-50]
            filedb=imgb[int(r/2)-int(200/binning):int(r/2)+int(200/binning),50:int(200/binning)+50]

            fnimga.append(np.array(imga,dtype=float)/sigma_clipped_stats(fileda)[0])
            fnimgb.append(np.array(imgb,dtype=float)/sigma_clipped_stats(filedb)[0])

        mflata = stats.sigma_clipped_stats(fnimga,axis=0,sigma=2.5)[1]
        mflatb = stats.sigma_clipped_stats(fnimgb,axis=0,sigma=2.5)[1]

        hdr = fits.Header()#
        #logger.info('the normalized filed gatea:['+str(int(r/2)-(100/binning))+':'+str(int(r/2)+(100/binning))+','+str(int(c/2)-(200/binning+100/binning))+':'+str(int(c/2)-(100/binning))+']')
        save_fit(hdr,opath,tid+'_master_flat_'+filterid+'_'+bins+'_o1.fits',mflata)
        save_fit(hdr,opath,tid+'_master_flat_'+filterid+'_'+bins+'_o2.fits',mflatb)
        logger.info('the master flat is ok') 





def master_flat(ipath,opath,flatfilelist,tid,filterid,bins):

        mbiasa=fits.open(opath+tid+'_'+'master_bias_'+bins+'_o1.fits')[1].data
        mbiasb=fits.open(opath+tid+'_'+'master_bias_'+bins+'_o2.fits')[1].data
        #rawfilelist = glob.glob(ipath+tid+'_Tflat_'+filterid+'*.fit')#'y50a_Tflat_jr*.fit'
        logger.info(f'the len of flat file is:{len(flatfilelist)}')
        filelist=[]
        nfile = len(filelist)
        fimga = []
        fimgb = []
        fnimga=[]
        fnimgb=[]
        # fmimga=[]
        # fmimgb=[]

    #check quality of flat
        for i in range(0,len(flatfilelist)):
            try:
                img = fits.open(flatfilelist[i])
                if(quality(img[0].data)):
                    filelist.append(np.array(img[0].data))
                else:
                    logger.info('the flat'+str(flatfilelist[i]+'is a bad flat'))
                    logger.info(np.median(img[0].data)) 
            except:
                pass
        nfile = len(filelist)

        if(nfile<5):
                logger.info('there many bad flat and the master flat can not be built')
                return 0
        if(nfile>20):
            nfile=20
        for i in range(0, nfile):
            #logger.info('===========================================')
            #logger.info('the raw flat is:',filelist[i] )
            #logger.info('the bias is','mbiasa_std=',np.std(mbiasa),'mbiasb_std=',np.std(mbiasb))
            r,c=np.shape(np.array(filelist[i].data))
            binning=int(bins[-1])
            #imga=filelist[i][:,500:int(c/2)]-mbiasa
            #imgb=filelist[i][:,int(c/2):int(c)-500]+-mbiasb
            imga=filelist[i][:,0:int(c/2)]+(-1)*mbiasb
            imgb=filelist[i][:,int(c/2):int(c)]+(-1)*mbiasb
            fimga.append(imga)
            fimgb.append(imgb)

            #fileda=imga[1484:1583,int(c/2)-144:int(c/2)-44]
            #filedb=imgb[1484:1583,44:144]
            #fileda=imga[int(r/2)-int(200/binning):int(r/2)+int(200/binning),int(c/2)-500-int(200/binning)-50:int(c/2)-500-50]
                
            #filedb=imgb[int(r/2)-int(200/binning):int(r/2)+int(200/binning),50:int(200/binning)+50]
            fileda=imga[int(r/2)-int(200/binning):int(r/2)+int(200/binning),int(c/2)-int(200/binning)-50:int(c/2)-50]
            filedb=imgb[int(r/2)-int(200/binning):int(r/2)+int(200/binning),50:int(200/binning)+50]
            fnimga.append(np.array(imga,dtype=float)/sigma_clipped_stats(fileda)[0])
            fnimgb.append(np.array(imgb,dtype=float)/sigma_clipped_stats(filedb)[0])

        mflata = stats.sigma_clipped_stats(fnimga,axis=0,sigma=2.5)[1]
        mflatb = stats.sigma_clipped_stats(fnimgb,axis=0,sigma=2.5)[1]

        hdr = fits.Header()#
        #logger.info('the normalized filed gatea:['+str(int(r/2)-(100/binning))+':'+str(int(r/2)+(100/binning))+','+str(int(c/2)-(200/binning+100/binning))+':'+str(int(c/2)-(100/binning))+']')
        save_fit(hdr,opath,tid+'_master_flat_'+filterid+'_'+bins+'_o1.fits',mflata)
        save_fit(hdr,opath,tid+'_master_flat_'+filterid+'_'+bins+'_o2.fits',mflatb)
        logger.info('the master flat is ok')


                        
