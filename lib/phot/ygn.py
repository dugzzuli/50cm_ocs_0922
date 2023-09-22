import math
import numpy as np
from astropy import stats
from astropy.table import Table
from astropy.stats import sigma_clipped_stats,SigmaClip
import os
from astropy.io import fits
import glob
import scipy.signal as signal

from loguru import logger

from lib.phot.ycalcom1 import savenpy


def qname(q):
	qname=''
	if (int(q)<10):
		qname='0'+str(q)
	else:
		qname=str(q)
	return qname
        

def GainNoise(Flat1,Flat2,Zero1,Zero2):
    #Gain = (Flatplus - Zeroplus) / (Var(Flatdif ) - Var(Zerodif))
    #Noise=Gain*std(eZero)/math.sqrt(2)
    r, c = np.shape(Flat1)
    Flatdif = Flat1 +(-1)*Flat2
    Zerodif = Zero1+ (-1)*Zero2
    Flatplus= Flat1 + Flat2
    Zeroplus = Zero1 + Zero2
    FZ=(sigma_clipped_stats(Flat1)[1]+sigma_clipped_stats(Flat2)[1])-(sigma_clipped_stats(Zero1)[1]+sigma_clipped_stats(Zero2)[1])
    sigma_Flatdif=stats.sigma_clipped_stats(Flatdif)[2]#std
    sigma_Zerodif=stats.sigma_clipped_stats(Zerodif)[2]#std
    gain=FZ/(sigma_Flatdif**2-sigma_Zerodif**2)
    noise= (gain * sigma_Zerodif) / math.sqrt(2)
    return gain, noise



def Noise(Zero1,Zero2,gain):
    r, c = np.shape(Zero1)
    
    Zerodif = Zero1+ (-1)*Zero2
    
    Zeroplus = Zero1 + Zero2
     
    
    sigma_Zerodif=stats.sigma_clipped_stats(Zerodif)[2]#std
     
    noise= (gain * sigma_Zerodif) / math.sqrt(2)
    return noise



def quality(flatimg):
        if(np.median(flatimg)<8000 or np.median(flatimg)>40000):
                return False
        return True


def checkflat(f1,f2):
    diff=np.median(f1)-np.median(f2)
    if (abs(diff)>4000):
        return True
    return False


#calculate gain and noise batch for one sub-ccd
#input: zeros_file_directory and flat_file_directory


def noise_gen(ipath,opath,Zerolist,tid,bins,gaina,gainb):
    imgb=fits.open(Zerolist[0])
    
    rf,cf=np.shape(np.array(imgb[0].data))
     
 
    noiselist=[]
    lennum= len(Zerolist)
    binning=int(bins[-1])
  
    zimga=[]
    zimgb=[]
    #gate a:
    r1=int(rf/2-1600/binning)
    r2=int(rf/2+1600/binning)
 
    c1=int(int(cf/2)-50-1200/binning)
    c2=int(cf/2)-50

    # gate b:
    c3= 50+int(cf/2)
    c4=int(50+1200/binning)+int(cf/2)
    #fileda=imga[r1:r2,c1:c2]
    #filedb=imgb[r1:r2,40:440]
    logger.info(r1,':',r2,',',c1,':',c2)
    for i in range(0,lennum):
         
        zimg = fits.open(Zerolist[i])

        
        zimga.append(zimg[0].data[r1:r2,c1:c2])
        zimgb.append(zimg[0].data[r1:r2,c3:c4])

     
    a_noiselist=[]
     
    b_noiselist=[]
    for i in range(0, lennum-1):
         
        na = Noise(zimga[i], zimga[i+1],gaina)
        logger.info('a gate:the noise=',na)
         
        a_noiselist.append(na)

        nb = Noise(zimgb[i], zimgb[i+1],gainb)
        logger.info('b gate:the noise=',nb)
         
        b_noiselist.append(nb)

    if(len(a_noiselist)>2 and len(b_noiselist)>2):
        
       a_noise=stats.sigma_clipped_stats(a_noiselist)[1]
        
       
       b_noise=stats.sigma_clipped_stats(b_noiselist)[1]
       gn=[[a_noise, b_noise],a_noiselist,b_noiselist]
       logger.info(gn)
    #hdr = fits.Header()

    # save_fit(hdr,opath,tid+'_'+filterid+'_'+'gn'+'_o01.fits',nflata)
    # save_fit(hdr,opath,tid+'_'+filterid+'_'+'gn'+'_o02.fits',nflatb)
       savenpy(opath,tid+'_'+bins+'_noise.npy',gn)

       logger.info(' noise of '+ tid+' : noise_a=',a_noise,'noise_b=',b_noise)
    else:
       logger.info('gain noise of can not generate for not enough qualify calibrations')


#ygn.gn_gen(ipath,opath,bfiles,ffiles,filterid,tid,bins)

def gn_gen(ipath,opath,Zerolist,Flatlist,filterid,tid,bins):
    #notice here may need try catch to ensure these files exists
    #Zerolist= glob.glob(ipath+'/'+tid+'*bias*.fit')
    #Flatlist= glob.glob(ipath+tid+'_Tflat_'+filterid+'*.fit')#'y50a_Tflat_jr*.fit'
    imgb=fits.open(Zerolist[0])
    imgf=fits.open(Flatlist[0])
    rb,cb=np.shape(np.array(imgb[0].data))
    rf,cf=np.shape(np.array(imgf[0].data))
    if(cb!=cf):
       logger.info('the shape of bias and flat files is not same!')
       return 0
    #Flatlist= glob.glob(ipath+date+'/'+tid+'*flat*.fit')
    gainlist=[]
    noiselist=[]
    lennum=min(len(Flatlist),len(Zerolist))
    binning=int(bins[-1])
    fimga=[]
    fimgb=[]
    zimga=[]
    zimgb=[]
    #gate a:
    r1=int(rf/2-1600/binning)
    r2=int(rf/2+1600/binning)
 
    c1=int(int(cf/2)-50-1200/binning)
    c2=int(cf/2)-50

    # gate b:
    c3= 50+int(cf/2)
    c4=int(50+1200/binning)+int(cf/2)
    #fileda=imga[r1:r2,c1:c2]
    #filedb=imgb[r1:r2,40:440]
    logger.info(r1,':',r2,',',c1,':',c2)
    for i in range(0,lennum):
        fimg = fits.open(Flatlist[i])
        zimg = fits.open(Zerolist[i])

        fimga.append(fimg[0].data[r1:r2,c1:c2])
        fimgb.append(fimg[0].data[r1:r2,c3:c4])
        zimga.append(zimg[0].data[r1:r2,c1:c2])
        zimgb.append(zimg[0].data[r1:r2,c3:c4])

    a_gainlist=[]
    a_noiselist=[]
    b_gainlist=[]
    b_noiselist=[]
    for i in range(0, lennum-1):
        if(quality(fimga[i]) or quality(fimgb[i])):
           logger.info('the flat is not qualify in this loop:{}'.format(i))
        if (checkflat(fimga[i], fimga[i+1]) or checkflat(fimgb[i], fimgb[i+1])):
            logger.info('this loop is out'+str(i)+'and'+str(i+1))
            continue

        ga, na = GainNoise(fimga[i], fimga[i+1], zimga[i], zimga[i+1])
        logger.info('a gate:the gain={}, noise= {}'.format(ga,na))
        a_gainlist.append(ga)
        a_noiselist.append(na)

        gb, nb = GainNoise(fimgb[i], fimgb[i+1], zimgb[i], zimgb[i+1])
        logger.info('b gate:the gain= {}, noise={}'.format(gb, nb))
        b_gainlist.append(gb)
        b_noiselist.append(nb)

    if(len(a_gainlist)>2 and len(b_gainlist)>2):
       a_gain=stats.sigma_clipped_stats(a_gainlist)[1]
       a_noise=stats.sigma_clipped_stats(a_noiselist)[1]
       ta=Table()
       b_gain=stats.sigma_clipped_stats(b_gainlist)[1]
       b_noise=stats.sigma_clipped_stats(b_noiselist)[1]
       gn=[[a_gain,a_noise],[b_gain,b_noise]]
    #hdr = fits.Header()

    # save_fit(hdr,opath,tid+'_'+filterid+'_'+'gn'+'_o01.fits',nflata)
    # save_fit(hdr,opath,tid+'_'+filterid+'_'+'gn'+'_o02.fits',nflatb)
       savenpy(opath,tid+'_'+filterid+'_'+bins+'_gn.npy',gn)

       logger.info('gain noise of '+ tid+' + '+filterid+' : gain_a=',a_gain,'noise_a=',a_noise,'gain_b=',b_gain,'noise_b=',b_noise)
    else:
       logger.info('gain noise of can not generate for not enough qualify calibrations')




