import math
import numpy as np
from astropy import stats
from astropy.table import Table
from astropy.stats import sigma_clipped_stats,SigmaClip
import os
from astropy.io import fits
import glob
import scipy.signal as signal
from execute.flatdiff import mkdir
import lib.phot.ybias as yb#creat_folder, single_fits_pack, creat_hdu_header,gen_overscan
import lib.phot.ylflat  as yf
import lib.phot.ydark as yd
import logging
import matplotlib


matplotlib.use('Agg')
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import pylab as pl
import os, sys
from loguru import logger



 




def savenpy(opath,filename,data):
    files = os.path.exists(opath+filename)  
    if files:
        os.remove(opath+filename)
    np.save(opath+filename,data)

def qname(q):
    qname=''
    if (int(q)<10):
        qname='0'+str(q)
    else:
        qname=str(q)
    return qname

def get_filename(fullfilename):
    fpath,fullname=os.path.split(fullfilename)
    fname,ext=os.path.splitext(fullname)
    return fpath,fname


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


# def splitfilename(filename):
# #    logger.info('====================================================',filename)
#     namelist=filename.split(".")[0].split('_')
#     #python存为字典
#     tid=namelist[0]
#     objectname=namelist[1]
#     file_id_filterid=namelist[2]
#     filterid=file_id_filterid[-2:]
#     file_id=file_id_filterid[1:5]
#     return tid,objectname,filterid,file_id

def splitfilename(filename):
    namelist=filename.split('_')
    #python存为字典
    tid=namelist[0]
    objectname=namelist[1]
    filterid=namelist[2]
    file_id=namelist[3]
    return tid,objectname,filterid,file_id

def quality(flatimg):
        if(np.median(flatimg)<16000 or np.median(flatimg)>35000):
                return False
        return True


def checkflat(f1,f2):
    diff=np.median(f1)-np.median(f2)
    if (abs(diff)>3000):
        return True
    return False


#calculate gain and noise batch for one sub-ccd
#input: zeros_file_directory and flat_file_directory



def diff_bias(ipath,opath,bfilelist,tid,bins,nsigma=3,gaina=1,gainb=1):

    #:logger.info('tid=',tid)
    figdir=opath+'figure/'
    bfilelist = glob.glob(ipath+tid+'_bias_'+'*.fit')
    nfile = len(bfilelist)
    aimga = []
    aimgb = []
    diffa_median=[]
    diffb_median=[]
    diff_std_a=[]
    diff_std_b=[]
    na=[]
    nb=[]
    nfile=10
    for i in range(0,nfile):
            img = fits.open(bfilelist[i])
            r,c=np.shape(np.array(img[0].data))
            imga=img[0].data[:,0:int(c/2)]
            imgb=img[0].data[:,int(c/2):int(c)]
            aimga.append(imga)
            aimgb.append(imgb)

            if (i>0):
                img0 = fits.open(bfilelist[i-1])
                logger.info(i)
                # diffimga=img0[0].data[:,0:int(c/2)]
                # diffimgb=img0[0].data[:,int(c/2):int(c)]
    
                # diff_a=sigma_clipped_stats(imga-diffimga)
                # diff_b=sigma_clipped_stats(imgb-diffimgb) 
                diff_a=sigma_clipped_stats(imga+(-1)*aimga[i-1])
                diff_b=sigma_clipped_stats(imgb+(-1)*aimgb[i-1]) 
                logger.info('diff_a=',diff_a)
                diff_std_a.append(diff_a[2])
                diff_std_b.append(diff_a[2])
                diffa_median.append(diff_a[1])
                diffb_median.append(diff_b[1])
                na.append((gaina * diff_a[2]) / math.sqrt(2) )
                nb.append((gainb * diff_b[2]) / math.sqrt(2) )

    amed    = np.median(na)
    bmed    = np.median(nb)
    logger.info('amed=',amed)
    logger.info('bmed={}'.format(bmed))
    listlen=np.linspace(1,nfile-1,nfile-1)
    xlim     = [0.5, nfile-0.5]
    pl.scatter(listlen, na, color="red", marker="o", s=6)
    pl.scatter(listlen, nb, color="blue", marker="*", s=6)
    pl.plot(xlim, [amed, amed],"r-",linewidth=1.0)
    pl.plot(xlim, [bmed, bmed],"b-",linewidth=1.0)
    # pl.plot(xlim, [zpMed-zpStd, zpMed-zpStd],"r--",linewidth=1.5)
    # pl.plot(xlim, [zpMed+zpStd, zpMed+zpStd],"r--",linewidth=1.5)
    pl.xlim(xlim)
    pl.xlabel("Bias_number", fontsize=10)
    pl.ylabel("Diff_Flux", fontsize=10)
    pl.title(tid+'_Diff_Bias',fontsize=10)
    pl.savefig(figdir+tid+'_diff_bias.pdf')
    pl.clf()
    pl.close()
     
    return na,nb,diff_std_a,diff_std_b,diffa_median,diffb_median
    logger.info('the bias diff check is ok')


#ygn.gn_gen(ipath,opath,bfiles,ffiles,filterid,tid,bins)


#def recommend_gain(gainlist):
    





def gn_gen(ipath,opath,Zerolist,Flatlist,filterid,tid,bins):
    #notice here may need try catch to ensure these files exists
    #Zerolist= glob.glob(ipath+'/'+tid+'*bias*.fit')
    #Flatlist= glob.glob(ipath+tid+'_Tflat_'+filterid+'*.fit')#'y50a_Tflat_jr*.fit'

    figdir=opath+'figure/'
    mkdir(figdir)
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
    r1=int(rf/2-1200/binning)
    r2=int(rf/2+1200/binning)
 
    c1=int(int(cf/2)-50-800/binning)
    c2=int(cf/2)-50

    # gate b:
    c3= 50+int(cf/2)
    c4=int(50+800/binning)+int(cf/2)
    #fileda=imga[r1:r2,c1:c2]
    #filedb=imgb[r1:r2,40:440]
    #logger.info(r1,':',r2,',',c1,':',c2)
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
            logger.info('the flat is not qualify in this loop:'+str(i))
        if (checkflat(fimga[i], fimga[i+1]) or checkflat(fimgb[i], fimgb[i+1])):
            logger.info('this loop is out',i,'and',i+1)
            continue

        ga, na = GainNoise(fimga[i], fimga[i+1], zimga[i], zimga[i+1])
        logger.info(f'a gate:the gain={ga} noise={na}')
        a_gainlist.append(ga)
        a_noiselist.append(na)

        gb, nb = GainNoise(fimgb[i], fimgb[i+1], zimgb[i], zimgb[i+1])
        logger.info(f'b gate:the gain={gb} noise={nb}')
        b_gainlist.append(gb)
        b_noiselist.append(nb)

    if(len(a_gainlist)>2 and len(b_gainlist)>2):
        a_gain=stats.sigma_clipped_stats(a_gainlist)[1]
        a_gain_std=stats.sigma_clipped_stats(a_gainlist)[2]
        a_noise=stats.sigma_clipped_stats(a_noiselist)[1]
        a_noise_std=stats.sigma_clipped_stats(a_noiselist)[2]
        ta=Table()
        b_gain=stats.sigma_clipped_stats(b_gainlist)[1]
        b_gain_std=stats.sigma_clipped_stats(b_gainlist)[2]
        b_noise=stats.sigma_clipped_stats(b_noiselist)[1]
        b_noise_std=stats.sigma_clipped_stats(b_noiselist)[2]
        gn=[[a_gain,a_noise,a_gain_std,a_noise_std],[b_gain,b_noise,b_gain_std,b_noise_std],a_gainlist,a_noiselist,b_gainlist,b_noiselist]

        plt.subplot(1,2,1)  #等效于plt.subplot(221)
        
        listlen=np.linspace(1,len(a_gainlist),len(a_gainlist))
        xlim     = [0.5, len(a_gainlist)+0.5]
        plt.scatter(listlen, a_gainlist, color="red", marker="o", s=6)
        #pl.scatter(listlen, b_gainlist, color="blue", marker="*", s=6)
        plt.plot(xlim, [a_gain, a_gain],"r-",linewidth=2.0)
        #pl.plot(xlim, [b_gain, b_gain],"b-",linewidth=2.0)
        plt.plot(xlim, [a_gain-a_gain_std, a_gain-a_gain_std],"r--",linewidth=1.5)
        plt.plot(xlim, [a_gain+a_gain_std, a_gain+a_gain_std],"r--",linewidth=1.5)
        plt.xlim(xlim)
        plt.xlabel("Number", fontsize=10)
        plt.ylabel("Gain", fontsize=10)
        plt.tick_params(labelsize=8)  
        plt.title(tid+'_left_gain='+str(round(a_gain,3))+'+/-'+str(round(a_gain_std,3)),fontsize=10)

        plt.subplot(1,2,2)
        listlen=np.linspace(1,len(b_gainlist),len(b_gainlist))
        xlim     = [0.5, len(b_gainlist)+0.5]
        plt.scatter(listlen, b_gainlist, color="red", marker="o", s=6)
        #pl.scatter(listlen, b_gainlist, color="blue", marker="*", s=6)
        plt.plot(xlim, [b_gain, b_gain],"r-",linewidth=2.0)
        #pl.plot(xlim, [b_gain, b_gain],"b-",linewidth=2.0)
        plt.plot(xlim, [b_gain-b_gain_std, b_gain-b_gain_std],"r--",linewidth=1.5)
        plt.plot(xlim, [b_gain+b_gain_std, b_gain+b_gain_std],"r--",linewidth=1.5)
        plt.xlim(xlim)   
        plt.xlabel("Number", fontsize=10)    
        plt.tick_params(labelsize=8)  
        #plt.ylabel("Gain", fontsize=5)
        plt.title(tid+'_right_gain='+str(round(b_gain,3))+'+/-'+str(round(b_gain_std,3)),fontsize=10)

        plt.savefig(figdir+tid+'_'+filterid+'_gain.pdf')
        #plt.clf()
        plt.close()

    # save_fit(hdr,opath,tid+'_'+filterid+'_'+'gn'+'_o01.fits',nflata)
    # save_fit(hdr,opath,tid+'_'+filterid+'_'+'gn'+'_o02.fits',nflatb)
        savenpy(opath,tid+'_'+filterid+'_'+bins+'_gn.npy',gn)
        logger.info(gn)

        logger.info('gain noise of '+ tid+' + '+filterid+' : gain_a='+str(a_gain)+'noise_a='+str(a_noise)+'gain_b='+str(b_gain)+'noise_b='+str(b_noise))
    else:   
        logger.info('gain noise of can not generate for not enough qualify calibrations')
    





def calibration(ipath,opath):
    mkdir(opath)    
    blist = glob.glob(ipath+'fileguide/'+'*bias*.npy')
    flist= glob.glob(ipath+'fileguide/'+'*flat*.npy')
    dlist= glob.glob(ipath+'fileguide/'+'*dark*.npy')
    logger.info('the calibration tid and filters=',flist)
    tidlist=['y50a','y50b']
    filterlist=['mu','mv','mg','mr','mi','mz']

    for item in tidlist:
        tid=item
        bins='bin2'
        if not os.path.exists(opath+tid+'_master_bias_'+bins+'_o1.fits'):
            logger.info('there is no mbias for '+tid+'  this day, begin to generate')
            try:
               bfiles=np.load(ipath+'fileguide/'+tid+'_bias_'+bins+'.npy')
               yb.masterbias(ipath,opath,bfiles,tid,bins,nsigma=3)
            except Exception as e:
               logger.info(e)
               logger.info('the bias genpro is failed ')
    for item in flist:
        tid,objname,filterid,bins=splitfilename(get_filename(item)[1])

        if not( os.path.exists(opath+tid+'_'+filterid+'_'+bins+'_gn.npy')):
            logger.info('there is no gn for this day, begin to generate gn')
            try:
                bfiles=np.load(ipath+'fileguide/'+tid+'_bias_'+bins+'.npy')
            except Exception as e:
                logger.info(e)
                bfiles=np.load(ipath+'fileguide/'+tid+'_Bias_'+bins+'.npy')
                ffiles=np.load(ipath+'fileguide/'+tid+'_flat_'+filterid+'_'+bins+'.npy')
            try:
               gn_gen(ipath,opath,bfiles,ffiles,filterid,tid,bins)
            except Exception as e:
               logger.info(e)
               logger.info('the gn genpro is failed ')

        if not os.path.exists(opath+tid+'_master_dark_'+bins+'_o1.fits'):
           logger.info('there is no mdark for '+tid+'  this day, begin to generate')
           try:
              dfiles=np.load(ipath+'fileguide/'+tid+'_dark_'+bins+'.npy')
              yd.masterdark(ipath,opath,tid,bins,nsigma=3)
           except Exception as e:
              logger.info(e)
              logger.info('the dark genpro is failed ')
        if not os.path.exists(opath+tid+'_master_flat_'+filterid+'_'+bins+'_o1.fits'):
            logger.info('there is no mflat for '+tid+'+'+filterid+'+'+bins+'  this day and this filter, begin to generate')
            try:
              ffiles=np.load(ipath+'fileguide/'+tid+'_flat_'+filterid+'_'+bins+'.npy')

              yf.master_flat(ipath,opath,ffiles,tid,filterid,bins)
            except Exception as e:
              logger.info(e)
              logger.info('the flat genpro is failed ')
    logger.info('calibration files is prepared')




def check_gn(ipath,opath):
    mkdir(opath)    
    flist= glob.glob(ipath+'fileguide/'+'*flat*.npy')
    logger.info('the calibration tid and filters=',flist)

    #logger.info(flist)

    
    for item in flist:
        logger.info(f'####{item}')
        tid,objname,filterid,bins=splitfilename(get_filename(item)[1])

        if not( os.path.exists(opath+tid+'_'+filterid+'_'+bins+'_gn.npy')):
            logger.info('there is no gn for this day, begin to generate gn')
        #if(1==1):
            try:
                bfiles=np.load(ipath+'fileguide/'+tid+'_bias_'+bins+'.npy')
            except Exception as e:
                logger.info(e)
                bfiles=np.load(ipath+'fileguide/'+tid+'_Bias_'+bins+'.npy')
            ffiles=np.load(ipath+'fileguide/'+tid+'_flat_'+filterid+'_'+bins+'.npy')
            try:
               gn_gen(ipath,opath,bfiles,ffiles,filterid,tid,bins)
            except Exception as e:
               logger.info(e)
               logger.info('the gn genpro is failed ')



def check_bias(ipath,opath):
    blist = glob.glob(ipath+'fileguide/'+'*bias*.npy')
    tidlist=['y50a','y50b']
    
    logger.info(blist[0])
    
    
    for item in blist:
         
        tid,objname,bins=item.rsplit('/',1)[1].rsplit('.')[0].split('_')
        logger.info(tid,objname,bins)
        
        try:
            bfiles=np.load(ipath+'fileguide/'+tid+'_bias_'+bins+'.npy')
        except Exception as e:
            logger.info(e)
            bfiles=np.load(ipath+'fileguide/'+tid+'_Bias_'+bins+'.npy')

        diff_bias(ipath,opath,bfiles,tid,bins,nsigma=3,gaina=1,gainb=1)


def pro(date):
      
    rootdir=yb.get_current_dir()
    rawpath=rootdir+'reception/'+str(date)+'/raw/'
    fileguide_raw_path=rawpath+'fileguide/'
    calpath=rootdir+'reception/'+str(date)+'/cal/'
    #logdir = rootdir+'run_pipe_log/'+str(date)+'/calibration/'
    
    check_gn(rawpath,calpath)

    check_bias(rawpath,calpath)


 



