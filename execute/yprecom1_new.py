import time
import traceback

import ccdproc
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
import time as ts
import re
import sys
import os
import glob
from execute.flatdiff import mkdir 
import lib.phot.ygn as ygn
import lib.phot.ycalcom1 as ygn
from numpy.core.overrides import array_function_from_dispatcher 
from photutils.segmentation import make_source_mask
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground,Background2D 
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
import gc
import lib.phot.ybias as yb#creat_folder, single_fits_pack, creat_hdu_header,gen_overscan
import lib.phot.ylflat  as yf
import lib.phot.ydark as yd
import logging
from loguru import logger
from astropy.io.fits import Header


import time
import traceback

import ccdproc
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
import time as ts
import re
import sys
import os
import glob
from execute.flatdiff import mkdir
from lib.gchelp import elapsed_time 
import lib.phot.ygn as ygn
import lib.phot.ycalcom1 as ygn
from numpy.core.overrides import array_function_from_dispatcher 
from photutils.segmentation import make_source_mask
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground,Background2D 
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
import gc
import lib.phot.ybias as yb #creat_folder, single_fits_pack, creat_hdu_header,gen_overscan
import lib.phot.ylflat  as yf
import lib.phot.ydark as yd
import logging
from loguru import logger




def flatclassify(soucedatadir):
    fflist = glob.glob(soucedatadir+'*flat*.fit*')
    typeset=set()  
    for i in range(0,len(fflist)):
        filepath,tempfilename = os.path.split(fflist[i])
        tid,objectname,filterid,file_id=ygn.splitfilename(tempfilename)
        #loggerpre.info('tid,objectname,filterid,file_id',tid,objectname,filterid,file_id)
        typeset.add((tid,filterid))
    return list(typeset)


@logger.catch
def backgroud(img):
    logger.info('background begin')
    mask1= make_source_mask(img, nsigma=2, npixels=5, dilate_size=8)
    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = SExtractorBackground()
    bkg = Background2D(img, (50, 50), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,mask=mask1)
    #bkg = Background2D(img, (128, 128), filter_size=(5, 5),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,mask=mask1)
    img_bkg= bkg.background
    img_bkg_rms=bkg.background_rms
    logger.info('background end')
    return img_bkg,img_bkg_rms   


@logger.catch
def calibration1(ipath,opath):
    mkdir(opath)    
    flist= glob.glob(ipath+'fileguide/'+'*flat*.npy')
    tidlist=['y50a','y50b']
    filterlist=['mu','mv','mg','mr','mi','mz']
    binslist=['bin1','bin2']
    
    for bins in binslist:
        bbins=glob.glob(ipath+'fileguide/'+'*bias_'+str(bins)+'.npy')
        fbins=glob.glob(ipath+'fileguide/'+'*flat_'+str(bins)+'.npy')
        if(len(bbins)==0 and len(fbins)==0):
            continue

        for item in tidlist:
            tid=item
            if not os.path.exists(opath+tid+'_master_bias_'+bins+'_o1.fits'):
                logger.info('there is no mbias for '+tid+'  this day, begin to generate')
                try:
                   bfiles=np.load(ipath+'fileguide/'+tid+'_bias_'+bins+'.npy')
                   yb.masterbias(ipath,opath,bfiles,tid,bins,nsigma=3)
                except Exception as e:
                   logger.error(str(traceback.format_exc()))
                   logger.error('the bias genpro is failed ')
        for item in flist:
            tid,objname,filterid,bins=(ygn.get_filename(item)[1]).split("_")

            if not( os.path.exists(opath+tid+'_'+filterid+'_'+bins+'_gn.npy')):
                logger.info('there is no gn for this day, begin to generate gn')
                try:
                    bfiles=np.load(ipath+'fileguide/'+tid+'_bias_'+bins+'.npy',allow_pickle=True)
                except Exception as e:
                    logger.error(str(traceback.format_exc()))
                    bfiles=np.load(ipath+'fileguide/'+tid+'_Bias_'+bins+'.npy',allow_pickle=True)
                ffiles=np.load(ipath+'fileguide/'+tid+'_flat_'+filterid+'_'+bins+'.npy',allow_pickle=True)
                try:
                   ygn.check_gn(ipath,opath)
                except Exception as e:
                   logger.error(str(traceback.format_exc()))
                   logger.error('the gn genpro is failed ')
            if not os.path.exists(opath+tid+'_master_bias_'+bins+'_o1.fits'):
                logger.info('there is no mbias for '+tid+'  this day, begin to generate')
                try:
                   bfiles=np.load(ipath+'fileguide/'+tid+'_bias_'+bins+'.npy',allow_pickle=True)
                   yb.masterbias(ipath,opath,bfiles,tid,bins,nsigma=3)
                except Exception as e:
                   logger.error(str(traceback.format_exc()))
                   logger.error('the bias genpro is failed ')
            if not os.path.exists(opath+tid+'_master_dark_*'+bins+'*_o1.fits'):
               logger.info('there is no mdark for '+tid+'  this day, begin to generate')
               try:
                  dfiles=np.load(ipath+'fileguide/'+tid+'_dark_'+bins+'.npy',allow_pickle=True)
                  yd.masterdark(ipath,opath,tid,bins,nsigma=3)
               except Exception as e:
                  logger.error(str(traceback.format_exc()))
                  logger.error('the dark genpro is failed ')
            if not os.path.exists(opath+tid+'_master_flat_'+filterid+'_'+bins+'_o1.fits'):
                logger.info('there is no mflat for '+tid+'+'+filterid+'+'+bins+'  this day and this filter, begin to generate')
                try:
                   ffiles=np.load(ipath+'fileguide/'+tid+'_flat_'+filterid+'_'+bins+'.npy',allow_pickle=True)
                #yf.master_lflat(ipath,opath,ffiles,tid,filterid,bins)
                   yf.master_flat(ipath,opath,ffiles,tid,filterid,bins)
                except Exception as e:
                   logger.error(traceback.format_exc())
                   logger.error('the flat genpro is failed ')
    logger.info('calibration files is prepared')


    
def ttfname_header(filename):
    hdr=fits.getheader(filename)
    #OBJECT  = 'y50b_HZ44_' 
    tt=hdr['OBJECT']
    try: 
        tid = hdr['TELEID'].strip().lower()
        target = hdr['OBJECT'].strip()
        filterid = 'm'+hdr['FILTER'].strip()
    except:
        tid,target = tt.split('_')[0],tt.split('_')[1]
        filterid='m'+hdr['FILTER']
    return tid,target,filterid


@elapsed_time
def combine_process(filename,ipath,calpath,opath,libpath,date):
    mkdir(opath)
    
    filepath=ipath+filename
    filename1=filename.split('.fit')[0]
    # if(os.path.exists(opath+ filename1+'_subbkg_o0.fits')):
    #     fileguidelist=[]
    #     fileguidelist.append(str(opath+str(filename1)+'_sciimg_o1.fits'))
    #     fileguidelist.append(str(opath+str(filename1)+'_sciimg_o2.fits'))
    #     fileguidelist.append(str(opath+str(filename1)+'_mskimg_o1.fits'))
    #     fileguidelist.append(str(opath+str(filename1)+'_mskimg_o2.fits'))
    #     fileguidelist.append(str(opath+str(filename1)+'_flgimg_o1.fits'))
    #     fileguidelist.append(str(opath+str(filename1)+'_flgimg_o2.fits'))
    # # fileguidelist.append(str(opath+str(filename1)+'_uncimg_o1.fits'))
    # # fileguidelist.append(str(opath+str(filename1)+'_uncimg_o2.fits'))
    #     fileguidelist.append(str(opath+str(filename1)+'_bkg_o1.fits'))
    #     fileguidelist.append(str(opath+str(filename1)+'_bkg_o2.fits'))

    #     fileguidelist.append(str(opath+str(filename1)+'_subbkg_o0.fits'))
    #     return fileguidelist
    
    try:
        tid,objectid,filterid=ttfname_header(filepath)
        logger.info('++++++++++++++++++++++++++++++++++++++++')
        logger.info("{},{},{}".format(tid,objectid,filterid))
     
        logger.info('++++++++++++++++++++++++++++++++++++++++')
    except Exception as e:
        logger.error(str(traceback.format_exc()))

    
    #logger.info(filepath)
    img=fits.open(filepath)[0].data
    r,c=np.shape(img)
    icountmax=np.nanmax(img)
    #data = CCDData(img, unit=u.adu)
    img1=img[:,0:int(c/2)]
    img2=img[:,int(c/2):int(c)]
    #logger.info(filepath)
    #=================================================================fits头读取和参数增加与设置
    hdr_ccd=fits.getheader(filepath)
    binfact=int(hdr_ccd['XBINNING'])
    exptime = float(hdr_ccd['EXPTIME'])
    #pierside=hdr_ccd['PIERSIDE']
    #logger.info(binfact)
    pscale  = 0.297 * binfact # pixel scale in unit of arcsec
    yearmonth=str(date)[:-2]
    #logger.info(calpath+tid+'_mg_bin'+str(binfact)+'_gn.npy')
    try:
        if os.path.exists(calpath+tid+'_mg_bin'+str(binfact)+'_gn.npy'):
            #logger.info(calpath+tid+'_mg_bin'+str(binfact)+'_gn.npy')
            gn=np.load(calpath+tid+'_mg_bin'+str(binfact)+'_gn.npy',allow_pickle=True)
            logger.info('**********************',len(gn))
            gain1 =gn[0][0]
            noise1=gn[0][1]
            gain2 =gn[1][0]
            noise2=gn[1][1]
        elif os.path.exists(calpath+tid+'_mr_bin'+str(binfact)+'_gn.npy'):
            gn=np.load(calpath+tid+'_mr_bin'+str(binfact)+'_gn.npy',allow_pickle=True) 
            gain1 =gn[0][0]
            noise1=gn[0][1]
            gain2 =gn[1][0]
            noise2=gn[1][1]
        elif os.path.exists(calpath+tid+'_mi_bin'+str(binfact)+'_gn.npy'):
            gn =np.load(calpath+tid+'_mi_bin'+str(binfact)+'_gn.npy',allow_pickle=True) 
            gain1 =gn[0][0]
            noise1=gn[0][1]
            gain2 =gn[1][0]
            noise2=gn[1][1]
        elif os.path.exists(libpath):
            print('there is no gn this day,change to libpath##########################3')
            #print(libpath+tid+'_mg_bin'+str(binfact)+'_gn_'+ yearmonth+'*.npy')
            # cgnrlist=glob.glob(libpath+tid+'_mg_bin'+str(binfact)+'_gn_'+ yearmonth+'*.npy')
            # cgnglist=glob.glob(libpath+tid+'_mr_bin'+str(binfact)+'_gn_'+ yearmonth+'*.npy')
            cgnrlist=glob.glob(libpath+tid+'_mg_bin'+str(binfact)+'_gn_*.npy')
            cgnglist=glob.glob(libpath+tid+'_mr_bin'+str(binfact)+'_gn_*.npy')
            print(cgnrlist)
            if(len(cgnrlist)>0):
                #logger.info('##########################',cgnrlist)
                cgnrname=sorted(cgnrlist,reverse=True)[0]
                gain1,noise1=np.load(cgnrname,allow_pickle=True)[0][0],np.load(cgnrname,allow_pickle=True)[0][1]
                gain2,noise2=np.load(cgnrname,allow_pickle=True)[1][0],np.load(cgnrname,allow_pickle=True)[1][1]
                 
            elif(len(cgnglist)>0):
                cgngname=sorted(cgnglist,reverse=True)[0]
                gain1,noise1=np.load(cgngname,allow_pickle=True)[0][0],np.load(cgngname,allow_pickle=True)[0][1]
                gain2,noise2=np.load(cgngname,allow_pickle=True)[1][0],np.load(cgngname,allow_pickle=True)[1][1]
        
        hdr_subbkg=hdr_ccd.copy()
    
        #Header.add_history()
        hdr_subbkg["GAIN1"]=gain1 
        hdr_subbkg["GAIN2"]=gain2 
        hdr_subbkg['GAIN']=((gain1+gain2)/2)
        hdr_subbkg['BUNIT']='ADU'
        hdr_subbkg["RNOISE1"]=noise1
        hdr_subbkg["RNOISE2"]=noise2
        hdr_subbkg.add_history('preprocessing v20230215',before='GAIN1')##fy20230215
        
        hdr1=hdr_ccd.copy()
        hdr2=hdr_ccd.copy()
        
        hdr1["GAIN"]= gain1
        hdr1["RNOISE"]= noise1
        hdr1.add_history('preprocessing v20230215',before='GAIN')##fy20230215
         
        hdr2["GAIN"]= gain2
        hdr2["RNOISE"]= noise2
        hdr2.add_history('preprocessing v20230215',before='GAIN')##fy20230215
        logger.info('$$$$$$$$$$$$$$calibration is ok $$$$$$$$$$$$$$$$$$$$')
        
    except Exception as e:
               logger.error(str(traceback.format_exc()))
               logger.error('the gain/noise load may be error!')

     
    
    try:
        if(os.path.exists(calpath+tid+'_master_bias_bin'+str(binfact)+'_o1.fits')):
             master_bias1=fits.open(calpath+tid+'_master_bias_bin'+str(binfact)+'_o1.fits')[1].data
             master_bias2=fits.open(calpath+tid+'_master_bias_bin'+str(binfact)+'_o2.fits')[1].data
        elif (os.path.exists(libpath)):
             #logger.info(tid+'_master_bias_bin'+str(binfact)+'_o1_'+yearmonth+'*.fits')
             #cblist=glob.glob(libpath+tid+'_master_bias_bin'+str(binfact)+'_o1_'+yearmonth+'*.fits')
             cblist=glob.glob(libpath+tid+'_master_bias_bin'+str(binfact)+'_o1_*.fits')

             if(len(cblist)>0):
                cbname=sorted(cblist,reverse=True)[0]
                
                master_bias1=fits.open(cbname)[1].data
                master_bias2=fits.open(cbname[:-16]+'o2'+cbname[-14:])[1].data
                logger.info('@@@@@@@@@@@@@ the bias calibration will choose in lib: '+ cbname + ' @@@@@@@@@@')
        
        if(os.path.exists(calpath+tid+'_master_flat_'+filterid+'_bin'+str(binfact)+'_o2.fits')):
            #  master_flat1=fits.open(calpath+tid+'_master_flat_'+filterid+'_bin'+str(binfact)+'_o1.fits')[1].data
            #  master_flat2=fits.open(calpath+tid+'_master_flat_'+filterid+'_bin'+str(binfact)+'_o2.fits')[1].data
            master_flat1=fits.getdata(calpath+tid+'_master_flat_'+filterid+'_bin'+str(binfact)+'_o1.fits') 
            master_flat2=fits.getdata(calpath+tid+'_master_flat_'+filterid+'_bin'+str(binfact)+'_o2.fits') 
        else:
             
             #cfflist=glob.glob(libpath+tid+'_master_flat_'+filterid+'_bin'+str(binfact)+'_o1_'+yearmonth+'*.fits')
             cfflist=glob.glob(libpath+tid+'_master_flat_'+filterid+'_bin'+str(binfact)+'_o1_*.fits')
            
             if(len(cfflist)>0):
                cffname=sorted(cfflist,reverse=True)[0]
                print('#############',cffname)

                #master_flat1=fits.open(cffname)[1].data
                master_flat1=fits.getdata(cffname)
                cffname2=cffname.replace('_o1','_o2')
                print('#############',cffname2)
                #master_flat2=fits.open(cffname[:-16]+'o2'+cffname[-14:])[1].data
                #master_flat2=fits.open(cffname2)[1].data
                master_flat2=fits.getdata(cffname2)  
                logger.info('@@@@@@@@@@@@@ the flat calibration will choose in lib: '+ cffname + ' @@@@@@@@@@')
    except Exception as e:
               logger.error(str(traceback.format_exc()))
               logger.error('the calibration(master flat/bias) files load may be error!')


    logger.info('gain and noise is {},{},{}'.format(gain1,noise1,gain2,noise2))
    
    try:
        
        flat_corrected1=(img1-master_bias1)/master_flat1
        flat_corrected2=(img2-master_bias2)/master_flat2
    except Exception as e:
            logger.error(str(traceback.format_exc()))
            logger.error('jump in exp loop1')
            try:
                fr,fc=np.shape(master_flat1)
                br,bc=np.shape(master_bias1)
                ir,ic=np.shape(img1)
                min_r=np.min([fr,br,ir])
                min_c=np.min([fc,bc,ic])

                figc=8176/binfact
                figr=6132/binfact
                bincut=int(1000/binfact)
                if(fc==bc):
                    if(fc<ic and fc==bc):
                        img1=img[:,bincut:int(c/2)]
                        img2=img[:,0:int(c)-bincut]
                    elif(ic<fc and ic<bc):
                        master_flat1=master_flat1[:,bincut:int(fc)]
                        master_flat2=master_flat2[:,0:int(fc)-bincut]
                        master_bias1=master_bias1[:,bincut:int(fc)]
                        master_bias2=master_bias2[:,0:int(fc)-bincut]
                    elif(ic<fc and ic==bc):
                        master_flat1=master_flat1[:,bincut:int(fc)]
                        master_flat2=master_flat2[:,0:int(fc)-bincut]
                    elif(ic<bc and fc==bc):
                        master_bias1=master_bias1[:,bincut:int(fc)]
                        master_bias2=master_bias2[:,0:int(fc)-bincut]
   
                flat_corrected1=(img1-master_bias1)/master_flat1
                flat_corrected2=(img2-master_bias2)/master_flat2
                 
            except Exception as e:
                   logger.error(str(traceback.format_exc()))
                   logger.error('the preprocessing in (master flat/bias) files may be error! please check the size of image!')
    # initialize mask and flagmark
    # bad pixel
    mask11=np.ones_like(img1)
    mask11[np.where((master_flat1>1.5) | (master_flat1<0.7))]=0
    flag11=np.zeros_like(img1)
    flag11[np.where((master_flat1>1.5) | (master_flat1<0.7))]=1
    
    #hot pixel mask  
    mask12=np.ones_like(mask11)
    mask12[np.where(img1>60000)]=0 
    mask12=mask11*mask12
    flag12=np.zeros_like(flag11)
    flag12[np.where(img1>60000)]=2
    flag12=flag11+flag12
 

    #检测cosmicray
    #fdimg = CCDData(flat_corrected1, unit=u.adu)
    logger.info('cosmicray begin')
    cr_cleaned1, crmask1 = ccdproc.cosmicray_lacosmic(flat_corrected1, sigclip=5, gain=gain1, readnoise=noise1) 
    
 
    #cr_cleaned, crmask 
    mask13=np.array(~crmask1, dtype=int)*1
    flag13=np.zeros_like(flag11)
    flag13[np.where(mask13==0)]=4
    mask13=mask12*mask13
    flag13=flag12+flag13
    
 
    #检测satellite
    # mask4=np.ones_like(mask1)
    # flag4=np.zeros_like(flag1) 
    # if(len(sa.Satellite_line(np.array(cr_cleaned)))>1):
    #     mask4=sa.Satellite_line(np.array(cr_cleaned))
    #     flag4[np.where(mask4<1)]=8 
    # else:         
    #     logger.info('未检测到satellite')
    # mask4=mask3*mask4
    # flag4=flag3+flag4
    
    mask21=np.ones_like(img2)
    mask21[np.where((master_flat2>1.5) | (master_flat2<0.7))]=0
    flag21=np.zeros_like(img2)
    flag21[np.where((master_flat2>1.5) | (master_flat2<0.7))]=1
    
    #hot pixel mask  
    mask22=np.ones_like(mask21)
    mask22[np.where(img2>60000)]=0 
    mask22=mask21*mask22
    flag22=np.zeros_like(flag21)
    flag22[np.where(img2>60000)]=2
    flag22=flag21+flag22
 
    
    #检测cosmicray
    #fdimg = CCDData(flat_corrected1, unit=u.adu)
    cr_cleaned2, crmask2 = ccdproc.cosmicray_lacosmic(flat_corrected2, sigclip=5, gain=gain2, readnoise=noise2) 
    #logger.info('the mask type is:',type(crmask1[1][1]))
    logger.info('coemicray end')

    #cr_cleaned, crmask
 
    mask23=np.array(~crmask2, dtype=int)*1
    flag23=np.zeros_like(flag21)
    flag23[np.where(mask23==0)]=4
    mask23=mask22*mask23
    flag23=flag22+flag23
    
    flat_corrected1[crmask1]=np.nan
    flat_corrected2[crmask2]=np.nan
    kernel = Gaussian2DKernel(x_stddev=2, y_stddev=2, x_size=7, y_size=7)
    flat_corrected1 = interpolate_replace_nans(flat_corrected1,kernel)
    flat_corrected2 = interpolate_replace_nans(flat_corrected2,kernel)
    #gain corrected
    #cr_cleaned=CCDData(cr_cleaned1, unit=u.adu)
       
    flat_corrected1=flat_corrected1*gain1
    flat_corrected2=flat_corrected2*gain2
    gain_mean=(gain1+gain2)/2
    flat_corrected1=flat_corrected1/gain_mean
    flat_corrected2=flat_corrected2/gain_mean
    
    ####################################change##########################
 
    flat_corrected1_s=flat_corrected1
    
    flat_corrected2_s=flat_corrected2
    bkg1,bkgrms1=backgroud(flat_corrected1_s)
    bkg_sub1= flat_corrected1_s-bkg1
 
    bkg2,bkgrms2=backgroud(flat_corrected2_s)
    bkg_sub2= flat_corrected2_s-bkg2
    

    r1=int(r/2-1600/binfact)
    r2=int(r/2+1600/binfact)
 
    c1=int(int(c/2)-50-1200/binfact)
    c2=int(c/2)-50

    # gate b:
    c3= 50  
    c4=int(50+1200/binfact) 

    center_bkg1=bkg1[r1:r2,c1:c2]
    center_bkg2=bkg2[r1:r2,c3:c4]
    bkg_ratio=round(np.nanmedian(center_bkg1)/np.nanmedian(center_bkg2),3)
 
    bkg_sub2_norm=bkg_sub2*bkg_ratio
    #bkg_ratio2=round(np.nanmedian(center_bkg2)/np.nanmedian(center_bkg1),3)
    bkg_sub1_norm=bkg_sub1#*bkg_ratio2
    

    logger.info('bkg_ratio={}'.format(bkg_ratio))
    logger.info('np.median(flat_corrected1)={}, np.median(bkg1)={}, np.median(subbkg1)={},subbkg1_norm={}'.format(
        np.nanmedian(flat_corrected1), np.nanmedian(bkg1), np.nanmedian(bkg_sub1), np.nanmedian(bkg_sub1_norm),
        np.nanmedian(bkg_sub1_norm)))

    logger.info(
        ('np.median(flat_corrected2)={}, np.median(bkg2)={}, np.median(subbkg2)={},subbkg2_norm={}').format(
            np.nanmedian(flat_corrected2), np.nanmedian(bkg2), np.nanmedian(bkg_sub2), np.nanmedian(bkg_sub2_norm)))

    #gain_corrected1=bkg_sub1* gain1
    #gain_corrected2=bkg_sub2* gain2
    bkg_sub=np.hstack((bkg_sub1_norm,bkg_sub2_norm))
    

    #写入fits头
    saturate1=(65535.0-sigma_clipped_stats(master_bias1)[1])-np.nanmedian(bkg1)
    saturate2=(65535.0-sigma_clipped_stats(master_bias2)[1])-np.nanmedian(bkg2)
    saturate1=round(saturate1,2)
    saturate2=round(saturate2,2) 
    
    hdr_subbkg["SATURAT1"] = saturate1
    hdr_subbkg["SATURAT2"] = saturate2
    hdr_subbkg["SATURATE"] = np.min([saturate1,saturate2])*0.9
    hdr_subbkg.add_history('background eastimation, v20230215',before='SATURATE')##fy20230215
    hdr_subbkg["SKYBKG1"] = round(np.nanmedian(bkg1),2)
    hdr_subbkg["SKYBKG2"] = round(np.nanmedian(bkg2),2)
    hdr_subbkg["SKYRMS1"] =round(sigma_clipped_stats(bkgrms1)[1],2)
    hdr_subbkg["SKYRMS2"] =round(sigma_clipped_stats(bkgrms2)[1],2)
    # hdr_subbkg["SUBBKG1"] = round(np.nanmedian(bkg_sub1),2)
    # hdr_subbkg["SUBBKG2"] = round(np.nanmedian(bkg_sub2),2)
    # hdr_subbkg["CSUBBKG2"] = round(np.nanmedian(bkg_sub2_norm),2)
    hdr_subbkg["SUBBKG"] = round(np.nanmedian(bkg_sub),2)
    #hdr_subbkg["BKG1D2"]=bkg_ratio
     
    hdr1["SATURATE"]=saturate1
    hdr1.add_history('background eastimation',before='SATURATE')
    hdr1["SKYBKG"]= hdr_subbkg["SKYBKG1"]
    hdr1["SKYRMS"]=hdr_subbkg["SKYRMS1"]
 
    hdr2["SATURATE"]=saturate2
    hdr2.add_history('background eastimation',before='SATURATE')
    hdr2["SKYBKG"]=hdr_subbkg["SKYBKG2"]
    hdr2["SKYRMS"]=hdr_subbkg["SKYRMS2"]

   #save_fit(hdr_ccd,opath,filename1+'_fsciimg_o'+str(ccdindex)+'.fits',flat_corrected)
    yb.save_fit(hdr1,opath,str(filename1)+'_sciimg_o1.fits',flat_corrected1)
    yb.save_fit(hdr1,opath,str(filename1)+'_mskimg_o1.fits',mask13)
    yb.save_fit(hdr1,opath,str(filename1)+'_flgimg_o1.fits',flag13)    #flag-img
    #save_fit(hdr_ccd,opath,str(filename1)+'_uncimg_o1.fits',np.array(data_with_deviation1))     
    yb.save_fit(hdr1,opath,str(filename1)+'_bkg_o1.fits',np.array(bkg1))    #img-bkg
    yb.save_fit(hdr1,opath,str(filename1)+'_cosmsk_o1.fits',np.array(mask13))

    yb.save_fit(hdr2,opath,str(filename1)+'_cosmsk_o2.fits',np.array(mask23))
    yb.save_fit(hdr2,opath,str(filename1)+'_sciimg_o2.fits',flat_corrected2)
    yb.save_fit(hdr2,opath,str(filename1)+'_mskimg_o2.fits',mask23)
    yb.save_fit(hdr2,opath,str(filename1)+'_flgimg_o2.fits',flag23)    #flag-img
    # yb.#save_fit(hdr_ccd,opath,str(filename1)+'_uncimg_o2.fits',np.array(data_with_deviation2))  
    yb.save_fit(hdr2,opath,str(filename1)+'_bkg_o2.fits',np.array(bkg2))    #img-bkg
    

    yb.save_fit(hdr_subbkg,opath,str(filename1)+'_subbkg_o0.fits',np.array(bkg_sub))    #img-bkg


    fileguidelist=[]
    fileguidelist.append(str(opath+str(filename1)+'_sciimg_o1.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_sciimg_o2.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_mskimg_o1.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_mskimg_o2.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_flgimg_o1.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_flgimg_o2.fits'))
   # fileguidelist.append(str(opath+str(filename1)+'_uncimg_o1.fits'))
   # fileguidelist.append(str(opath+str(filename1)+'_uncimg_o2.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_bkg_o1.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_bkg_o2.fits'))

    fileguidelist.append(str(opath+str(filename1)+'_subbkg_o0.fits'))
    
    #释放局部变量，但是不hi到会不会有异常
    logger.info(locals().keys())
    for x in list(locals().keys()):
        try:
            del locals()[x]
        except:
            continue
    gc.collect()
    
    return fileguidelist

@logger.catch
@elapsed_time
def pro_combine(tid_target_filter,date):

    t0=ts.time()
    rootdir=yb.get_current_dir()
    
    #libpath=rootdir+'lib/'
    if(int(date)<20211001):
       libpath=rootdir+'lib/cal1/'
    else: 
       libpath=rootdir+'lib/cal3/'

    rawpath=rootdir+'reception/'+str(date)+'/raw/'
    fileguide_raw_path=rawpath+'fileguide/'
    calpath=rootdir+'reception/'+str(date)+'/cal/'
    scipath=rootdir+'reception/'+str(date)+'/sci/'+tid_target_filter+'/'
    mkdir(scipath)
    fileguide_sci_path=scipath+'fileguide/'
    mkdir(fileguide_sci_path)
    logdir = rootdir+'run_pipe_log/'+str(date)+'/'+tid_target_filter+'/'
  
    mkdir(logdir)
    calibration1(rawpath,calpath)
    # set logging
    skey = str(date) + '_' + tid_target_filter + '_prep'
    loggn = logdir + skey + ".log"
    if os.path.exists(loggn): os.system("rm %s" % loggn)
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logFile = logging.FileHandler(loggn)
    logFile.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger(skey).addHandler(logFile)
    loggerpre = logging.getLogger(skey)

    loggerpre.info("^_^ Telescope_filter_target: %s" % (tid_target_filter))

    loggerpre.info("[_yprep] "+"^_^ Telescope_filter_target: %s"%(tid_target_filter))

    t1=ts.time()
    tcal=t1-t0
    loggerpre.info(f"[_yprep] "+f'the calibration is done,the time cost is:{str(tcal)}')
    loggerpre.info(f"[_yprep] "+f'the calibration is done,the time cost is:{str(tcal)}')
    loggerpre.info(f"[_yprep] "+f"the calibration is done,the time cost is:{str(tcal)}")
    loggerpre.info(f"[_yprep] "+f"the calibration is done,the time cost is: is:{str(tcal)}")
   
    if not (os.path.exists(fileguide_raw_path+'sci_'+tid_target_filter+'.npy')):#fy20230217
        loggerpre.info("There is no such observation: %s" % (tid_target_filter)) 
        return None
    #filelist= glob.glob(rawinpath+tid+'_'+target+'_'+filterid+'_'+'*.fit')
    filelist=np.load(fileguide_raw_path+'sci_'+tid_target_filter+'.npy') #fy20230217
    nfile=len(filelist)#fy20230217
     
    #record the fileguidelist for different type 
    fgl_sciimg1=[]
    fgl_mskimg1=[]
    fgl_flgimg1=[]
    fgl_uncimg1=[]
    fgl_bkg1=[]
   
    fgl_sciimg2=[]
    fgl_mskimg2=[]
    fgl_flgimg2=[]
    fgl_uncimg2=[]
    fgl_bkg2=[]
    
    fgl_subbkg=[] 
    
    for i in range(0,nfile):
        dirname,filename=os.path.split(filelist[i])
        try:
            dirname,filename=os.path.split(filelist[i])
            listc=combine_process(filename,dirname+'/',calpath,scipath,libpath,date)
            #list2=combine_process(filename,dirname+'/',calpath,scipath,2)
            fgl_sciimg1.append(listc[0])
            fgl_sciimg2.append(listc[1])

            fgl_mskimg1.append(listc[2])
            fgl_mskimg2.append(listc[3])

            fgl_flgimg1.append(listc[4])
            fgl_flgimg2.append(listc[5])

           # fgl_uncimg1.append(listc[6])
           # fgl_uncimg2.append(listc[7])

            fgl_bkg1.append(listc[6])
            fgl_bkg2.append(listc[7])

            fgl_subbkg.append(listc[8])
            loggerpre.info("[_yprep] "+'****************** '+filename+' preprocessing is ok!')
            #filename=str(listc[8]).rsplit('/',1)[1]
            loggerpre.info("[_yprep] "+"The preprocessing is ok!: %s" %(filename))
            loggerpre.info("[_yprep] "+"The preprocessing is ok!: %s" %(filename))
        except Exception as e:
            loggerpre.error(str(traceback.format_exc()))
            #filename=str(listc[8]).rsplit('/',1)[1]
            loggerpre.error("XXXXXXXXXXXXXXXXX The preprocessing is failed: %s" %(filename))
            loggerpre.info("XXXXXXXXXXXXXXXXX The preprocessing is failed: %s" %(filename))
        finally:
            for k in list(locals().keys()):
            # if locals[k] is np.nan:
                try:
                    del locals[k]
                except:
                    continue
            gc.collect()

    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_sciimg_o1.npy',fgl_sciimg1)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_mskimg_o1.npy',fgl_mskimg1)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_mskimg_o1.npy',fgl_mskimg1)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_mskimg_o1.npy',fgl_mskimg1)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_mskimg_o1.npy',fgl_mskimg1)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_mskimg_o1.npy',fgl_mskimg1)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_mskimg_o1.npy',fgl_mskimg1)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_flgimg_o1.npy',fgl_flgimg1)
    #savenpy(fileguide_sci_path,tid_target_filter+'_uncimg_o1.npy',fgl_uncimg1)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_bkg_o1.npy',fgl_bkg1)
    
    
    

    
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_sciimg_o2.npy',fgl_sciimg2)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_mskimg_o2.npy',fgl_mskimg2)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_flgimg_o2.npy',fgl_flgimg2)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_bkg_o2.npy',fgl_bkg2)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_subbkg_o0.npy',fgl_subbkg)
    
    # try:
    #     del fgl_sciimg1
    #     del fgl_mskimg1
    #     del fgl_mskimg1
    #     del fgl_mskimg1
    #     del fgl_mskimg1
    #     del fgl_mskimg1
    #     del fgl_mskimg1
    #     del fgl_flgimg1
    #     del fgl_bkg1
    #     del fgl_sciimg2
    #     del fgl_mskimg2
    #     del fgl_flgimg2
    #     del fgl_bkg2
    #     del fgl_subbkg
    # except:
    #     pass
    
    for k in list(locals().keys()):
        # if locals[k] is np.nan:
        try:
            del locals[k]
        except:
            continue
        
    gc.collect()
    
    
    t2=ts.time()
    tpro=t2-t1
    loggerpre.info("[_yprep] "+'the preprcess of'+ tid_target_filter+ 'is done, the time cost is:{}'.format(str(tpro)))