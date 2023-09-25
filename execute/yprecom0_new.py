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

from execute.yphotutils import datename,d2hms
from execute.yphotutils import make_get_request,datename,d2hms,mad,wds9reg,read_param,crossmatch,reg,ref_expt_select,readfits,hdr_edit_fits,re_hdr_edit_fits,read_list,timeplus,HJDTrans,get_scamp_head,overlapRect,pointRect,_spherical_to_cartesian,_great_circle_distance,mkdir
from execute import yphotutils



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


            hdark = glob.glob(opath+tid+'_master_dark_'+bins+'*_o1.fits')
            if not (hdark):
            #if not os.path.exists(opath+tid+'_master_dark_'+bins+'*_o1.fits'):
                logger.info('there is no mdark for '+tid+'  this day, begin to generate')
                try:
                    dfiles=np.load(ipath+'fileguide/'+tid+'_dark_'+bins+'.npy',allow_pickle=True)
                    yd.masterdark(ipath,opath,dfiles,tid,bins,nsigma=3)
                    #yd.masterdark(ipath,opath,tid,bins,nsigma=3)
                except Exception as e:
                    logger.error(str(traceback.format_exc()))
                    logger.error('the dark genpro is failed ')

        for item in flist:
            tid,objname,filterid,bins=(ygn.get_filename(item)[1]).split("_")
            if not( os.path.exists(opath+tid+'_'+filterid+'_'+bins+'_gn.npy')):
                logger.info('there is no gn for this day, begin to generate gn')
                temp_npy=ipath+'fileguide/'+tid+'_bias_'+bins+'.npy'
                if os.path.exists(temp_npy):
                    bfiles=np.load(ipath+'fileguide/'+tid+'_bias_'+bins+'.npy',allow_pickle=True)
                
                ffiles=np.load(ipath+'fileguide/'+tid+'_flat_'+filterid+'_'+bins+'.npy',allow_pickle=True)
                try:
                    ygn.check_gn(ipath,opath)
                except Exception as e:
                    logger.error(f"the gn genpro is failed :{str(traceback.format_exc())}")


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

 
    



def combine_process(filename,ipath,calpath,opath,libpath,date,repro=0):
    mkdir(opath)
    
    filepath=ipath+filename
    filename1=filename.split('.fit')[0]
    if(filename1[-1]=='_'):
        filename1 = filename1[:-1]
    logger.info('filename1 is: '+filename1)

    if( repro==1 and os.path.exists(opath+ filename1+'_sciimg.fits') ):
        fileguidelist=[]
        fileguidelist.append(str(opath+str(filename1)+'_fcimg.fits'))
     
        fileguidelist.append(str(opath+str(filename1)+'_mskimg.fits'))
     
        fileguidelist.append(str(opath+str(filename1)+'_flgimg.fits'))
         
    # fileguidelist.append(str(opath+str(filename1)+'_uncimg_o1.fits'))
    # fileguidelist.append(str(opath+str(filename1)+'_uncimg_o2.fits'))
        fileguidelist.append(str(opath+str(filename1)+'_bkg.fits'))
     

        fileguidelist.append(str(opath+str(filename1)+'_sciimg.fits'))
        return fileguidelist
    
    try:
        tid,objectid,filterid=yphotutils.ttfname_header(filepath)
        logger.info('++++++++++++++++++++++++++++++++++++++++')
        logger.info("{},{},{}".format(tid,objectid,filterid))
     
        logger.info('++++++++++++++++++++++++++++++++++++++++')
    except Exception as e:
        logger.error(str(traceback.format_exc()))
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": filepath, "reason": "these are no parameters of object ra and dec in header of this fits: -------"+str(traceback.format_exc()),"stage": 0}
            response = make_get_request(url, params)
            return None 
        except:
            pass


        
    #logger.info(filepath)
    img=fits.open(filepath)[0].data
    
    r,c=np.shape(img)
    hdr_ccd=fits.getheader(filepath)
    binfact=int(hdr_ccd['XBINNING'])
     
    if(r<6132/binfact or c<8176/binfact):
        logger.error(f' 6132/binfact={6132/binfact} or  8176/binfact={8176/binfact}')
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": filepath, "reason": f"the file shape may be error! -------row={r},col={c}", "stage": 0}
            response = make_get_request(url, params)
            return None 
        except:
            pass
        
    img1=img[:,0:int(c/2)]
    img2=img[:,int(c/2):int(c)]
    #logger.info(filepath)
    logger.info(np.shape(img1))
    logger.info(np.shape(img2))
    #=================================================================fits头读取和参数增加与设置
    
    try:
        if os.path.exists(calpath+tid+'_mr_bin'+str(binfact)+'_gn.npy'):
            #logger.info(calpath+tid+'_mg_bin'+str(binfact)+'_gn.npy')
            gn=np.load(calpath+tid+'_mr_bin'+str(binfact)+'_gn.npy',allow_pickle=True)
            logger.info('**********************',len(gn))
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
            cgnrlist=glob.glob(libpath+tid+'_mr_bin'+str(binfact)+'_gn_*.npy')
            cgnglist=glob.glob(libpath+tid+'_mi_bin'+str(binfact)+'_gn_*.npy')
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
        
        hdr_sciimg=hdr_ccd.copy()
    
        #Header.add_history()
        hdr_sciimg["GAIN1"]=gain1 
        hdr_sciimg["GAIN2"]=gain2 
        hdr_sciimg['GAIN']=((gain1+gain2)/2)
        hdr_sciimg['BUNIT']='ADU'
        hdr_sciimg["RNOISE1"]=noise1
        hdr_sciimg["RNOISE2"]=noise2
        hdr_sciimg.add_history('preprocessing v20230215',before='GAIN1')##fy20230215
        
         
        logger.info('$$$$$$$$$$$$$$calibration is ok $$$$$$$$$$$$$$$$$$$$')
        
    except Exception as e:
        logger.error(str(traceback.format_exc()))
        logger.error('the gain/noise load may be error!')
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": filepath, "reason": "the gain/noise load may be error! -------"+str(traceback.format_exc()),"stage": 0}
            response = make_get_request(url, params)
            return None 
        except:
            pass

     
    
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
                
                
        logger.info('gain and noise is {},{},{}'.format(gain1,noise1,gain2,noise2))
        flat_corrected1=(img1-master_bias1)/master_flat1
        flat_corrected2=(img2-master_bias2)/master_flat2
        
    except Exception as e:
        
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": filepath, "reason": "the calibration(master flat/bias) files load may be error! -------"+str(traceback.format_exc()),"stage": 0}
            response = make_get_request(url, params)
            return None 
        except:
            pass
        raise(traceback.format_exc())


    


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
    flat_corrected_o0=np.hstack((flat_corrected1,flat_corrected2))
    
    mask_o0=np.hstack((np.array(mask13,dtype='int16'),np.array(mask23,dtype='int16')))
    flag_o0=np.hstack((np.array(flag13,dtype='int16'),np.array(flag23,dtype='int16')))
    if(c>8000 and binfact==1) or (c> 4000 and binfact ==2) :
        mask_o0= mask_o0[:,int(1000/binfact):c-int(1000/binfact)]
        flag_o0= flag_o0[:,int(1000/binfact):c-int(1000/binfact)]
        flat_corrected_o0 = flat_corrected_o0[:,int(1000/binfact):c-int(1000/binfact)]
        logger.info(f"#################imgcut################## :{np.shape(flat_corrected_o0)}")
           
         
    flat_corrected_o0=np.array(flat_corrected_o0,dtype='float32') 
    yb.save_fit(hdr_sciimg,opath,str(filename1)+'_mskimg.fits',mask_o0)
    yb.save_fit(hdr_sciimg,opath,str(filename1)+'_flgimg.fits',flag_o0)  
    #gain corrected
    #cr_cleaned=CCDData(cr_cleaned1, unit=u.adu)
       
    flat_corrected1=flat_corrected1*gain1
    flat_corrected2=flat_corrected2*gain2
    gain_mean=(gain1+gain2)/2
    gc1=flat_corrected1/gain_mean
    gc2=flat_corrected2/gain_mean
    
    
    bkg1,bkgrms1=backgroud(gc1)
    bkg2,bkgrms2=backgroud(gc2)
    
    bkg_o0=np.hstack((bkg1,bkg2))
    bkg_rms_o0=np.hstack((bkgrms1,bkgrms2))
 

    #写入fits头
    saturate1=(65535.0-sigma_clipped_stats(master_bias1)[1])-np.nanmedian(bkg1)
    saturate2=(65535.0-sigma_clipped_stats(master_bias2)[1])-np.nanmedian(bkg2)
    saturate1=round(saturate1,2)
    saturate2=round(saturate2,2)
    hdr_sciimg["SATURAT1"] = saturate1
    hdr_sciimg["SATURAT2"] = saturate2
    hdr_sciimg["SATURATE"] = np.min([saturate1,saturate2])
    hdr_sciimg["SKYBKG1"] = round(np.nanmedian(bkg1),2)
    hdr_sciimg["SKYBKG2"] = round(np.nanmedian(bkg2),2)
    hdr_sciimg["SKYRMS1"] =round(sigma_clipped_stats(bkgrms1)[1],2)
    hdr_sciimg["SKYRMS2"] =round(sigma_clipped_stats(bkgrms2)[1],2)
    bkg_o0=np.array(bkg_o0,dtype='float32')
    if(c>8000 and binfact==1) or (c> 4000 and binfact ==2) :
        bkg_o0= bkg_o0[:,int(1000/binfact):c-int(1000/binfact)]
    yb.save_fit(hdr_sciimg,opath,str(filename1)+'_bkg.fits',bkg_o0)

    binfact=int(hdr_sciimg['XBINNING'])
    pscale  = 0.297 * binfact # pixel scale in unit of arcsec
    ira='0'
    idec='0'
    try:
        if("OBJCTRA" in hdr_sciimg.keys()):
            ira  = ":".join(hdr_sciimg["OBJCTRA"].split())
            idec = ":".join(hdr_sciimg["OBJCTDEC"].split())
            logger.warning(hdr_sciimg["OBJCTDEC"])
            ira, idec = d2hms(ira, idec, conv=1)
        elif("RA" in hdr_sciimg.keys()):
            ira  = ":".join(hdr_sciimg["RA"].split())
            idec = ":".join(hdr_sciimg["DEC"].split())
            ira, idec = d2hms(ira, idec, conv=1)
 
    except:
        logger.info('these are no parameters of object ra and dec in header of this fits')
        logger.error(str(traceback.format_exc()))
          
        try:
            url = 'http://12.12.12.251:8888/process_request'
            params = {"eid": filepath, "reason": "these are no parameters of object ra and dec in header of this fits! -------"+str(traceback.format_exc()),"stage": 0}
            response = make_get_request(url, params)
            return None 
        except:
            pass
        
    
    hdr_sciimg['EPOCH']=2000.0
    hdr_sciimg["CTYPE1"]  = "RA---TAN"
    hdr_sciimg["CTYPE2"]  = "DEC--TAN"
    hdr_sciimg["EQUINOX"] = 2000.0
    hdr_sciimg["RADESYS"] = "ICRS"
    hdr_sciimg["CUNIT1"]  = "deg"
    hdr_sciimg["CUNIT2"]  = "deg"
    hdr_sciimg["CRVAL1"]  = ira
    hdr_sciimg["CRVAL2"]  = idec
    hdr_sciimg["CRPIX1"]  = c/2
    hdr_sciimg["CRPIX2"]  = r/2
    # #ihdr["SATURATE"] = np.min([65535.0,icountmax])

    # # determine the rotation matrix

    hdr_sciimg["CD1_1"]   = -pscale/3600.0
    hdr_sciimg["CD1_2"]   = 0.0
    hdr_sciimg["CD2_1"]   = 0.0
    hdr_sciimg["CD2_2"]   = -pscale/3600.0
    # hdr_sciimg.add_history('backgroud estimation',before='SATURAT1')##fy20230215

 
    wcsfile=ipath+filename[:-5]+'.wcs'
     
    if(os.path.exists(wcsfile)):
        try:
            logger.info(wcsfile)
            
            hdr_wcs = fits.getheader(wcsfile)
            logger.info(hdr_wcs["CD1_1"])
            logger.info(hdr_wcs["CD1_2"])
            logger.info(hdr_wcs["CD2_1"])
            logger.info(hdr_wcs["CD2_2"])

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
            
            

        except:
            logger.info('these are no wcs of in header of this fits')
            logger.error(str(traceback.format_exc()))
            
            
    
    yb.save_fit(hdr_sciimg,opath,str(filename1)+'_fcimg.fits',flat_corrected_o0)

    sciimg_o0=np.hstack((gc1,gc2)) 
    sciimg_o0=np.array(sciimg_o0,dtype='float32')
    if(c>8000 and binfact==1) or (c> 4000 and binfact ==2) :
        sciimg_o0= sciimg_o0[:,int(1000/binfact):c-int(1000/binfact)]
    yb.save_fit(hdr_sciimg,opath,str(filename1)+'_sciimg.fits',sciimg_o0)
    logger.info('sciimg has generted!')
    #logger.info(hdr_sciimg)
    fileguidelist=[]
    fileguidelist.append(str(opath+str(filename1)+'_fcimg.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_mskimg.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_flgimg.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_bkg.fits'))
    fileguidelist.append(str(opath+str(filename1)+'_sciimg.fits'))
    
    # #释放局部变量，但是不知道会不会有异常
    # logger.info(locals().keys())
    # for x in list(locals().keys()):
    #     try:
    #         del locals()[x]
    #     except:
    #         continue
    # gc.collect()
    
    return fileguidelist

@logger.catch
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
    #set logging
    skey = str(date) + '_' + tid_target_filter + '_prep'
    loggn = logdir + skey + ".log"
    if os.path.exists(loggn): os.system("rm %s" % loggn)
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logFile = logging.FileHandler(loggn)
    logFile.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger(skey).addHandler(logFile)
    loggerpre = logging.getLogger(skey)

    #loggerpre.info("^_^ Telescope_filter_target: %s" % (tid_target_filter))
    loggerpre.info("[_yprep] "+"^_^ Telescope_filter_target: %s"%(tid_target_filter))

    t1=ts.time()
    tcal=t1-t0
    loggerpre.info(f"[_yprep] "+'the calibration is done,the time cost is:{str(tcal)}')
    loggerpre.info(f"[_yprep] "+'the calibration is done,the time cost is:{str(tcal)}')
    loggerpre.info(f"[_yprep] "+"the calibration is done,the time cost is:{str(tcal)}")
    loggerpre.info(f"[_yprep] "+"the calibration is done,the time cost is: is:{str(tcal)}")
    filelist=np.load(fileguide_raw_path+'sci_'+tid_target_filter+'.npy')
    #filelist= glob.glob(rawinpath+tid+'_'+target+'_'+filterid+'_'+'*.fit')
    nfile=len(filelist)
    
    #record the fileguidelist for different type 
    fgl_fcimg=[]
    fgl_mskimg=[]
    fgl_flgimg=[]
    fgl_uncimg=[]
    fgl_bkg=[]
    fgl_sciimg=[]
    
    
    for i in range(0,nfile):
        dirname,filename=os.path.split(filelist[i])
        subpath = datename(filelist[i]) 
        tscidir = scipath +  subpath +'/'
        mkdir(tscidir)
        try:
            dirname,filename=os.path.split(filelist[i])
            listc=combine_process(filename,dirname+'/',calpath,tscidir,libpath,date)
            #list2=combine_process(filename,dirname+'/',calpath,scipath,2)
            fgl_fcimg.append(listc[0])
            fgl_mskimg.append(listc[1])
            fgl_flgimg.append(listc[2])
            

            fgl_bkg.append(listc[3])
            fgl_sciimg.append(listc[4])

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
 
    
    
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_fcimg.npy',fgl_fcimg)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_mskimg.npy',fgl_mskimg)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_flgimg.npy',fgl_flgimg)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_bkg.npy',fgl_bkg)
    ygn.savenpy(fileguide_sci_path,tid_target_filter+'_sciimg.npy',fgl_sciimg)
    
    
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
    

