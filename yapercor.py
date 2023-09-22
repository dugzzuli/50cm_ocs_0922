from astropy.io import fits
import numpy as np
import glob
import os, warnings, sys
import matplotlib
import logging
from lib.gchelp import elapsed_time
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pickle
#import cal_airmass
from lib.phot import cal_airmass
#from symfit import Poly, variables, parameters, Model, Fit 
from astropy.modeling import models, fitting
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import vstack, hstack, Table
from astropy.stats import sigma_clip,sigma_clipped_stats
from loguru import logger
from lib.LogInstance import Logger
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures
import math

def get_current_dir():
    fp=os.getcwd()
    filepath,tempfilename = os.path.split(fp)
    filepath=filepath+'/'
    return filepath


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径




def save_object(obj, filename):
        with open(filename, 'wb') as outp:  # Overwrites any existing file.
                pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)


def aper_cor(ttfname,date,scipath,figdir,filename):
   #files = glob.glob(scipath+ttfname+'/'+filename+'*sexcat.fits')
   
#files = glob.glob(ipath+date+'/sci/'+field+'/'+filename+'*.fits')
#print(ipath+date+'/sci/'+field+'/'+filename+'*.fits')
   ttf=ttfname.split('_')
   tid=ttf[0]
   targetid=ttf[1]
   band=ttf[2][-1]
   print('the band is', band)
   try:
         res = fits.open(scipath+filename)
         hdr = res[2].header
         cat = res[2].data
   except:
      print ("Error: 没有找到文件或读取文件失败")
      return None

# pp = PdfPages('aper/'+date+'/'+field+'/'+filename+'_apcorr.pdf')
   
   pp = PdfPages(figdir+filename[:-5]+'_apcorr.pdf')

    

   ind = np.where((np.abs(cat['SPREAD_MODEL']) < 0.06) & (cat['FLAGS'] < 1))
   cat = cat[ind]

   magpsf = cat['MAG_PSF']
   magpsf_err = cat['MAGERR_PSF']
   magauto = cat['MAG_AUTO']
   magauto_err = cat['MAGERR_AUTO']
   logger.info(f"[_maper]:magerr={magauto_err}")
   plt.cla()
   plt.close("all")
   try:
      fig1 = plt.figure()    
      plt.figure(figsize=(14,8))
      plt.plot(magpsf,magpsf_err,'b.')
      #print('ind=',ind)
      ind = np.where((magpsf_err > 0.003) & (magpsf_err < 0.0075))
      maglim1 = 0.003
      maglim2 = 0.0075
      if len(ind) < 150:
         ind = np.where((magpsf_err > 0.002) & (magpsf_err < 0.009))
         maglim1 = 0.002
         maglim2 = 0.009
      if band == 'u':
         ind = np.where((magpsf_err > 0.002) & (magpsf_err < 0.015))
         maglim1 = 0.002
         maglim2 = 0.015

      plt.plot([-17.5,-7.5],[maglim1,maglim1],'r-')
      plt.plot([-17.5,-7.5],[maglim2,maglim2],'r-')
      plt.xlim(-17.5,-7.5)
      plt.ylim(0,0.015)
      plt.grid()
      plt.xlabel('Instrumental mag_psf',fontsize=15)
      plt.ylabel('mag_err',fontsize=15)
      plt.yticks(fontsize=15)
      plt.xticks(fontsize=15)
      plt.savefig(pp, format='pdf')
   except:
      print('Instrumental mag_psf+mag_err figure cannot be generated')

    

   xp = cat['XWIN_IMAGE']
   yp = cat['YWIN_IMAGE']

   fig1 = plt.figure()    
   plt.plot(xp,yp,'ro')
   plt.xlabel('X (pixel)')
   plt.ylabel('Y (pixel)')
   print('Total number of stars selected: ',len(cat))

   "select representative stars across the whole camera:"
   index = []
   ind = np.where((xp < 1500) & (yp < 1500))
   if(len(ind[0])>0):    
      index.append(ind[0][0])
   ind = np.where((xp < 1500) & (yp > 1500))
   if(len(ind[0])>0):  
      #print('ind[0][0]',ind[0][0])  
      index.append(ind[0][0])
   #index.append(ind[0][0])
   ind = np.where((xp > 2500) & (yp > 1500))
   if(len(ind[0])>0):     
      index.append(ind[0][0])
   #index.append(ind[0][0])
   ind = np.where((xp > 2500) & (yp < 1500))
   if(len(ind[0])>0):     
      index.append(ind[0][0])
   #index.append(ind[0][0])
   ind = np.where((xp > 1300) & (xp < 2000) & (yp > 1100) & (yp < 1800))
   if(len(ind[0])>0):     
      index.append(ind[0][0])
   #index.append(ind[0][0])
   ind = np.where((xp > 2100) & (xp < 2800) & (yp > 1100) & (yp < 1800))
   if(len(ind[0])>0): 
      index.append(ind[0][0])
   #index.append(ind[0][0])
   index = np.array(index)
   if(len(index)>0):
      plt.plot(xp[index[0]],yp[index[0]],'b*',markersize=15)
   if(len(index)>1):
      plt.plot(xp[index[1]],yp[index[1]],'k*',markersize=15)
   if(len(index)>2):    
      plt.plot(xp[index[2]],yp[index[2]],'g*',markersize=15)
   if(len(index)>3):    
      plt.plot(xp[index[3]],yp[index[3]],'m*',markersize=15)
   if(len(index)>4):
      plt.plot(xp[index[4]],yp[index[4]],'c*',markersize=15)
   if(len(index)>5):
      plt.plot(xp[index[5]],yp[index[5]],'y*',markersize=15)
   plt.savefig(pp, format='pdf')


   aper = [3,4.5,6,7.5,9,10.5,12,13.5,15,16.5,18,19.5,21,22.5,24,25.5,27,28.5,30,31.5,33,34.5,36,37.5,39,40.5,42,43.5,45,46.5]
   aper = np.array(aper)

   m1 = cat['MAG_APER'][index[0]]
   dm1  = m1[1:29] - m1[0:28]

   fig1 = plt.figure()    
   plt.figure(figsize=(14,8))
   plt.plot(aper[0:28],dm1,'b*',markersize=8)
   plt.plot([1,46.5],[0.0,0.0])
   plt.grid()
   plt.xlabel('Aperture size (pixel)',fontsize=15)
   plt.ylabel('m_(k+1) - m_(k)',fontsize=15)
   plt.yticks(fontsize=15)
   plt.xticks(fontsize=15)
   plt.xlim(2,48)
   plt.ylim(-0.46,0.025)


   if(len(index)>1): 
      m1 = cat['MAG_APER'][index[1]]
      dm2  = m1[1:29] - m1[0:28]
      plt.plot(aper[0:28],dm2,'k*',markersize=8)
   if(len(index)>2):
      m1 = cat['MAG_APER'][index[2]]
      dm3  = m1[1:29] - m1[0:28]
      plt.plot(aper[0:28],dm3,'g*',markersize=8)
   if(len(index)>3):
      m1 = cat['MAG_APER'][index[3]]
      dm4  = m1[1:29] - m1[0:28]
      plt.plot(aper[0:28],dm4,'m*',markersize=8)
   if(len(index)>4):
      m1 = cat['MAG_APER'][index[4]]
      dm5  = m1[1:29] - m1[0:28]
      plt.plot(aper[0:28],dm5,'c*',markersize=8)
   if(len(index)>5):
      m1 = cat['MAG_APER'][index[5]]
      dm6  = m1[1:29] - m1[0:28]
      plt.plot(aper[0:28],dm6,'y*',markersize=8)
   plt.savefig(pp, format='pdf')

   "best aperture size"
   fig1 = plt.figure()    
   plt.figure(figsize=(14,8))
   if(len(index)>0):
      em1 = cat['MAGERR_APER'][index[0]]
      plt.plot(aper,em1,'b*-')
   if(len(index)>1):
      em1 = cat['MAGERR_APER'][index[1]]
      plt.plot(aper,em1,'k*-')
   if(len(index)>2):
      em1 = cat['MAGERR_APER'][index[2]]
      plt.plot(aper,em1,'g*-')
   if(len(index)>3):
      em1 = cat['MAGERR_APER'][index[3]]
      plt.plot(aper,em1,'m*-')
   if(len(index)>4):
      em1 = cat['MAGERR_APER'][index[4]]
      plt.plot(aper,em1,'c*-')
   if(len(index)>5):
      em1 = cat['MAGERR_APER'][index[5]]
      plt.plot(aper,em1,'y*-')
   plt.ylim(0,0.025)
   plt.xlabel('Aperture size (pixel)',fontsize=15)
   plt.ylabel('mag_err',fontsize=15)
   plt.yticks(fontsize=15)
   plt.xticks(fontsize=15)
   plt.grid()
   mfwhm = np.median(cat['FWHM_IMAGE'])
   n = len(cat)
   allap = []
   alldx = []
   for i in range(n):
      magerr = cat['MAGERR_APER'][i]
      ind = np.where(np.abs(magerr) == min(np.abs(magerr)))
      allap.append(aper[ind])
      alldx.append(ind[0][0])
   baper = round(np.median(allap))
   bindx = round(np.median(alldx))
   plt.plot([aper[bindx],aper[bindx]],[0,0.025],'r-')

   print('The bset aperture size (pixel):',baper)
   print('The median FWHM (pixel):',mfwhm)
   print('The ratio between bset aperture size and median FWHM:',baper/mfwhm)
   #print(bindx)
   #print(aper[bindx])
   plt.savefig(pp, format='pdf')

   fig1 = plt.figure()    
   plt.hist(np.array(allap), bins=30)
   plt.xlabel('Aperture size (pixel)')
   plt.ylabel('Number')
   plt.xlim(2, 15)
   plt.savefig(pp, format='pdf')


   magdif = []
   magdif2 = []
   magdif3 = []
   for i in range(n):
      mag = cat['MAG_APER'][i]
      magpsf = cat['MAG_PSF'][i]
      magd = mag[22] - mag[bindx]
      magdif.append(np.array(magd))
      magd2 = mag[bindx] - magpsf
      magdif2.append(np.array(magd2))
      magd3 = mag[22] - magpsf 
      magdif3.append(np.array(magd3))

   magdif = np.array(magdif)
   magdif2 = np.array(magdif2)
   magdif3 = np.array(magdif3)
   #print('magdif=',np.median(magdif))
   ind = np.where((magdif > np.median(magdif) -0.35) & (magdif < np.median(magdif) +0.35))
   magdif = magdif[ind]
   magdif2 = magdif2[ind]
   magdif3 = magdif3[ind]
   xp = xp[ind]
   yp = yp[ind]
   #print(len(xp), len(yp), len(magdif))
   fig1 = plt.figure()    
   plt.scatter(xp, yp, s=255, c=magdif, cmap='rainbow')
   plt.xlabel('X (pixel)')
   plt.ylabel('Y (pixel)')
   plt.colorbar(label='mag_ref - mag_best')
   plt.savefig(pp, format='pdf')
   #print('cccccccccccccc')
   xdata = xp
   ydata = yp
   zdata = magdif
# Fit the data using astropy.modeling
   p_init = models.Polynomial2D(degree=3)
   fit_p = fitting.LevMarLSQFitter()

   with warnings.catch_warnings():
# Ignore model linearity warning from the fitter
      warnings.simplefilter('ignore')
      p = fit_p(p_init, xdata, ydata, zdata)
      resd = zdata-p(xdata, ydata)
      sigma = sigma_clipped_stats(resd)[2]
      ind = np.where(abs(resd-np.median(resd)) < 2.75*sigma)
      xdata = xdata[ind]
      ydata = ydata[ind]
      zdata = zdata[ind]
      magdif2 = magdif2[ind]
      magdif3 = magdif3[ind]
      p = fit_p(p_init, xdata, ydata, zdata)
      sigma = sigma_clipped_stats(zdata-p(xdata, ydata))[2]
      resd = zdata-p(xdata, ydata)
      ind = np.where(abs(resd-np.median(resd)) < 2.75*sigma)
      xdata = xdata[ind]
      ydata = ydata[ind]
      zdata = zdata[ind]
      magdif2 = magdif2[ind]
      magdif3 = magdif3[ind]
      p = fit_p(p_init, xdata, ydata, zdata)

   # print(sigma_clipped_stats(zdata)[2])
   # print(sigma_clipped_stats(zdata-p(xdata, ydata))[2])
   # print(sigma_clipped_stats(magdif2+p(xdata, ydata))[2])
   # print(np.median(magdif2+p(xdata, ydata)))
   fig1 = plt.figure()    
   plt.scatter(xdata, ydata, s=255, c=zdata-p(xdata, ydata), cmap='rainbow')
   plt.xlabel('X (pixel)')
   plt.ylabel('Y (pixel)')
   plt.colorbar(label='mag_ref - mag_best - mfit')
   plt.savefig(pp, format='pdf')

   fig1 = plt.figure()    
   plt.figure(figsize=(14,8))
   plt.plot(xdata,zdata-p(xdata, ydata),'b.')
   plt.ylim(-0.08,0.08)
   plt.xlabel('X (pixel)')
   plt.grid()
   plt.ylabel('mag_ref - (mag_best - mag_corr)')
   plt.savefig(pp, format='pdf')

   fig1 = plt.figure()    
   plt.scatter(xdata, ydata, s=255, c=magdif2+p(xdata, ydata), cmap='rainbow',vmin=-0.05,vmax=0.05)
   plt.xlabel('X (pixel)')
   plt.ylabel('Y (pixel)')
   plt.colorbar(label='mag_best - mag_corr - mag_psf')
   plt.savefig(pp, format='pdf')

   fig1 = plt.figure()    
   plt.figure(figsize=(14,8))
   plt.plot(xdata,magdif2+p(xdata, ydata),'b.')
   plt.ylim(-0.035,0.035)
   plt.xlabel('X (pixel)')
   plt.grid()
   plt.ylabel('mag_best - mag_corr - mag_psf')
   plt.savefig(pp, format='pdf')

   fig1 = plt.figure()
   plt.figure(figsize=(14,8))
   plt.plot(xdata,magdif3,'b.')
   plt.ylim(-0.085,0.085)
   plt.grid()
   plt.xlabel('X (pixel)')
   plt.ylabel('mag_ref - mag_psf')
   plt.savefig(pp, format='pdf')

#save the index of best aperture and fitting results
   res = [bindx,p]
   #save_object(res, 'aper/'+date+'/'+field+'/'+filename+'_apcorr.pkl')
   #scipath+ttfname+'/'+filename[:-5]+'_apcorr.fits'
   save_object(res, scipath+filename[:-5]+'_apcorr.pkl')
   pp.close()

#apply corrections to the entire field:
   res = fits.open(scipath+filename)
   cat = res[2].data
   cat = Table(cat)

   xp = cat['X_IMAGE']
   yp = cat['Y_IMAGE']
   ra  = np.array(cat['ALPHA_J2000'])
   dec = np.array(cat['DELTA_J2000'])

   #  n = len(cat)
   #  apmag = []
   #  apmag_err = []
   #  for i in range(n):
   #      mag = cat['MAG_APER'][i]
   #      magerr = cat['MAGERR_APER'][i]
   #      apmag.append(mag[bindx])
   #      apmag_err.append(magerr[bindx])
   #  apmag = np.array(apmag)
   #  apmag_err = np.array(apmag_err)
   #  mag_faper = apmag + np.array(p(xp,yp))
   fitsname=str(scipath+filename)
   fitsname_sciimg=fitsname.split("_sexcat")[0]
   hdr = fits.getheader(fitsname_sciimg+'.fits')
        
 
   dig="_sexcat"
   fitsname_sciimg=fitsname.split(dig)[0] 
   hdr = fits.getheader(fitsname_sciimg+'.fits')
 
   date_obs = hdr['DATE-OBS']
   expt = float(hdr['EXPTIME'])


   n = len(cat)
   apmag = []
   apmag_err = []
   magerr_faper_s=[]
   aflux=[]
   afluxerr =[]
   for i in range(n):
      mag = cat['MAG_APER'][i]
      magerr = cat['MAGERR_APER'][i]
      flux = cat['FLUX_APER'][i]
      fluxerr = cat['FLUXERR_APER'][i]
      apmag.append(mag[bindx])
      apmag_err.append(magerr[bindx])
      aflux.append(flux[bindx])
      afluxerr.append(fluxerr[bindx])
      iamagerr_s = 1.0857 * ( float(fluxerr[bindx])/float(flux[bindx]))
      magerr_faper_s.append(iamagerr_s)
   apmag = np.array(apmag)
   apmag_err = np.array(apmag_err)
   mag_faper = apmag + np.array(p(xp, yp)) #+ math.log10(expt) * 2.5
   mag_faper_s = mag_faper + math.log10(expt) * 2.5
   #magerr_faper_s  = 1.0857 * ( afluxerr/aflux)
        
     
   airmass = []
   jdm = []
   for i in range(n):
      amss,jd_mid,tbd_mid = cal_airmass.cam(date_obs,expt,ra[i],dec[i])
      airmass.append(amss)
      jdm.append(jd_mid)
    
   dat = Table({'MAG_CORR': np.array(p(xp, yp)), 'MAG_FAPER': mag_faper, 'MAGERR_FAPER': apmag_err, 'MAG_FAPER_S': mag_faper_s, 'MAGERR_FAPER_S': magerr_faper_s,'AIRMASS': airmass,
                    'JD_MID': jdm})
   zat = hstack([cat, dat])
   files = os.path.exists(scipath + '/' + filename[:-5] + '_apcorr.fits')
   if files:
      os.remove(scipath + '/' + filename[:-5] + '_apcorr.fits')
   zat.meta['JD_MID'] = np.median(jd_mid)
   zat.meta['TBD_MID'] = tbd_mid
   res[2].data =np.array(zat)
   
   #hdul.writeto(ildacn, overwrite=True)
   res.writeto(scipath + '/' + filename[:-5] + '_apcorr.fits', overwrite=True)
   res.close()
   #zat.write(scipath + '/' + filename[:-5] + '_apcorr.fits')
   #loggerloguru.info(f"[_maper]:{np.median(dat['MAG_FAPER'])}")
 

def single_apercor_pro( date, ttfname, scipath,figdir,filename1):
   rootpath=get_current_dir()
   loggerloguru = Logger(date, '', "_maper").get_logger
   
   #figdir = rootpath + "figures/" + date + "/" + date + "_" + ttfname + "/"
   #os.mkdir(figdir)
   #loggerloguru.info("[_maper]:"+'figdir='+str( figdir))
   files = scipath + filename1
   #loggerloguru.info("[_maper]:"+f'the aperture correction will begin, the sexcat files are:{files}')
   try:
   #if(1==1):
      filename = str(files).rsplit("/", 1)[1]
      aper_cor(ttfname, date, scipath, figdir, filename)
      loggerloguru.info("[_maper]:"+'the aperture correction is ok!')
   except Exception as e:
      
      loggerloguru.info("[_maper]:"+str(files) + ' is broken')
      loggerloguru.info(f"[_maper]:{e}")
   

@elapsed_time
def pro_old(ttfname,date):
    rootpath=get_current_dir()
    upath=rootpath+'reception/'+str(date)+'/'
    scipath = upath+'sci/'+ttfname+'/'
    print(scipath)
    #print(lpathsci[0:5])
    lpathsciimg=np.load(scipath+'fileguide/'+ttfname+'_sciimg.npy')
    lpathsciimg=np.sort(lpathsciimg)
    imgdir  = rootpath + "images/"
    confdir = rootpath + "config/"
    ancdir  = rootpath + "ancData/"
    figdir  = rootpath + "figures/"+date+"_"+ttfname+"/"
    pdfdir = rootpath+'50cmpy/'
    mkdir(figdir)
    print('figdir=', figdir)
    ttf=ttfname.split('_')
    tid=ttf[0]
    targetid=ttf[1]
    files = glob.glob(scipath+'*sexcat.fits')
    print('the aperture correction will begin, the sexcat files are:',files)
    for item in files:
        try:
           filename=str(item).rsplit("/",1)[1]
           newf=str(item)[:-5]+'_apcorr.fits'
           # if (os.path.exists(newf)):
           #    continue
           # print(filename)
           aper_cor(ttfname,date,scipath,figdir,filename)
        except:
             print(filename+' is broken')
             continue
    print('the aperture correction is processed!')


def pro(date,ttfname):
      rootpath=get_current_dir()

      upath = rootpath + 'reception/' + str(date) + '/'
      scipath = upath + 'sci/' + ttfname + '/'
      figdir = rootpath + "figures/" + date + "/" + date + "_" +ttfname + "/"
      mkdir(figdir)
      loggerloguru = Logger(date, '', "_maper").get_logger
      try:
         tid,objid,filterid = str(ttfname).split('_')
      except:
         tid,objid1,objid2,filterid = str(ttfname).split('_')
         objid = objid1+'_'+objid2
      files = glob.glob(scipath + '*sexcat.fits')
   
      
      #file = glob.glob(fileguide_raw_path+'sci_'+tid+'_*'+objid+'_'+filterid+'.npy')
      sexcatlist = glob.glob(scipath+'*sciimg_sexcat.fits')
      sublist = os.path.join(scipath, str(date)+'T*')
   
      subdirs = [d for d in glob.glob(sublist) if os.path.isdir(d)]

      
      loggerloguru.info("[_maper]:"+'the aperture correction will begin, the sexcat files are:', files)
      if(subdirs):
         with ProcessPoolExecutor(max_workers=5) as executor:
               taskaper=[]
               for j in range(0,len(subdirs)):
                  filenames = glob.glob(subdirs[j]+'/'+tid+'*sciimg_sexcat.fits')
                  if filenames:
                     if(os.path.exists(filenames[0][:-5] + '_apcorr.fits')):
                           loggerloguru.info(filenames[0][:-5] + '_apcorr.fits'+'has been apercor processed')
                           continue
                  if(filenames):
                     filename = str(filenames[0]).rsplit("/", 1)[1]
                     
                     tscipath = subdirs[j]+'/'
                     # if os.path.exists(subdirs[j]+'/'+filename[:-5] + '_apcorr.fits') and config["log"]["record"]:
                     #     loggerloguru.info("[_maper]:"+filename[:-5] + '_apcorr.fits'+' was been apercor processd! ')
                     #     continue
                     # self.aper_cor(self.ttfname, self.date, scipath, figdir, filename)
                     tempaper=executor.submit(aper_cor,ttfname, date, tscipath, figdir, filename)
                     loggerloguru.info("[_maper]:"+f"aper_cor,ttfname:{ttfname}, date:{date}, scipath:{scipath}, figdir:{figdir}, filename:{filename}")
                     loggerloguru.info("[_maper]:"+filename[:-5] + '_apcorr.fits'+' was been apercor processd! ')
                     taskaper.append(tempaper)
      
      elif files:
         with ProcessPoolExecutor(max_workers=10) as executor:
               taskaper=[]
               for item in files:
                  filename = str(item).rsplit("/", 1)[1]
                  if os.path.exists(scipath + '/' + filename[:-5] + '_apcorr.fits') and config["log"]["record"]:
                     loggerloguru.info("[_maper]:"+filename[:-5] + '_apcorr.fits'+' was been apercor processd! ')
                     continue
                  # self.aper_cor(self.ttfname, self.date, scipath, figdir, filename)
                  tempaper=executor.submit(aper_cor,ttfname, date, scipath, figdir, filename)
                  loggerloguru.info("[_maper]:"+f"aper_cor,ttfname:{ttfname}, date:{date}, scipath:{scipath}, figdir:{figdir}, filename:{filename}")
                  loggerloguru.info("[_maper]:"+filename[:-5] + '_apcorr.fits'+' was been apercor processd! ')
                  taskaper.append(tempaper)
      
                  
               if taskaper:
                  
                  concurrent.futures.wait(taskaper)
               
               executor.shutdown()
      else:
         loggerloguru.info(f"[_maper]: len(files):{len(files)}")


      loggerloguru.info("[_maper]:"+'the aperture correction is ok!')
        
if __name__ == "__main__":
    import os
    import argparse
    import time
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default="20230129", help='输入处理日期')
    parser.add_argument('--tid_target_filter', type=str, default="y50a_SC8002_mr", help='')
    args = parser.parse_args()  # 参数解析
    starttime = time.time()
    loggerloguru = Logger(args.date, '', "_yaper").get_logger
    pro(args.tid_target_filter,args.date)
    endtime=time.time()
    loggerloguru.info("[_yapercor] "+"time spending:{}s".format(endtime-starttime))
