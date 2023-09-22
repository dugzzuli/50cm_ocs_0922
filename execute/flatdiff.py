from loguru import logger
from scipy import ndimage, misc
from astropy import stats
from astropy.io import fits
import numpy as np

import os,sys,io
import matplotlib
from lib.phot.ybias import get_current_dir

matplotlib.use('agg')
import matplotlib.pyplot as plt


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径

def readfits(file1,file2,filename):
        img1=fits.getdata(file1)
        img2=fits.getdata(file2)
        img=np.hstack((img1,img2))
        hdr=fits.getheader(file1)
        #logger.info(img)
        fits.writeto(filename,img,header=hdr,overwrite=True)
        return img

def pro(tid,filtername,date):

    rootpath=get_current_dir()
    upath=rootpath+'reception/'+str(date)+'/'+'cal/'
    file1=upath+tid+'_master_flat_'+filtername+'_bin2_o1.fits'
    file2=upath+tid+'_master_flat_'+filtername+'_bin2_o2.fits'
    if (os.path.exists(file1) ):
        filename=upath+tid+'_combine_mflat_'+date+'_'+filtername+'.fits'
        img=readfits(file1,file2,filename)
        plt.figure()
        plt.imshow(img,origin='lower',vmin=0.9,vmax=1.1,cmap='rainbow')
        plt.colorbar()
        #plt.plot([1544,1544],[0,3060],'y')
        plt.xlabel('Xpixel',size=15)
        plt.ylabel('Ypixel',size=15)
        plt.title(tid+'_master_flat_'+filtername+'_'+date)
        plt.savefig(upath+tid+'combine_mflat_'+date+'_'+filtername+'.pdf')
        plt.close()
        logger.info('this combine flat process is done')
    else:
        logger.info(tid+'_master_flat_'+filtername +' is not exist')


def diffdate(tid,filtername,date1,date2):
    rootpath=get_current_dir()
    upath=rootpath+'reception/'+str(date1)+'/'+'cal/'
    upath1=rootpath+'reception/'+str(date1)+'/'+'cal/'
    upath2=rootpath+'reception/'+str(date2)+'/'+'cal/'
    #pro(tid,filtername,date1)
    #pro(tid,filtername,date2)
    filename1=upath1+tid+'_combine_mflat_'+date1+'_'+filtername+'.fits'
    filename2=upath2+tid+'_combine_mflat_'+date2+'_'+filtername+'.fits'
    if (os.path.exists(filename1) and os.path.exists(filename2)):
        f1=fits.getdata(filename1)
        f2=fits.getdata(filename2)
        diff_f=f1/f2
        filename_diff=upath+tid+'_diff_mflat_'+filtername+'_'+date1+'-'+date2+'.fits'
        fits.writeto(filename_diff,diff_f,overwrite=True)
         
        plt.figure()
        plt.imshow(diff_f,origin='lower',vmin=0.9,vmax=1.1,cmap='rainbow')
        plt.colorbar()
        #plt.plot([1544,1544],[0,3060],'y')
        plt.xlabel('Xpixel',size=15)
        plt.ylabel('Ypixel',size=15)
        plt.title(tid+'_diff_mflat_'+filtername+'_'+date1+'-'+date2)
        plt.savefig(upath+tid+'_diff_mflat_'+filtername+'_'+date1+'-'+date2+'.pdf')
        plt.close()
        logger.info(tid+ 'of '+ filtername +' in two days is ok')
    else:
        logger.info('there is no co-exist  '+ tid+ '_'+ filtername +' in two days')

def diffplot(date1,date2):
    filtern=['mu','mv','mg','mr','mi','mz']
     
    for i in range(0,6):
        filtername=filtern[i]
        tid='y50a'
        pro(tid,filtername,date1)
        pro(tid,filtername,date2)
        diffdate(tid,filtername,date1,date2)
    for j in range(0,6):
        filtername=filtern[j]
        tid='y50b'
        pro(tid,filtername,date1)
        pro(tid,filtername,date2)
        diffdate(tid,filtername,date1,date2)


# date1=sys.argv[1]
# date2= sys.argv[2]
# diffplot(date1,date2)
