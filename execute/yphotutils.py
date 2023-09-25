
import datetime
import gc
import logging
import math
import os
import sys
import time
import traceback
import warnings
from collections import Counter
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

import lib.phot.PlotFun as plotFun
from execute.flatdiff import mkdir
from lib.LogInstance import Logger

from lib.phot.yrefsed import refSed_gaia
import subprocess
import json
import requests

def make_get_request(url, params):
    logger.error(f'make_get_request({url}, {params})')


    headers = {'Content-Type': 'application/json'}
    timeout = 10
    
    try:
        # 将数据转换为JSON格式
        json_data = json.dumps(params)
        
        # 发送POST请求，并指定超时时间
        response = requests.post(url, data=json_data, timeout=timeout, headers={'Content-Type': 'application/json'})

        if response.status_code == 200:
            logger.error(response.content)
            return response.content
        else:
            logger.error(response.status_code)
            return 'An error occurred:', response.status_code
    
    except requests.exceptions.Timeout:
        # 捕获超时异常，并返回错误消息
        return 'Request timeout'
    
def mkdir(path):
    folder = os.path.exists(path)
    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径

def ttfname_split(ttfname):
    
    ttflist=ttfname.split('_')
    logger.info(ttflist)
   
    lenttf = len(ttflist)
    logger.info(lenttf)
    tid= ttflist[0]
    filterid = ttflist[lenttf-1]
    if(lenttf>3):
        for i in range (1,lenttf-2):
            target = ttflist[i]+'_'
    else:
        target = ttflist[1]
    return tid,target,filterid

 

def ttfname_header(filename):
    hdr=fits.getheader(filename)
    #OBJECT  = 'y50b_HZ44_'
    # logger.error(hdr)     ###202309121713 zhangjh 
    tt=hdr['OBJECT']
    try: 
        tid = hdr['TELEID'].strip().lower()
        target = hdr['OBJECT'].strip()
        filterid = 'm'+hdr['FILTER'].strip()
    except:
        
        # logger.error(traceback.format_exc())     ###202309121713 zhangjh
        logger.info(f"tt:{tt}")
        if str(filename).lower().__contains__('50a'):
            tid="y50a"
            target =tt
        elif str(filename).lower().__contains__('50b'):
            tid="y50b"
            target =tt
        else:
            tid= tt.split('_')[0]
            target =tt.split('_')[1]    
        
        filterid='m'+hdr['FILTER']
        logger.info(f"{tid.replace(' ',''),target.replace(' ',''),filterid.replace(' ','')}")
    return tid.replace(' ',''),target.replace(' ',''),filterid.replace(' ','')



def astronet_comd(astrocomd):
     
    astroComdlist=[]
    astroComdlist=astrocomd.split(" ")
    #try:
    if(1==1):
        # 执行命令并设置超时时间为300秒
        completed_process_astrocomd = subprocess.run(astroComdlist, timeout=300, shell=False, capture_output=True)
        #completed_process_astrocomd = subprocess.call(astroComdlist, timeout=300, shell=False, capture_output=True)
        #logger.info(iscampComdlist)           
        # 处理命令执行结果
        if completed_process_astrocomd.returncode == 0:
           
            logger.info(completed_process_astrocomd.stdout)
        else:
            logger.info(completed_process_astrocomd.stderr)
            logger.info('XXXXXXXXXXXXXXXastro failedXXXXXXXXXXXXXXXX')
    
    # except subprocess.TimeoutExpired:
    #     logger.info("Command scamp_comd execution timed out.")
    #     return None, None


def iwcs(iobjFrame_o0):
    wcsfile=iobjFrame_o0[:-5]+'.wcs'
    
    if os.path.exists(iobjFrame_o0[:-5]+'.head'):
        ihead = iobjFrame_o0[:-5]+'.head'
        try:
            os.system("rm %s" % ihead)
        except:
            logger.info(ihead +' delete failed')
    if os.path.exists(iobjFrame_o0[:-5]+'.ahead'):
        ahead = iobjFrame_o0[:-5]+'.ahead'
        try:
            os.system("rm %s" % ahead)
        except:
            logger.info(ahead +' delete failed')
    iimgMat,hdr_sciimg = fits.getdata(iobjFrame_o0,header=True)
    if "PV1_1" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_1"]
    if "PV1_2" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_2"]
    if "PV1_0" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_0"]
    if "PV1_4" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_4"]
    if "PV1_5" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_5"]
    if "PV1_6" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_6"]
    if "PV1_7" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_7"]
    if "PV1_8" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_8"]
    if "PV1_9" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_9"]
    if "PV1_10" in hdr_sciimg.keys():
        del hdr_sciimg["PV1_10"]
    if "PV2_1" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_1"]
    if "PV2_2" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_2"]
    if "PV2_0" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_0"]
    if "PV2_4" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_4"]
    if "PV2_5" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_5"]
    if "PV2_6" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_6"]
    if "PV2_7" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_7"]
    if "PV2_8" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_8"]
    if "PV2_9" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_9"]
    if "PV2_10" in hdr_sciimg.keys():
        del hdr_sciimg["PV2_10"]
        
        
    if(os.path.exists(wcsfile)):
        try:
            logger.info(wcsfile)
            
            hdr_wcs = fits.getheader(wcsfile)
            logger.info(hdr_wcs["CD1_1"])
            logger.info(hdr_wcs["CD1_2"])
            logger.info(hdr_wcs["CD2_1"])
            logger.info(hdr_wcs["CD2_2"])

            hdr_sciimg["CTYPE1"] = "RA---TPV"
            hdr_sciimg["CTYPE2"] = "DEC--TPV"
            hdr_sciimg["RADESYS"] = "ICRS"
            hdr_sciimg["CUNIT1"] = "deg"
            hdr_sciimg["CUNIT2"] = "deg"
            
            hdr_sciimg["CRVAL1"]  = hdr_wcs["CRVAL1"]
            hdr_sciimg["CRVAL2"]  = hdr_wcs["CRVAL2"]
            
            hdr_sciimg["CRPIX1"]  = hdr_wcs["CRPIX1"]
            hdr_sciimg["CRPIX2"]  = hdr_wcs["CRPIX2"]
            

            hdr_sciimg["CD1_1"]   = hdr_wcs["CD1_1"]
            hdr_sciimg["CD1_2"]   = hdr_wcs["CD1_2"]
            hdr_sciimg["CD2_1"]   = hdr_wcs["CD2_1"]
            hdr_sciimg["CD2_2"]   = hdr_wcs["CD2_2"]
            
            fits.writeto(iobjFrame_o0, iimgMat, hdr_sciimg, overwrite=True)

        except:
            logger.info('these are no wcs of in header of this fits')
            logger.error(str(traceback.format_exc()))
 
def solveWcs(x, y, mag, ra, dec, pwd,  width, height, config_cfg,order = 1, pixelScale = 0.429):
            '''
            solve_wcs(x,y,mag,ra,dec,pwd,width=4096,height=4096,pixelScale=1,PV=False)
            x,y,mag        : array float32
            ra,dec         : float32
            pwd            : temprary file path
            width          : image width
            height         : image height
            pixelScale     : guess Pixel Scale
            PV             : use PV standard (for swarp) or not
            sip_tpv calculate the PV 
            '''
            # logger.info('Start Astrometry with Astrometry.net')
            cfgPath = config_cfg
            #CfgName = os.path.join(cfgPath, 'SFOV.cfg' )
            # logger.info('cfgName %s',CfgName)
            tmpcat = Table()
            tmpcat.add_columns([
                Column(data = x,name = 'X'),
                Column(data = y,name = 'Y'),
                Column(data = np.exp(-(mag-25) / 2.5),name = 'FLUX')])
            phdu = fits.PrimaryHDU()
            phdu.header['IMAGEW'] = width
            phdu.header['IMAGEH'] = height
            phdu.header['ANXCOL'] = 'X'
            phdu.header['ANYCOL'] = 'Y'
            # logger.info('IMAGEW %i',width)
            # logger.info('IMAGEH %i',height)
            # logger.info('ANXCOL X')
            # logger.info('ANYCOL Y')
            tmpname = 'tmp.xy'
            wcsname = 'tmp.wcs'
            fits.HDUList([phdu,fits.BinTableHDU(tmpcat)]).writeto(os.path.join(pwd,tmpname),overwrite=True)
            scale_low = pixelScale*0.75
            scale_high = pixelScale*1.5
            # logger.info('Scale Low %.2f Scale High %.2f ',scale_low,scale_high)
            cmd = ['solve-field', '--no-plots', 
                    '--config', cfgPath,
                    '--x-column', 'X', '--y-column', 'Y',
                    '--sort-column','FLUX' 
                    # '--no-remove-lines', '--uniformize', '0',
                    # only work on brightest sources
                    #'--objs', '10000','--no-verify', 
                    '--width', str(width), '--height', str(height),           
                    # '--keep-xylist', sexcat,
                    # '--depth', '50,150,200,250,300,350,400,450,500',
                    '--tweak-order','%i'%order,'--crpix-center',
                    '--scale-low', str(scale_low),'--scale-high', str(scale_high), '--scale-units', 'app',
                    '--overwrite', os.path.join(pwd,tmpname),'--wcs', os.path.join(pwd,wcsname), '--temp-dir', pwd]
             
            if ra != '' and dec != '':
                cmd += [ '--ra', str(ra), '--dec',str(dec), '--radius', str(pixelScale*max(width,height)*3 / 3600.)]
         
            if os.path.isfile(os.path.join(pwd,wcsname)):
                # logger.info('Found WCS File at %s',os.path.join(pwd,wcsname))
                wcshdu = fits.open(os.path.join(pwd,wcsname))[0]
                if 'A_ORDER' in wcshdu.header:
                    logger.info('Success')
                else:
                    logger.info('Failed')
                return wcshdu 
            else:
                logger.info('Failed')
                return fits.PrimaryHDU([]) 


def wcs(iobjFrame_o0,config_cfg):
    ihdr = fits.getheader(iobjFrame_o0)
    if os.path.exists(iobjFrame_o0[:-5]+'.head'):
        ihead = iobjFrame_o0[:-5]+'.head'
        try:
            os.system("rm %s" % ihead)
        except:
            logger.info(ihead +' delete failed')
    if os.path.exists(iobjFrame_o0[:-5]+'.ahead'):
        ahead = iobjFrame_o0[:-5]+'.ahead'
        try:
            os.system("rm %s" % ahead)
        except:
            logger.info(ahead +' delete failed')
    
 
    if "TEL_RA" in ihdr.keys():
        tra= str(ihdr['TEL_RA']) 
        tdec=str(ihdr['TEL_DEC'])    
    # astroComd = "time solve-field %s --scale-low 0.4 --scale-high 0.5  --scale-units app --config %s  --use-sextractor --no-plots --overwrite"    
    # isastroComd = astroComd%(iobjFrame_o0,config_cfg)#,str(tra),str(tdec))
    
    ##
    
    astroComd = "time solve-field %s --scale-low 0.4 --scale-high 0.6 --scale-units app --config %s --ra %s --dec %s --radius 3 --use-sextractor --no-plots --overwrite"    
    isastroComd = astroComd%(iobjFrame_o0,config_cfg,tra,tdec)
     
    ##
    logger.info(isastroComd)
    astronet_comd(isastroComd)
    
     
    try:
            os.system("rm %s" % iobjFrame_o0[:-5]+'.match')
            os.system("rm %s" % iobjFrame_o0[:-5]+'.new')
            os.system("rm %s" % iobjFrame_o0[:-5]+'.rdls')
            os.system("rm %s" % iobjFrame_o0[:-5]+'.axy')
            os.system("rm %s" % iobjFrame_o0[:-5]+'.corr')
            os.system("rm %s" % iobjFrame_o0[:-5]+'-indx.xyls')
            os.system("rm %s" % iobjFrame_o0[:-5]+'.solved')
            
    except:
            logger.error('there is no any pdf file')
        
    wcsfile=iobjFrame_o0[:-5]+'.wcs'
    iimgMat,hdr_sciimg = fits.getdata(iobjFrame_o0,header=True)
    #logger.info(hdr_sciimg)
    if(os.path.exists(wcsfile)):
        try:
            
            logger.info(wcsfile)
            
            hdr_wcs = fits.getheader(wcsfile)
            logger.info(hdr_wcs["CD1_1"])
            logger.info(hdr_wcs["CD1_2"])
            logger.info(hdr_wcs["CD2_1"])
            logger.info(hdr_wcs["CD2_2"])
            logger.info(f'sci_CRPIX1='+str(hdr_sciimg["CRPIX1"]))
            logger.info(f'sci_CRPIX2='+str(hdr_sciimg["CRPIX2"]))
            logger.info(f'wcs_CRPIX1='+str(hdr_wcs["CRPIX1"]))
            logger.info(f'wcs_CRPIX2='+str(hdr_wcs["CRPIX2"]))
            hdr_sciimg['EQUINOX']=2000.0
            
            
            
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
            
            fits.writeto(iobjFrame_o0, iimgMat, hdr_sciimg, overwrite=True)

        except:
            logger.info('these are no wcs of in header of this fits')
            logger.error(str(traceback.format_exc()))
  

def read_list(ipath,filelistpath):
    filelist=np.load(filelistpath)
    imglist=[]
    for i in range (0,len(filelist)):
        img=fits.open(ipath+filelist[i])[1].data
        imglist.append(img)
    return imglist

def timeplus(stime,sec):
    #today = '2021-04-29 19:44:27'
    a=str(stime[:-9])
    b=str(stime[11:])
    st = a+' '+b
    st=datetime.datetime.strptime(st,'%Y-%m-%d %H:%M:%S')
    # 计算偏移量
    offset = datetime.timedelta(seconds=sec)
    # 获取修改后的时间并格式化
    re_date = (st + offset)
    re=str(re_date)
    aa=re[:-9]
    bb=re[11:]
    re=aa+'T'+bb
    return re 
 
def HJDTrans(obj_time,exphalf,tra,tdec):
    objtime=timeplus(obj_time,exphalf)
    t=Time(objtime)
    tjd=Time(np.array(str(t.jd)),format='jd')
    logger.info(tjd)
    gmg=coord.EarthLocation.from_geodetic(100.8667,26.7089,height=3200)#东经100 °01′51″，北纬26 °42′32″，海拔3200米
    c=coord.SkyCoord(tra*u.degree,tdec*u.degree,frame='icrs')
    ltt_helio=tjd.light_travel_time(c,'heliocentric',location=gmg)
    tb=tjd.tdb+ltt_helio
    return tb

def get_scamp_head(headfilename):
    headList = open(headfilename,"r").read().splitlines()
    listh=[]
    for i in range(4,len(headList)-1):
        name=headList[i].split('=')
        name2=name[0].strip()
        value=name[1].strip()
        listh.append((name2,value.split('/')[0]))
    dic=dict(listh)
    return dic

def overlapRect(ra, rb):
    """
    If two rectangles are overlapped, return the 'area' of the 
    overlapping region;
    or return None

    Parameters:
    ra, rb: rectangle edges
            (x0, y0, width_x0, width_y0)
            (x0, y0, x1, y1)
    """
    rax0, rax1 = ra[0], ra[2]
    ray0, ray1 = ra[1], ra[3]

    rbx0, rbx1 = rb[0], rb[2]
    rby0, rby1 = rb[1], rb[3]

    dx = min(rax1, rbx1) - max(rax0, rbx0)
    dy = min(ray1, rby1) - max(ray0, rby0)
    if (dx>=0) and (dy>=0): return dx*dy

def pointRect(p, rt):
    """
    If a point in the rectangle region, return True;
    or return False

    Parameters:
    p: a point -- e.g. (x0,y0)
    rt: a rectangle -- e.g. 
        rt = (x0, y0, x1, y1)
    """
    px, py = p[0], p[1]

    rx0, rx1 = rt[0], rt[2]
    ry0, ry1 = rt[1], rt[3]

    dx = (px-rx0) * (px-rx1)
    dy = (py-ry0) * (py-ry1)

    if (dx<=0) and (dy<=0):
        return True
    else:
        return False


def _spherical_to_cartesian(ra, dec):
    """
    (Private internal function)
    Inputs in degrees. Outputs x,y,z
    """
    rar = np.radians(ra)
    decr = np.radians(dec)

    x = np.cos(rar) * np.cos(decr)
    y = np.sin(rar) * np.cos(decr)
    z = np.sin(decr)

    return x, y, z

def _great_circle_distance(ra1, dec1, ra2, dec2):
    """
    (Private internal function)
    Returns great ciircle distance. Inputs in degrees.

    Uses vicenty distance formula - a bit slower than others, but
    numerically stable.
    """
    lambs = np.radians(ra1)
    phis = np.radians(dec1)
    lambf = np.radians(ra2)
    phif = np.radians(dec2)

    dlamb = lambf - lambs

    numera = np.cos(phif) * np.sin(dlamb)
    numerb = np.cos(phis)*np.sin(phif) - np.sin(phis)*np.cos(phif)*np.cos(dlamb)
    numer = np.hypot(numera, numerb)
    denom = np.sin(phis)*np.sin(phif) + np.cos(phis)*np.cos(phif)*np.cos(dlamb)
    return np.degrees(np.arctan2(numer, denom))*3600.0 # convert to arcsecond

def d2hms(ra,dec, conv=0):
    """
    convert ra&dec in (degree, degree) to (hhmmss, ddmmss),
    or (hhmmss, ddmmss) to (degree, degree)

    Parameters:
    ra, dec: 
       if (ra, dec) is in (deg, deg), float
       if in (hms, dms), string
    conv: 0 or 1
       0: (degree, degree) to (hhmmss, ddmmss)
       1: (hhmmss, ddmmss) to (degree, degree)
    
    Example:
    d2hms(1.0, 1.0, conv=0)
    d2hms(00:00:00.0, 00:00:00.0, conv=1)
    """
    if conv==0:
        rah = ra/15.0
        ram = (rah - int(rah))*60.0
        ras = (ram - int(ram))*60.0
        if rah < 10:
            rah = "0%d:"%int(rah)
        else:
            rah = "%d:"%int(rah)

        if ram < 10:
            ram = "0%d:"%int(ram)
        else:
            ram = "%d:"%int(ram)

        if ras < 10:
            ras = "0%.4f"%float(ras)
        else:
            ras = "%.4f"%float(ras)

        sra = rah+ram+ras

        decabs = abs(dec)
        dech = int(decabs)
        decm = (decabs - dech)*60.0
        decs = (decm - int(decm))*60.0
        if dech < 10:
            dech = "0%d:"%int(dech)
        else:
            dech = "%d:"%int(dech)

        if decm < 10:
            decm = "0%d:"%int(decm)
        else:
            decm = "%d:"%int(decm)

        if decs < 10:
            decs = "0%.4f"%float(decs)
        else:
            decs = "%.4f"%float(decs)

        sdec=dech+decm+decs
        if dec < 0.0:
            sdec = "-"+sdec
        else:
            sdec = "+"+sdec
        return sra, sdec

    elif conv==1:
         
        decSign = str(dec)[0]
        sra = np.array(ra.split(":"), dtype=float)
        sdec = np.array(dec[1:].split(":"), dtype=float)
        sra = ((sra[-1]/60.0+sra[1])/60.0 + sra[0])*15.0
        if decSign=="-":
            sdec = -((sdec[-1]/60.0+sdec[1])/60.0 + abs(sdec[0]))
        elif decSign=="+":
            sdec = (sdec[-1]/60.0+sdec[1])/60.0 + sdec[0]
        else:
            raise ValueError("!!! Give a right dec value")
        return sra, sdec

def mad(data):
    """
    Median absolute deviation, which is defined as
    MAD = median(abs(data-median(data)))

    If data follows normal distributon, then the relationship between 
    the standard deviation (std) of the data and MAD is
    std = 1.4826*MAD

    Return: Normal like-MAD
    """
    mm = np.median(data)
    mx = np.median(abs(data-mm))
    return 1.4826*mx

def wds9reg(x,y,flag=None,radius=5.0,unit="arcsec",color="green",outfile="out.reg"):
    """
    Write ds9 region file.

    Parameters:
    coordinate: 2D array
       coordinate to be written.  It could be image coordinates or RA/Dec.
       Former, set unit="pixel"
       Later, set unit="arcsec"
    radius: float
       in unit of pixels or arcsec (when 'unit' keyword is set)
    unit: string
        pixel: write region file in unit of pixels
        arcsec (default): write region file in unit of RA and Dec
    color: string
       to specify which color to use.  Default is green
    outfile: string
       name of output region file

    Return:
       "outfile": can be read by ds9

    Example:
        pos = [100, 200]
        wds9reg(pos,outfile="pos.reg")
    """

    if not unit in ["arcsec","pixel"]:
        raise ValueError("!! Please set 'unit' as 'arcsec' or 'pixel'")

    fileobj = open(outfile, "w")
    note0 = "# Region file for DS9\n"
    global_pro1 = "global color=%s font='helvetica 10 normal' "%color
    global_pro2 = "select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n"
    fileobj.write(note0)
    fileobj.write(global_pro1+global_pro2)

    if unit == "arcsec":
        fileobj.write("fk5\n")
        fmt = 'circle(%10.6f,%10.6f,%5.2f")\n'
        if flag is not None: fmt='circle(%10.6f,%10.6f,%5.2f") # text={%d}\n'
    if unit == "pixel":
        fileobj.write("image\n")
        fmt = 'circle(%10.6f,%10.6f,%5.2f)\n'
        if flag is not None: fmt='circle(%10.6f,%10.6f,%5.2f) # text={%d}\n'

    for i in range(len(x)):
        if flag is not None:
            ds9row = fmt%(x[i],y[i],radius,flag[i])
        else:
            ds9row = fmt%(x[i], y[i], radius)
        fileobj.write(ds9row)
    fileobj.close()

    return

def read_param(filename):
    pfile = open(filename).readlines()
    nn = len(pfile)
    param = {} # define dictionary structure
    for i in range(nn):
        rowlist = pfile[i].split()
        if len(rowlist)<=1: continue # blank row
        if not "#" in rowlist:
            if len(rowlist)==2:
                key, value = rowlist[0:2]
                param.update({key:value})
            else:
                logger.info("!! Something is wrong with parameter '%s'."%rowlist[0])
                return
        elif rowlist.index("#")==2:
            key, value = rowlist[0:2]
            param.update({key:value})
        elif rowlist.index("#")==0:
            continue # annotation
        else:
            logger.info("!! Something is wrong with parameter '%s'."%rowlist[0])
            return
    return param

def crossmatch1(ra1, dec1, ra2, dec2, aperture=1.0):
    """
    Match two sets of on-sky coordinates to each other.
    I.e., find nearest neighbor of one that's in the other.
    """
    """
    Finds matches in one catalog to another.
    
    Parameters
    ra1 : array-like
          Right Ascension in degrees of the first catalog
    dec1 : array-like
          Declination in degrees of the first catalog (shape of array must match `ra1`)
    ra2 : array-like
          Right Ascension in degrees of the second catalog
    dec2 : array-like
          Declination in degrees of the second catalog (shape of array must match `ra2`)
    aperture : cross-matching aperture, float, default 1.0"
                How close (in arcseconds) a match has to be to count as a match.
    Returns
    -------
    idx1 : int array
           Indecies into the first catalog of the matches. Will never be
           larger than `ra1`/`dec1`.
    idx2 : int array
           Indecies into the second catalog of the matches. Will never be
           larger than `ra1`/`dec1`.
    """
    ra1 = np.array(ra1, copy=False)
    dec1 = np.array(dec1, copy=False)
    ra2 = np.array(ra2, copy=False)
    dec2 = np.array(dec2, copy=False)

    if ra1.shape != dec1.shape:
        raise ValueError('!! ra1 and dec1 do not match!')
    if ra2.shape != dec2.shape:
        raise ValueError('!! ra2 and dec2 do not match!')

    nobj1, nobj2 = len(ra1), len(ra2)
    if nobj1 > nobj2:
        ra1, ra2 = ra2, ra1
        dec1, dec2 = dec2, dec1

    x1, y1, z1 = _spherical_to_cartesian(ra1.ravel(), dec1.ravel())
    coords1 = np.empty((x1.size, 3))
    coords1[:, 0], coords1[:, 1], coords1[:, 2] = x1, y1, z1

    x2, y2, z2 = _spherical_to_cartesian(ra2.ravel(), dec2.ravel())
    coords2 = np.empty((x2.size, 3))
    coords2[:, 0], coords2[:, 1], coords2[:, 2] = x2, y2, z2

    kdt = ckdt(coords2)
    idxs2 = kdt.query(coords1)[1]

    ds = _great_circle_distance(ra1, dec1, ra2[idxs2], dec2[idxs2]) # in arcsecond
    idxs1 = np.arange(ra1.size)

    msk = ds < aperture
    idxs1, idxs2 = idxs1[msk], idxs2[msk]
    ds = ds[msk]
    # Usually, there is duplicate ID in idxs2, here we only keep the one with smaller distance
    dupid = [xx for xx, yy in Counter(idxs2).items() if yy > 1]
    badid = np.array([])
    if len(dupid) > 0:
        for k in dupid:
            kid = np.where(idxs2==k)[0]
            nkid = np.delete(kid,ds[kid].argmin())
            badid = np.append(badid,nkid)
    
    if len(badid)>0: idxs1, idxs2 = np.delete(idxs1,badid), np.delete(idxs2,badid)

    if nobj1 > nobj2:
        newid = np.argsort(idxs2)
        idxs1, idxs2 = idxs2[newid], idxs1[newid]

    return idxs1, idxs2

def crossmatch(ra1, dec1, ra2, dec2, aperture=1.0):
    """
    Match two sets of on-sky coordinates to each other.
    I.e., find nearest neighbor of one that's in the other.
    """
    """
    Finds matches in one catalog to another.
    Parameters
    ra1 : array-like
          Right Ascension in degrees of the first catalog
    dec1 : array-like
          Declination in degrees of the first catalog (shape of array must match `ra1`)
    ra2 : array-like
          Right Ascension in degrees of the second catalog
    dec2 : array-like
          Declination in degrees of the second catalog (shape of array must match `ra2`)
    aperture : cross-matching aperture, float, default 1.0"
                How close (in arcseconds) a match has to be to count as a match.
    Returns
    -------
    idx1 : int array
           Indecies into the first catalog of the matches. Will never be
           larger than `ra1`/`dec1`.
    idx2 : int array
           Indecies into the second catalog of the matches. Will never be
           larger than `ra1`/`dec1`.
    """
    ra1 = np.array(ra1, copy=False)
    dec1 = np.array(dec1, copy=False)
    ra2 = np.array(ra2, copy=False)
    dec2 = np.array(dec2, copy=False)
    if ra1.shape != dec1.shape:
        raise ValueError('!! ra1 and dec1 do not match!')
    if ra2.shape != dec2.shape:
        raise ValueError('!! ra2 and dec2 do not match!')
    nobj1, nobj2 = len(ra1), len(ra2)
    if nobj1 > nobj2:
        ra1, ra2 = ra2, ra1
        dec1, dec2 = dec2, dec1
    x1, y1, z1 = _spherical_to_cartesian(ra1.ravel(), dec1.ravel())
    coords1 = np.empty((x1.size, 3))
    coords1[:, 0], coords1[:, 1], coords1[:, 2] = x1, y1, z1
    x2, y2, z2 = _spherical_to_cartesian(ra2.ravel(), dec2.ravel())
    coords2 = np.empty((x2.size, 3))
    coords2[:, 0], coords2[:, 1], coords2[:, 2] = x2, y2, z2
    kdt = ckdt(coords2)
    idxs2 = kdt.query(coords1)[1]
    ds = _great_circle_distance(ra1, dec1, ra2[idxs2], dec2[idxs2]) # in arcsecond
    idxs1 = np.arange(ra1.size)
    msk = ds < aperture
    idxs1, idxs2 = idxs1[msk], idxs2[msk]
    ds = ds[msk]
    # Usually, there is duplicate ID in idxs2, here we only keep the one with smaller distance
    dupid = [xx for xx, yy in Counter(idxs2).items() if yy > 1]
    badid = np.array([])
    if len(dupid) > 0:
        for k in dupid:
            kid = np.where(idxs2==k)[0]
            nkid = np.delete(kid,ds[kid].argmin())
            badid = np.append(badid,nkid)
    badid = np.array(badid, dtype=int)
     
    #print(badid)
           
 
            
    if len(badid)>0: idxs1, idxs2 = np.delete(idxs1,badid), np.delete(idxs2,badid)
    if nobj1 > nobj2:
        newid = np.argsort(idxs2)
        idxs1, idxs2 = idxs2[newid], idxs1[newid]
    return idxs1, idxs2




    
def reg(ildacn,regname):
    #ildac:输入的星表文件（路径+文件名） regname：输出的region文件（路径+文件名）
    ildac_b      = fits.getdata(ildacn,ext=2)
    bpx, bpy=ildac_b["X_IMAGE"],ildac_b["Y_IMAGE"]
    try:
        mag = ildac_b["MAG_AUTO_S"]
    except:
        mag =ildac_b["MAG_AUTO"]
 
    regListn   = regname#ildac_b[:-5]+ "_select.reg"
    regList    = open(regListn,"w")
    regList.write('# Region file for DS9'+"\n")
    regList.write('global color=green font=\'helvetica 10 normal\' select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'+"\n")
    regList.write('image'+"\n")
    for i in range(0,len(bpx)):
        sline='circle('+str(bpx[i])+','+ str(bpy[i])+', 8.00)'  +'#  font="helvetica 8 normal" text={' + 'imag:'+str(round(mag[i],2))+'}'
        #logger.info(sline)
        regList.write(sline+"\n")
    regList.close()
    logger.info('the reg is ok')


def ref_expt_select(confdir,filterid,bins,expti):
    if(bins == 1):
        limmag_txt = confdir + 'limmag_bin1.txt'
    if(bins == 2):
        limmag_txt = confdir + 'limmag_bin2.txt'
    if(len(filterid)>1):
        filterid=filterid[-1]
        logger.info(filterid)
    limmag_data = np.loadtxt(limmag_txt,skiprows = 1)
    filter_dic = {'u':1,'v':2,'g':3,'r':4,'i':5,'z':6}

    col = filter_dic[filterid]
    mag = limmag_data[:,col]  
    exptime = np.array(limmag_data[:,0],dtype=float)

    maglim_dic = dict(zip(exptime,mag))
    idx = abs(exptime - expti).argmin()
    maglim_max = maglim_dic[exptime[idx]]
    maglim_max = math.ceil(maglim_max )
    logger.info(f'the exptime is: {exptime[idx]}')
    logger.info(f'the maglim_max is: {maglim_max}')

    return maglim_max

def readfits(filename):
    hdu=fits.open(filename)
    index=len(hdu)
    #print(len(hdu),'######################################',hdu.info())
    data=hdu[index-1].data
    #print(np.shape(data))
    hdr=hdu[0].header
    return hdr,data


def re_hdr_edit_fits(filename):
     
    try:
        iimgMat,ihdr=fits.getdata(filename,header=True)
        
        r,c=np.shape(iimgMat)
    except:
        logger.info(f'the fits file is broken:{filename}')
        logger.info("XXXXXXXXXXXXXXXXX the fits file is broken: %s"%(filename))
        return None, None
    
    #logger.info('#####################44')
    ihdr = fits.getheader(filename)
    # logger.info(filename)
    # logger.info(f'the fits header:{ihdr}')
    # logger.info(filename)
    binfact=int(ihdr['XBINNING'])
    pscale  = 0.297 * binfact # pixel scale in unit of arcsec
    ira='0'
    idec='0'
    try:
        
        if("OBJCTRA" in ihdr.keys()):
            ira  = ":".join(ihdr["OBJCTRA"].split())
            idec = ":".join(ihdr["OBJCTDEC"].split())
            ira, idec = d2hms(ira, idec, conv=1)
        elif("RA" in ihdr.keys()):
            ira  = ":".join(ihdr["RA"].split())
            idec = ":".join(ihdr["DEC"].split())
            ira, idec = d2hms(ira, idec, conv=1)
    except:
        logger.info('these are no parameters of object ra and dec in header of this fits')
        logger.error(str(traceback.format_exc()))
    
     
        ihdr['EPOCH']=2000.0
        ihdr["CTYPE1"]  = "RA---TAN"
        ihdr["CTYPE2"]  = "DEC--TAN"
        ihdr["EQUINOX"] = 2000.0
        ihdr["RADESYS"] = "ICRS"
        ihdr["CUNIT1"]  = "deg"
        ihdr["CUNIT2"]  = "deg"
        ihdr["CRVAL1"]  = ira
        ihdr["CRVAL2"]  = idec
        ihdr["CRPIX1"]  = c/2
        ihdr["CRPIX2"]  = r/2
        ihdr["CD1_1"]   = -pscale/3600.0
        ihdr["CD1_2"]   = 0.0
        ihdr["CD2_1"]   = 0.0
        ihdr["CD2_2"]   = -pscale/3600.0
    #img=fits.getdata(filename)
    fits.writeto(filename,iimgMat,ihdr,overwrite= True)
    
    return ira, idec 


def hdr_edit_fits(filename):
    
    #ihdr = fits.getheader(filename)
     
     
    # hdr['OBJ_RA']  = '05h52m27.51s'       
    # hdr['OBJ_DEC'] = '+15d53m16.6s'       
    # logger.info('#####################33')
    # logger.info(filename)
    try:
        iimgMat,ihdr=fits.getdata(filename,header=True)
        
        r,c=np.shape(iimgMat)
    except:
        logger.info(f'the fits file is broken:{filename}')
        logger.info("XXXXXXXXXXXXXXXXX the fits file is broken: %s"%(filename))
        return None, None
    
    #logger.info('#####################44')
    ihdr = fits.getheader(filename)
    # logger.info(filename)
    # logger.info(f'the fits header:{ihdr}')
    # logger.info(filename)
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
        #logger.info(f'ira={ira}, idec={idec}')
    except:
        logger.info('these are no parameters of object ra and dec in header of this fits')
        logger.error(str(traceback.format_exc()))
    
    # ihdr["SATURAT1"] = ihdr["SATURAT1"]
    # ihdr["SATURAT2"] = ihdr["SATURAT2"]
    # ihdr["SATURATE"] = ihdr["SATURATE"]
    # ihdr["SKYBKG1"] = ihdr["SKYBKG1"]
    # ihdr["SKYBKG2"] =ihdr["SKYBKG2"]
    # ihdr["SKYRMS1"] =ihdr["SKYRMS1"] 
    # ihdr["SKYRMS2"] =ihdr["SKYRMS2"]
    # ihdr.add_history('backgroud estimation',before='SATURAT1')##fy20230215
                
    #ihdr.add_history('preprocessing v20230215',before='SKYRMS2')##fy20230215
    # ihdr.append ('EPOCH',2000.0)
    # ihdr.append ('CTYPE1','RA---TAN')
    # ihdr.append ('CTYPE2','DEC--TAN')
    # ihdr.append ('EQUINOX',2000.0)
    # ihdr.append ('RADESYS','ICRS')
    # ihdr.append ('CUNIT1','deg')
    # ihdr.append ('CUNIT2','deg')
    # ihdr.append ('CRVAL1',ira)
    # ihdr.append ('CRVAL2', idec)
    # ihdr.append ('CRPIX1',c/2)
    # ihdr.append ('CRPIX2',r/2)
    # ihdr.append ('CD1_1',-pscale/3600.0)
    # ihdr.append ('CD1_2',0.0)
    # ihdr.append ('CD2_1',0.0)
    # ihdr.append ('CD2_2',-pscale/3600.0)


    #ihdr["SATURATE"] = np.min([65535.0,icountmax])

    # determine the rotation matrix
    if not "CD1_1" in ihdr.keys():
        ihdr['EPOCH']=2000.0
        ihdr["CTYPE1"]  = "RA---TAN"
        ihdr["CTYPE2"]  = "DEC--TAN"
        ihdr["EQUINOX"] = 2000.0
        ihdr["RADESYS"] = "ICRS"
        ihdr["CUNIT1"]  = "deg"
        ihdr["CUNIT2"]  = "deg"
        ihdr["CRVAL1"]  = ira
        ihdr["CRVAL2"]  = idec
        ihdr["CRPIX1"]  = c/2
        ihdr["CRPIX2"]  = r/2
        ihdr["CD1_1"]   = -pscale/3600.0
        ihdr["CD1_2"]   = 0.0
        ihdr["CD2_1"]   = 0.0
        ihdr["CD2_2"]   = -pscale/3600.0
    #img=fits.getdata(filename)
    fits.writeto(filename,iimgMat,ihdr,overwrite= True)
    
    return ira, idec 
    

def datename(i_obj_file):
    try:
        date=fits.getheader(i_obj_file)['DATE-OBS']
        filename = str(date).replace(':','').replace('-','')
        return filename
    except Exception as e:
        logger.error(str(traceback.format_exc()))
        return None
