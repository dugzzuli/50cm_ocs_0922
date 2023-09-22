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
from lib.gchelp import elapsed_time
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
from astropy.table import Table
from astropy.time import Time
from loguru import logger
from scipy.optimize import leastsq
from scipy.spatial import cKDTree as ckdt

import lib.phot.PlotFun as plotFun
from execute.flatdiff import mkdir
from lib.LogInstance import Logger

from lib.phot.yrefsed import refSed_gaia
from lib.phot import yrefsed




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
        decSign = dec[0]
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
    
    if len(badid)>0: idxs1, idxs2 = np.delete(idxs1,badid), np.delete(idxs2,badid)

    if nobj1 > nobj2:
        newid = np.argsort(idxs2)
        idxs1, idxs2 = idxs2[newid], idxs1[newid]

    return idxs1, idxs2

def reg(ildac_b,regname):
    bpx, bpy=ildac_b["X_IMAGE"],ildac_b["Y_IMAGE"]
    regListn   = regname#ildac_b[:-5]+ "_select.reg"
    regList    = open(regListn,"w")
    regList.write('# Region file for DS9'+"\n")
    regList.write('global color=green font=\'helvetica 10 normal\' select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'+"\n")
    regList.write('image'+"\n")
    logger.info(regname)

    for i in range(0,len(bpx)):
      
        sline='circle('+str(bpx[i])+','+ str(bpy[i])+', 8.00)'  
        #logger.info(sline)
        regList.write(sline+"\n")
    regList.close()
    logger.info('the reg is ok')


def ref_expt_select(confdir,filterid,bins,expti):
    if(bins == 1):
        limmag_txt = confdir + 'limmag_bin1.txt'
    if(bins == 2):
        limmag_txt = confdir + 'limmag_bin2.txt'

    limmag_data = np.loadtxt(limmag_txt,skiprows = 1)
    filter_dic = {'u':1,'v':2,'g':3,'r':4,'i':5,'z':6}

    col = filter_dic[filterid]
    mag = limmag_data[:,col]  
    exptime = np.array(limmag_data[:,0],dtype=float)

    maglim_dic = dict(zip(exptime,mag))
    idx = abs(exptime - expti).argmin()
    maglim_max = maglim_dic[exptime[idx]]
    maglim_max = math.ceil(maglim_max )
    logger.info('the exptime is',exptime[idx])
    logger.info('the maglim_max is',maglim_max)

    return maglim_max

def readfits(filename):
    hdu=fits.open(filename)
    index=len(hdu)
    #print(len(hdu),'######################################',hdu.info())
    data=hdu[index-1].data
    #print(np.shape(data))z
    hdr=hdu[0].header
    return hdr,data

#fy_new 20230131
def reget_refsed(ancdir,objID,fra,fdec):
    #refpath = '%sGaiaDR3_%s.txt'%(refsavepath,targetname)
    refpath=yrefsed.RefCatDownload(ancdir,objID,fra,fdec)
    scamp_ldac=yrefsed.ScampRefCat(refpath,blimmag=0,dlimmag=25)
    return scamp_ldac
     
    

def photometry(rootpath,ttfname,date):
    upath=rootpath+'reception/'+str(date)+'/'
    scipath = upath+'sci/'+ttfname+'/'
    logger.info(scipath+'fileguide/'+ttfname+'_subbkg_o0.npy')
    lpathsubbkg=np.load(scipath+'fileguide/'+ttfname+'_subbkg_o0.npy')
    logger.info(lpathsubbkg)
    lpathsubbkg=np.sort(lpathsubbkg)
    #basic setting
    #dirname, filename = os.path.split(os.path.abspath(__file__))
    imgdir  = rootpath + "images/"
    confdir = rootpath + "config/"
    ancdir  = rootpath + "ancData/"
    figdir  = rootpath + "figures/"+str(date)+'/'+str(date)+'_'+ttfname+"/"
    pdfdir = rootpath+'50cmpy/'
    logdir = rootpath+'run_pipe_log/'+str(date)+'/'+ttfname+'/'
     
    mkdir(logdir)
    mkdir(figdir)   
    ttf=ttfname.split('_')
    tid=ttf[0]
    targetid=ttf[1]
    try:
        
        hrs,imgs=readfits(lpathsubbkg[0])
        logger.info(lpathsubbkg[0])
        ira  = ":".join(hrs["OBJCTRA"].split())
        idec = ":".join(hrs["OBJCTDEC"].split())
    except:
        logger.info ("Error: 没有找到文件或读取文件失败")
        for i in range(1,len(lpathsubbkg)):
            try:
                hrs=fits.getheader(lpathsubbkg[i])
                
                if(len(hrs["OBJCTRA"])>1):
                    ira  = ":".join(hrs["OBJCTRA"].split())
                    idec = ":".join(hrs["OBJCTDEC"].split())
                    ira, idec = d2hms(ira, idec, conv=1)
                    hrs,imgs=readfits(lpathsubbkg[i])
                break
            except:
                logger.info('these are no parameters of object ra and dec in header of this fits')
                continue
              
    try:
        ira  = ":".join(hrs["OBJCTRA"].split())
        idec = ":".join(hrs["OBJCTDEC"].split())
        ira, idec = d2hms(ira, idec, conv=1)
    except:

        logger.info('these are no parameters of object ra and dec in header of this fits')
    
     
    r,c=np.shape(imgs)
    filtID  = ttf[2][-1]
    objID   = ttf[1] 

    binfact = int(hrs['XBINNING'])
    pscale  = 0.297 * binfact # pixel scale in unit of arcsec
    logger.info(pscale)
    saturate_value=hrs['SATURATE']
    exptlast = hrs['EXPTIME']
    refCatn  = ancdir + "GaiaStar_%s_%s.ldac"%(objID,date)
    refCat2n = ancdir + "GaiaStar_%s_%s.cat"%(objID,date)
    #if not (os.path.exists(refCatn)):
    
    maxmaglim=ref_expt_select(confdir,filtID,binfact,exptlast)
    refSed_gaia(rootpath,ttfname,date,ira,idec,10,maxmaglim)
    refCat2  = Table.read(refCat2n, format="ascii")
    raRef    = refCat2["ra"]
    decRef   = refCat2["dec"]
    try:
        # magRef   = refCat2["mag%s"%filtID.upper()]
         magRef  = refCat2["mag"]
    except:
         logger.info('the filteris is '+filtID.upper()+' and use mag_nofilter only')
         magRef  = refCat2["mag"]    
    nrefStar = len(raRef)

    scipath=rootpath+'reception/'+str(date)+'/sci/'+ttfname+'/'
    tscidir = scipath 
    mkdir(tscidir)


    # 3S configuration sexParam1  = confdir + "default.param"
    sexParam1  = confdir + "default.param"
 
    sexParam2  = confdir + "default.paramNew"
    sexConf    = confdir + "default.sex"
    swarpConf  = confdir + "default.swarp"
    scampConf  = confdir + "default.scamp"
    scampConf2  = confdir + "second.scamp"
    psfexConf  = confdir + "default.psfex"
 
     

    # set logging
     

    miniarea = str(6/binfact)
    # logger.info(miniarea)
    sexComd1   = "sex %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s "
    sexComd2   = "-DETECT_MINAREA " +miniarea+ " -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 "
    sexComd3   = "-SATUR_LEVEL %8.1f -CHECKIMAGE_TYPE NONE "
    sexComd4   = "-PSF_NAME %s "
    sexComd    = sexComd1 + sexComd2 + sexComd3
    sexComdPSF = sexComd + sexComd4
    scampComd1 = "scamp %s -c %s -ASTREFCAT_NAME %s -MOSAIC_TYPE LOOSE " #UNCHANGED
    scampComd2 = "-FWHM_THRESHOLDS 2,20 -SN_THRESHOLDS 5,1000 "
    scampComd  = scampComd1 + scampComd2
    psfexComd  = "psfex %s -c %s -PSF_SIZE %d,%d -PSF_DIR %s "
 

    # set logging
    skey  = str(date)+'_'+ttfname+'_phos'
    loggn = logdir + skey + ".log"
    if os.path.exists(loggn): os.system("rm %s"%loggn)
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logFile = logging.FileHandler(loggn)
    logFile.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger(skey).addHandler(logFile)
    loggerpho = logging.getLogger(skey)
    loggerpho.info("^_^ Telescope_filter_target: %s"%(ttfname))
    ldacListn   = tscidir +ttfname +'_o0_ldac.list'
    ldacList    = open(ldacListn,"w")
    subImgListn = tscidir +ttfname+'_o0_imageSub.list'
    subImgListn1 = tscidir +ttfname+'_o1_imageSub.list'
    subImgListn2 = tscidir +ttfname+'_o2_imageSub.list'
    subImgList  = open(subImgListn,"w")
    subImgList1  = open(subImgListn1,"w")
    subImgList2 = open(subImgListn2,"w")
    #fgl_ldac=[]
    fgl_ldac=[]

    #hdrlist=[]


    for i,file_name in enumerate(lpathsubbkg):
        filei=lpathsubbkg[i].rsplit('/',1)[1]
        #logger.info('####### the precessed fits is: '+ filei)
        loggerpho.info('####### the processed fits is: %s'%(filei))
        #对每一个frame run sextractor
        iobjFrameNew = lpathsubbkg[i]#[:-3] + "sub.fit"
        iobjFrame_o0=iobjFrameNew[:-14]+date+'_subbkg_o0.fits'
        ipsfFrameNew = iobjFrame_o0[:-4] + "psf"
        try:
            ihdr,iimgMat= readfits(iobjFrameNew)
            r,c=np.shape(iimgMat)
        except:
            loggerpho.info(f'the fits file is broken:{iobjFrameNew}')
            loggerpho.info("XXXXXXXXXXXXXXXXX the fits file is broken: %s"%(iobjFrameNew))
            continue
 
        binfact=int(ihdr['XBINNING'])
        pscale  = 0.297 * binfact # pixel scale in unit of arcsec

        try:
            ira  = ":".join(ihdr["OBJCTRA"].split())
            idec = ":".join(ihdr["OBJCTDEC"].split())
            ira, idec = d2hms(ira, idec, conv=1)
        except:
            loggerpho.info('these are no parameters of object ra and dec in header of this fits')
        
            logger.error(str(traceback.format_exc()))
            continue
            
         
        if(c>8000 and binfact==1) or (c> 4000 and binfact ==2) :

            iimgMat = iimgMat[:,int(1000/binfact):c-int(1000/binfact)]
            r,c=np.shape(iimgMat)
            loggerpho.info(f"#################imgcut################## :{np.shape(iimgMat)}")
           
        fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True)
         

        ihdr['BUNIT']="AUD/s"
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
        #ihdr["SATURATE"] = np.min([65535.0,icountmax])

        # determine the rotation matrix
        ihdr["CD1_1"]   = pscale/3600.0
        ihdr["CD1_2"]   = 0.0
        ihdr["CD2_1"]   = 0.0
        ihdr["CD2_2"]   = pscale/3600.0
        
        fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True)
        
        exptime=ihdr['EXPTIME']
        # hdr1=fits.getheader(iobjFrame_o0)
        # logger.info(hdr1)

     
        #logger.info('rewrite fits')
        # run sextractor to extract objects
        ildacn = iobjFrame_o0[:-4] + "ldac"
        ildac_select_reg=iobjFrame_o0[:-5] + "_select.reg"
        saturate_value=ihdr['SATURATE']*0.9
        isexComd = sexComd%(iobjFrame_o0,sexConf,sexParam1,ildacn,saturate_value)
        os.system(isexComd)
        ldac=fits.getdata(ildacn,ext=2)
        mean_fwhm=np.median(ldac['FWHM_IMAGE'])
        mean_ell = np.median(ldac['ELLIPTICITY'])
        loggerpho.info(f'mean_fwhm ={mean_fwhm}')
        loggerpho.info(f'mean_ell ={mean_ell}')
        # construct PSF model
        
        loggerpho.info("################# psfexComd##################")
  
        ipsfComd = psfexComd%(ildacn, psfexConf, 36/binfact, 36/binfact, tscidir) #21
        os.system(ipsfComd)
        os.system("mv *pdf %s"%(figdir))
        os.system("mv *fits %s"%tscidir)
 
        loggerpho.info("   Run Scamp to refine the astrometry")
        #warnings.filterwarnings('error') 
         ##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&        
        # run scamp in a limited area
        #warnings.filterwarnings('error') 
        #try:
        
        ildac=fits.getdata(ildacn,ext=2)
            
        if(filtID =='v' or filtID =='u'):
            ixcenloc=np.where((ildac["X_IMAGE"]>c/2+10) | (ildac["X_IMAGE"]<c/2-10))
            #ixyloc=np.where((iximg<c-(1000/binfact))&(iximg>1000/binfact))#&(iyimg<6100)&(iyimg>30))

            ildac =ildac[ixcenloc]
        hdul=fits.open(ildacn)###
        hdul[1].data=ildac

        reg(ildac,ildac_select_reg)
        #isnr         = ildac["SNR_WIN"]
        #ildac["SNR_WIN"]=isnr*exptime
        hdul.close()
            
            
        iscampComd = scampComd%(ildacn,scampConf,refCatn)
        os.system(iscampComd)
        ihead   = ildacn[:-4] + "head"
        ihead1 = ildacn[:-4] + "ahead"
            

        os.system("mv %s %s"%(ihead,ihead1))
        iscampComd2 = scampComd%(ildacn,scampConf2,refCatn)
        os.system(iscampComd2)


        ihead   = ildacn[:-4] + "head"
        iimgMat, ihdr = fits.getdata(iobjFrame_o0,header=True)
        iheadList = open(ihead,"r").read().splitlines()[3:-1]
        nhead = len(iheadList)
         
        for ih in range(nhead):
            iheadx = iheadList[ih].split("/")[0].split("=")
            try:
                ihdr[iheadx[0]] = float(iheadx[1])
            except:
                ihdr[iheadx[0]] = str(iheadx[1])

        ihdr["CTYPE1"] = "RA---TPV"
        ihdr["CTYPE2"] = "DEC--TPV"
        ihdr["RADESYS"] = "ICRS"
        ihdr["CUNIT1"]  = "deg"
        ihdr["CUNIT2"]  = "deg"
        fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True)

        os.system("rm %s"%ihead)
        os.system("mv *pdf %s"%(figdir))
 
       # generate the final catalog
        loggerpho.info("############ Generate final SExtractor catalog")
        isexComd = sexComdPSF%(iobjFrame_o0,sexConf,sexParam2,ildacn,saturate_value,ipsfFrameNew)
        os.system(isexComd)

        ## some diagnosis charts and files
        
        #logger.info(np.array(ildac["MAG_AUTO"])[1]) 
        hdul=fits.open(ildacn)
        hdul[2].data=ildac
        #logger.info(np.array(hdul[1].data[0][0][5]),math.log10(exptime/10)*2.5,exptime)
        
        ildac["FLUX_AUTO"]=ildac["FLUX_AUTO"]/float(exptime)
        ildac["FLUX_PSF"]=ildac["FLUX_PSF"]/float(exptime)
        ildac["FLUX_APER"]=ildac["FLUX_APER"]/float(exptime)

        ildac["FLUXERR_AUTO"]=ildac["FLUXERR_AUTO"]/float(exptime)
        ildac["FLUXERR_PSF"]=ildac["FLUXERR_PSF"]/float(exptime)
        ildac["FLUXERR_APER"]=ildac["FLUXERR_APER"]/float(exptime)
         
        ildac["MAG_APER"]=ildac["MAG_APER"]+math.log10(exptime)*2.5
        ildac["MAG_AUTO"]=ildac["MAG_AUTO"]+ math.log10(exptime)*2.5
        ildac["MAG_PSF"] = ildac["MAG_PSF"]+math.log10(exptime)*2.5

        ildac["MAGERR_AUTO"] = 1.0857 * (ildac["FLUXERR_AUTO"]/ildac["FLUX_AUTO"])
        ildac["MAGERR_PSF"] = 1.0857 * (ildac["FLUXERR_PSF"]/ildac["FLUX_PSF"])
        ildac["MAGERR_APER"] = 1.0857 * (ildac["FLUXERR_APER"]/ildac["FLUX_APER"])
        
 
        hdul.writeto(ildacn,overwrite=True)
        hdul.close()

        ildac        = fits.getdata(ildacn,ext=2)
        iflag        = ildac["FLAGS"]
        isnr         = ildac["SNR_WIN"]
        iximg        = ildac["X_IMAGE"]
        iyimg        = ildac["Y_IMAGE"]
        ira          = ildac["ALPHA_J2000"]
        idec         = ildac["DELTA_J2000"]
        imagAUTO     = ildac["MAG_AUTO"] 
        imagPSF      = ildac["MAG_PSF"] 
        imagErrPSF   = ildac["MAGERR_PSF"]
        ifluxAUTO    = ildac["FLUX_AUTO"]
        ifluxErrAUTO = ildac["FLUXERR_AUTO"]    
        ifluxPSF     = ildac["FLUX_PSF"]
        ifluxErrPSF  = ildac["FLUXERR_PSF"]
        ifwhm        = ildac["FWHM_IMAGE"]
        iaimg        = ildac["AWIN_IMAGE"]
        ibimg        = ildac["BWIN_IMAGE"]
        itheta       = ildac["THETAWIN_IMAGE"]
        iflux_radius= ildac["FLUX_RADIUS"]

         
        # 1) region file
        iregn  = iobjFrame_o0[:-4] + "reg"
        wds9reg(ira,idec,flag=None,radius=8.0,unit="arcsec",color="green",outfile=iregn)
        psf_id    = iflux_radius>1.5
        # select high-snr stars
        gid        = (iflag==0) & (isnr>30.0)
        ira        = ira[gid]
        idec       = idec[gid]
        imagAUTO   = imagAUTO[gid]
        imagPSF    = imagPSF[gid]
        imagErrPSF = imagErrPSF[gid]
        ifwhm      = ifwhm[gid]
        iximg      = iximg[gid]
        iyimg      = iyimg[gid]
        iaimg      = iaimg[gid]
        ibimg      = ibimg[gid]
        itheta     = itheta[gid]

         

        idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)

        # 2) determine the astrometric accuracy
        idra  = 3600.0*(ira[ids]-raRef[idRef])
        iddec = 3600.0*(idec[ids]-decRef[idRef])
        iraSig, idecSig = mad(idra), mad(iddec)
         
        loggerpho.info("Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))
        
        #  ##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& new fy20230131    
        # if(len(idRef))<10:
         
        #     scamp_ldac=reget_refsed(ancdir,objID,ira,idec)
        #     loggerpho.info("   try online refsed to refine the astrometry")
         
        #     ildac=fits.getdata(ildacn,ext=2)
        #     scamp_ldac=reget_refsed(ancdir,objID,ira,idec)
        #     if(filtID =='v' or filtID =='u'):
        #         ixcenloc=np.where((ildac["X_IMAGE"]>c/2+10) | (ildac["X_IMAGE"]<c/2-10))
        #         #ixyloc=np.where((iximg<c-(1000/binfact))&(iximg>1000/binfact))#&(iyimg<6100)&(iyimg>30))

        #         ildac =ildac[ixcenloc]
        #     hdul=fits.open(ildacn)###
        #     hdul[1].data=ildac

        #     reg(ildac,ildac_select_reg)
        
        #     hdul.close()
                
                
        #     iscampComd = scampComd%(ildacn,scampConf,scamp_ldac)
        #     os.system(iscampComd)
        #     ihead   = ildacn[:-4] + "head"
        #     ihead1 = ildacn[:-4] + "ahead"
                

        #     os.system("mv %s %s"%(ihead,ihead1))
        #     iscampComd2 = scampComd%(ildacn,scampConf2,scamp_ldac)
        #     os.system(iscampComd2)


        #     ihead   = ildacn[:-4] + "head"
        #     iimgMat, ihdr = fits.getdata(iobjFrame_o0,header=True)
        #     iheadList = open(ihead,"r").read().splitlines()[3:-1]
        #     nhead = len(iheadList)
            
        #     for ih in range(nhead):
        #         iheadx = iheadList[ih].split("/")[0].split("=")
        #         try:
        #             ihdr[iheadx[0]] = float(iheadx[1])
        #         except:
        #             ihdr[iheadx[0]] = str(iheadx[1])

        #     ihdr["CTYPE1"] = "RA---TPV"
        #     ihdr["CTYPE2"] = "DEC--TPV"
        #     ihdr["RADESYS"] = "ICRS"
        #     ihdr["CUNIT1"]  = "deg"
        #     ihdr["CUNIT2"]  = "deg"
        #     fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True)

        #     os.system("rm %s"%ihead)
        #     os.system("mv *pdf %s"%(figdir))

    
            
            
        # # generate the final catalog
        #     loggerpho.info("############ Generate final SExtractor catalog")
        #     isexComd = sexComdPSF%(iobjFrame_o0,sexConf,sexParam2,ildacn,saturate_value,ipsfFrameNew)
        #     os.system(isexComd)

        #     ## some diagnosis charts and files
        #     ildac        = fits.getdata(ildacn,ext=2)
        #     iflag        = ildac["FLAGS"]
        #     isnr         = ildac["SNR_WIN"]
        #     iximg        = ildac["X_IMAGE"]
        #     iyimg        = ildac["Y_IMAGE"]
        #     ira          = ildac["ALPHA_J2000"]
        #     idec         = ildac["DELTA_J2000"]
        #     imagAUTO     = ildac["MAG_AUTO"]#+math.log10(exptime)*2.5
        #     imagPSF      = ildac["MAG_PSF"]#+math.log10(exptime)*2.5
        #     imagErrPSF   = ildac["MAGERR_PSF"]
        #     ifluxAUTO    = ildac["FLUX_AUTO"]
        #     ifluxErrAUTO = ildac["FLUXERR_AUTO"]
        #     ifluxPSF     = ildac["FLUX_PSF"]
        #     ifluxErrPSF  = ildac["FLUXERR_PSF"]
        #     ifwhm        = ildac["FWHM_IMAGE"]
        #     iaimg        = ildac["AWIN_IMAGE"]
        #     ibimg        = ildac["BWIN_IMAGE"]
        #     itheta       = ildac["THETAWIN_IMAGE"]
        #     iflux_radius= ildac["FLUX_RADIUS"]
        #     #logger.info(np.array(ildac["MAG_AUTO"])[1]) 
        #     hdul=fits.open(ildacn)
        #     hdul[2].data=ildac
        #     #logger.info(np.array(hdul[1].data[0][0][5]),math.log10(exptime/10)*2.5,exptime)
        #     ildac["MAG_AUTO"]=ildac["MAG_AUTO"]#+math.log10(exptime)*2.5
        #     ildac["MAG_PSF"] = ildac["MAG_PSF"]#+math.log10(exptime)*2.5
        #     ildac["FLUX_AUTO"]=ildac["FLUX_AUTO"]#/float(exptime)
        #     ildac["FLUX_PSF"]=ildac["FLUX_PSF"]#/float(exptime)
        #     hdul.writeto(ildacn,overwrite=True)
        #     hdul.close()

            
        #     # 1) region file
        #     iregn  = iobjFrame_o0[:-4] + "reg"
        #     wds9reg(ira,idec,flag=None,radius=8.0,unit="arcsec",color="green",outfile=iregn)
        #     psf_id    = iflux_radius>1.5
        #     # select high-snr stars
        #     gid        = (iflag==0) & (isnr>30.0)
        #     ira        = ira[gid]
        #     idec       = idec[gid]
        #     imagAUTO   = imagAUTO[gid]
        #     imagPSF    = imagPSF[gid]
        #     imagErrPSF = imagErrPSF[gid]
        #     ifwhm      = ifwhm[gid]
        #     iximg      = iximg[gid]
        #     iyimg      = iyimg[gid]
        #     iaimg      = iaimg[gid]
        #     ibimg      = ibimg[gid]
        #     itheta     = itheta[gid]

            

        #     idRef, ids = crossmatch(raRef, decRef, ira, idec, aperture=1.5)

        #     # 2) determine the astrometric accuracy
        #     idra  = 3600.0*(ira[ids]-raRef[idRef])
        #     iddec = 3600.0*(idec[ids]-decRef[idRef])
        #     iraSig, idecSig = mad(idra), mad(iddec)
            
        #     loggerpho.info("Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))
        # ##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& new end fy20230131            
            
            
        
            

        # 3) determine the zeropoint
        try:
            # izpfig = figdir + "zpoint_%s.pdf"%(ildacn.split("/")[-1][:-5])
            # ipsffig = figdir + "psfPhot_%s.pdf"%(ildacn.split("/")[-1][:-5])
            
            izpArr = magRef[idRef] - imagPSF[ids]
            imask  = sigma_clip(izpArr, sigma=3.0, maxiters=3, masked=True)
            izpArr = izpArr[~imask.mask]
            izp    = np.median(izpArr)
            izpErr = mad(izpArr)
       
            loggerpho.info("    Initial zeropoint is %6.3f mag determined by %4d stars"%(izp, len(izpArr)))
            if(len(izpArr)<10):
               nameonj=iobjFrame_o0.rsplit('/',1)[1]
               loggerpho.info('XXXXXXXXX'+nameonj+':  the astrometric measurement may not success, please check all logger.info information and source fits!')
            # plotFun.zpPlot(izpArr, imagPSF[ids][~imask.mask], figout=izpfig)
            # bkgNoise = np.sqrt(np.pi)*np.median(ifwhm) * (np.min([ihdr["SKYRMS1"],ihdr["SKYRMS2"]]))
            # plotFun.mag2Err(ifluxPSF[psf_id],ifluxErrPSF[psf_id],ifluxAUTO[psf_id],ifluxErrAUTO[psf_id],bkgNoise=bkgNoise,zpoint=izp,gain=ihdr["GAIN"],figout=ipsffig)
            # 4) show the Whisker plot
        # try:
        #     iwskfig     = figdir + "ewhiskerEll_%s.pdf"%(ildacn.split("/")[-1][:-5])
        #     iwskfwhmfig = figdir + "ewhiskerFWHM_%s.pdf"%(ildacn.split("/")[-1][:-5])
        #     plotFun.ewhisker(iximg, iyimg, iaimg, ibimg, itheta, figout=iwskfig)
        #     #plotFun.ewhiskerEllipse(iximg, iyimg, iaimg, ibimg, itheta, scale=30.0, figout=iwskfig)
        #     plotFun.ewhiskerFWHM(iximg,iyimg,ifwhm,itheta,scale=20,figout=iwskfwhmfig)
        # except Exception as e:
        #     logger.info(str(e))
        #     logger.info('the Whisker figure can not be generated!') 
        #     logger.info('XXXXXXXXXXXXXXXXX the Whisker figure can not be generated!')      
        #     continue    

            # 5) update image header
            # ihdr["ZPMAG"] = izp
            # ihdr["ZPERR"] = izpErr
            # fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True)


            # 6）split ccd and update the zpmaag and zperr

            dateobs,exptime=ihdr["DATE-OBS"],ihdr["EXPTIME"]
            ora=ihdr['OBJCTRA']
            odec=ihdr['OBJCTDEC']
            ldacdata=fits.open(ildacn)
            ldacdata.writeto(iobjFrame_o0[:-5] + '_'+date+"_sexcat.fits",overwrite=True)
            
        
            #hdrlist.append([dateobs,exptime])
            ldacList.write(ildacn+"\n")
            subImgList.write(iobjFrame_o0+"\n")

            #os.system('rm %s'%(iobjFrameNew))
            #os.system('rm %s'%(iobjFrameNew[-4]+'ldac'))
            filename=iobjFrame_o0.rsplit('/',1)[1]
            loggerpho.info(filename+' has been successfully processed!')
        except Exception as e:
            loggerpho.info(str(e))
            loggerpho.info('the zpPlot figure can not be generated!')         
        finally:
            for k in list(locals().keys()):
                # if locals[k] is np.nan:
                try:
                    del locals[k]
                except:
                    continue
            gc.collect()
          
        
    ldacList.close()
    subImgList.close()
    
    
    loggerpho.info("^_^ ALL astrometry and photometry processes are ok!")
@elapsed_time
def pro(tid_target_filter,date):
    rootpath=get_current_dir()
    photometry(rootpath,tid_target_filter,date)
    logger.info('the astrometric process is done')

# tid_target_filter=sys.argv[1]
# date = sys.argv[2]
# pro(tid_target_filter,date)

if __name__ == "__main__":
    import argparse
    import time
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20230129', help='输入处理日期')
    parser.add_argument('--tid_target_filter', type=str, help='')
    args = parser.parse_args()  # 参数解析像·
    loggerloguru = Logger(args.date, '', "_ycom").get_logger
    starttime=time.time()
    pro(args.tid_target_filter,args.date)
    endtime=time.time()
    loggerloguru.info("[_yphostcom1] "+"time spending:{}s".format(endtime-starttime))
