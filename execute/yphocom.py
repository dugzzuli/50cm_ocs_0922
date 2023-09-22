import time
from execute.flatdiff import mkdir, readfits
from lib.phot.ygn import savenpy
import numpy as np
from astropy.io import fits
from astropy import wcs
import os, sys
from collections import Counter
from scipy.spatial import cKDTree as ckdt
from astropy.stats import sigma_clip
from astropy.io import fits
import numpy as np
import os, sys
from scipy.optimize import leastsq
#import matplotlib
#import pylab as pl
from astropy import wcs
from astropy import coordinates as coord, units as u
import datetime
from astropy.time import Time
import ccdproc
from astropy.nddata import CCDData
from astropy.table import Table
from loguru import logger
from lib.phot.yrefsed import refSed_gaia 
import lib.phot.PlotFun as plotFun
from lib.phot.ybias import get_current_dir
import warnings

#from execute yphocom1 import save_npy  
warnings.filterwarnings('error') 





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


def update_scamp_head(fits_filename,head_filename):
    ihead  = head_filename
    headList = open(ihead,"r").read().splitlines()
    try:
        scidhr,sciimgMat= readfits(fits_filename)
    except:
        logger.info('the fits file is broken',fits_filename)
    dic=get_scamp_head(ihead)
    scidhr["CTYPE1"]='RA---TPV'
    scidhr["CTYPE2"]='DEC--TPV' 
    scidhr["CRVAL1"]=float(dic["CRVAL1"])
    scidhr["CRVAL2"]=float(dic["CRVAL2"])
    scidhr["CRPIX1"]=float(dic["CRPIX1"])
    scidhr["CRPIX2"]=float(dic["CRPIX2"])
    scidhr["CD1_1"]=float(dic["CD1_1"])
    scidhr["CD1_2"]=float(dic["CD1_2"])
    scidhr["CD2_1"]=float(dic["CD2_1"])
    scidhr["CD2_2"]=float(dic["CD2_2"])
    scidhr["PV1_0"]=float(dic["PV1_0"])
    scidhr["PV1_1"]=float(dic["PV1_1"])
    scidhr["PV1_2"]=float(dic["PV1_2"])
    scidhr["PV1_4"]=float(dic["PV1_4"])
    scidhr["PV1_5"]=float(dic["PV1_5"])
    scidhr["PV1_6"]=float(dic["PV1_6"])
    scidhr["PV1_7"]=float(dic["PV1_7"])
    scidhr["PV1_8"]=float(dic["PV1_8"])
    scidhr["PV1_9"]=float(dic["PV1_9"])
    scidhr["PV1_10"]=float(dic["PV1_10"])
    scidhr["PV2_0"]=float(dic["PV2_0"])
    scidhr["PV2_1"]=float(dic["PV2_1"])
    scidhr["PV2_2"]=float(dic["PV2_2"])
    scidhr["PV2_4"]=float(dic["PV2_4"])
    scidhr["PV2_5"]=float(dic["PV2_5"])
    scidhr["PV2_6"]=float(dic["PV2_6"])
    scidhr["PV2_7"]=float(dic["PV2_7"])
    scidhr["PV2_8"]=float(dic["PV2_8"])
    scidhr["PV2_9"]=float(dic["PV2_9"])
    scidhr["PV2_10"]=float(dic["PV2_10"])
    scidhr["FLXSCALE"]=1.0
    return scidhr,sciimgMat

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
        sra = np.array(ra.split(":"), dtype=float)
        sdec = np.array(dec.split(":"), dtype=float)
        sra = ((sra[-1]/60.0+sra[1])/60.0 + sra[0])*15.0
        if sdec[0]<0.0:
            sdec = -((sdec[-1]/60.0+sdec[1])/60.0 + abs(sdec[0]))
        else:
            sdec = (sdec[-1]/60.0+sdec[1])/60.0 + sdec[0]
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



def match_pos(iobjFrameNew,ancdir,objID):
    ildacn = iobjFrameNew[:-4] + "ldac"
    # hdr=fits.getheader(iobjFrameNew)
    # saturate_value=hdr['SATURATE']
    # isexComd = sexComd%(iobjFrameNew,confdir,ildacn,saturate_value)
    # os.system(isexComd)
    try:
        ildac = fits.getdata(ildacn,ext=2)
        iximg, iyimg = ildac["X_IMAGE"], ildac["Y_IMAGE"]
        ixMin, ixMax = iximg.min(), iximg.max()
        iyMin, iyMax = iyimg.min(), iyimg.max()
    except:
        logger.info('the fits file is broken',iobjFrameNew)
        return 

    iimgMat, ihdr = fits.getdata(iobjFrameNew,header=True)
    iwcs     = wcs.WCS(ihdr)
    irefTMPn = iobjFrameNew[:-5] + "_ref.cat"
    irefTMP  = open(irefTMPn,"w")
    refCatn  = ancdir + "GaiaStar_%s.ldac"%objID
    refCat2n = ancdir + "GaiaStar_%s.cat"%objID
    refCat2  = Table.read(refCat2n, format="ascii")
    raRef    = refCat2["ra"]
    decRef   = refCat2["dec"]
    nrefStar = len(raRef)
    # refine the WCS
    iflag, imag = ildac["FLAGS"], ildac["MAG_AUTO"]
    igid = iflag==0
    iximgG, iyimgG, imagG = iximg[igid], iyimg[igid], imag[igid]
    isciTMPn = iobjFrameNew[:-5] + "_sci.cat"
    isciTMP  = open(isciTMPn,"w")
    for k in range(len(imagG)):
        isciTMP.write("%4d %9.3f %9.3f %6.2f\n"%(k+1,iximgG[k],iyimgG[k],imagG[k]))
    isciTMP.close()

    # pixel coordinates of the reference stars
    
    
    for k in range(nrefStar):
        kimgpos = iwcs.all_world2pix(raRef[k],decRef[k],1)
        kxpos   = kimgpos[0].flatten()[0]
        kypos   = kimgpos[1].flatten()[0]
        if kxpos<ixMin or kxpos>ixMax: continue
        if kypos<iyMin or kypos>iyMax: continue
        irefTMP.write("%4d %9.3f %9.3f %6.2f\n"%(k+1,kxpos,kypos,magRef[k]))
    irefTMP.close()

    try:
    # match two star catalogs
        imatchFile = iobjFrameNew[:-4] + "_match"
        imatch = "match %s 1 2 3 %s 1 2 3 nobj=150 outfile=%s> %s"%(isciTMPn,irefTMPn,iobjFrameNew[:-5],imatchFile)
        os.system(imatch)
        imatchRes = open(imatchFile).read().splitlines()[0]
        imatchRes = imatchRes.split()[1:]
        ixOffset  = float(imatchRes[0].split("=")[-1])
        iyOffset  = float(imatchRes[3].split("=")[-1])
        irot1     = float(imatchRes[1].split("=")[-1])
        irot2     = float(imatchRes[2].split("=")[-1])
        irot      = np.arctan(irot2/irot1)*180.0/np.pi

        ixScale   = float(imatchRes[1].split("=")[-1])
        iyScale   = float(imatchRes[5].split("=")[-1])
        ixFlip    = np.sign(ixScale)
        iyFlip    = np.sign(iyScale)
        if ixFlip!=1 or iyFlip!=1:
            ihdr["CD1_1"] = ihdr["CD1_1"] * ixFlip
            ihdr["CD2_2"] = ihdr["CD2_2"] * ixFlip
        
            iwcs = wcs.WCS(ihdr)
            irefTMPn = iobjFrameNew[:-5] + "_ref.cat"
            irefTMP  = open(irefTMPn,"w")
            for k in range(nrefStar):
                kimgpos = iwcs.all_world2pix(raRef[k],decRef[k],1)
                kxpos   = kimgpos[0].flatten()[0]
                kypos   = kimgpos[1].flatten()[0]
                if kxpos<0 or kxpos>ximg: continue
                if kypos<0 or kypos>yimg: continue
                irefTMP.write("%4d %9.3f %9.3f %6.2f\n"%(k+1,kxpos,kypos,magRef[k]))
            irefTMP.close()
    #match two star catalogs
            imatchFile = iobjFrameNew[:-4] + "_match"
            imatch = "match %s 1 2 3 %s 1 2 3 outfile=%s nobj=150 > %s"%(isciTMPn,irefTMPn,iobjFrameNew[:-5],imatchFile)
            os.system(imatch)
            imatchRes = open(imatchFile).read().splitlines()[0]
            imatchRes = imatchRes.split()[1:]
            ixOffset  = float(imatchRes[0].split("=")[-1])
            iyOffset  = float(imatchRes[3].split("=")[-1]) 
        
            irot1     = float(imatchRes[1].split("=")[-1])      
            irot2     = float(imatchRes[2].split("=")[-1])
            irot      = np.arctan(irot2/irot1)*180.0/np.pi
        logger.info(" Refine result: (x_offset, y_offset, rotation) = (%4d, %4d, %.3f)"%(int(ixOffset),int(iyOffset),irot))
        
        os.system("rm %s.mtA %s.mtB %s.unA %s.unB"%(iobjFrameNew[:-5],iobjFrameNew[:-5],iobjFrameNew[:-5],iobjFrameNew[:-5]))
        os.system("rm %s %s %s"%(imatchFile,irefTMPn,isciTMPn))

    #update the header
        ihdr["CRPIX1"] = ihdr["CRPIX1"] - ixOffset
        ihdr["CRPIX2"] = ihdr["CRPIX2"] - iyOffset
    except:
        logger.info(iobjFrameNew+':  match process is failed')
    fits.writeto(iobjFrameNew,iimgMat,ihdr,overwrite=True) 


def ldac_gen(iobjFrame_o0):
    # 3S configuration
    sexParam1  = confdir + "default.param"
    sexParam2  = confdir + "default.paramNew"
    sexConf    = confdir + "default.sex"
    swarpConf  = confdir + "default.swarp"
    scampConf  = confdir + "default.scamp"
    psfexConf  = confdir + "default.psfex"

    sexComd1   = "sex %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s "
    sexComd2   = "-DETECT_MINAREA 3 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 "
    sexComd3   = "-SATUR_LEVEL %8.1f -CHECKIMAGE_TYPE NONE "
    sexComd4   = "-PSF_NAME %s "
    sexComd    = sexComd1 + sexComd2 + sexComd3
    sexComdPSF = sexComd + sexComd4

    scampComd1 = "scamp %s -c %s -ASTREFCAT_NAME %s -MOSAIC_TYPE LOOSE "
    scampComd2 = "-FWHM_THRESHOLDS 2,20 -SN_THRESHOLDS 20,1000 "
    scampComd  = scampComd1 + scampComd2

    swarpComd1 = "swarp @%s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s "
    swarpComd2 = "-CENTER_TYPE ALL -PIXEL_SCALE %3.1f -RESAMPLE_DIR %s "
    swarpComd  = swarpComd1 + swarpComd2

    psfexComd  = "psfex %s -c %s -PSF_SIZE %d,%d -PSF_DIR %s "
    



    try:
        ihdr,iimgMat= readfits(iobjFrameNew)
    except:
        logger.info('the fits file is broken',iobjFrameNew)
   # run sextractor to extract objects
    ildacn = iobjFrame_o0[:-4] + "ldac"
    saturate_value=ihdr['SATURATE']
    isexComd = sexComd%(iobjFrame_o0,sexConf,sexParam1,ildacn,saturate_value)
    os.system(isexComd)

    # construct PSF model
    ipsfComd = psfexComd%(ildacn, psfexConf, 21, 21, tscidir)
    os.system(ipsfComd)
    os.system("mv *pdf %s"%(figdir))
    os.system("mv *fits %s"%tscidir)
    try:
        ildac = fits.getdata(ildacn,ext=2)
        iximg, iyimg = ildac["X_IMAGE"], ildac["Y_IMAGE"]
        ixMin, ixMax = iximg.min(), iximg.max()
        iyMin, iyMax = iyimg.min(), iyimg.max()
    except:
        logger.info('the fits file is broken',iobjFrame_o0)

    # refine the WCS
    iflag, imag = ildac["FLAGS"], ildac["MAG_AUTO"]
    igid = iflag==0
    iximgG, iyimgG, imagG = iximg[igid], iyimg[igid], imag[igid]
    isciTMPn = iobjFrame_o0[:-5] + "_sci.cat"
    isciTMP  = open(isciTMPn,"w")
    for k in range(len(imagG)):
        isciTMP.write("%4d %9.3f %9.3f %6.2f\n"%(k+1,iximgG[k],iyimgG[k],imagG[k]))
    isciTMP.close()


  #iimgMat, ihdr = fits.getdata(iobjFrameNew,header=True)
    iwcs     = wcs.WCS(ihdr)
    irefTMPn = iobjFrame_o0[:-5] + "_ref.cat"
    irefTMP  = open(irefTMPn,"w")
    for k in range(nrefStar):
        kimgpos = iwcs.all_world2pix(raRef[k],decRef[k],1)
        kxpos   = kimgpos[0].flatten()[0]
        kypos   = kimgpos[1].flatten()[0]
        if kxpos<ixMin or kxpos>ixMax: continue
        if kypos<iyMin or kypos>iyMax: continue
        irefTMP.write("%4d %9.3f %9.3f %6.2f\n"%(k+1,kxpos,kypos,magRef[k]))
    irefTMP.close()

  


def photometry(rootpath,ttfname,date):
    upath=rootpath+'reception/'+str(date)+'/'
    scipath = upath+'sci/'+ttfname+'/'
    #logger.info(lpathsci[0:5])
    lpathsubbkg=np.load(scipath+'fileguide/'+ttfname+'_subbkg_o0.npy')
    lpathsubbkg=np.sort(lpathsubbkg)
    #logger.info('the preocessed fits are:', lpathsubbkg)
    ################lpathsubbkg=lpathsubbkg[102:119]   
    #basic setting
    #dirname, filename = os.path.split(os.path.abspath(__file__))
    imgdir  = rootpath + "images/"
    confdir = rootpath + "config/"
    ancdir  = rootpath + "ancData/"
    figdir  = rootpath + "figures/"+date+"_"+ttfname+"/"
    pdfdir = rootpath+'50cmpy/'
    mkdir(figdir)   
    ttf=ttfname.split('_')
    tid=ttf[0]
    targetid=ttf[1]
    try:
        hrs,imgs=readfits(lpathsubbkg[0])
    except:
        logger.info("Error: 没有找到文件或读取文件失败")
        hrs,imgs=readfits(lpathsubbkg[1])
    


    try:
            ira  = ":".join(hrs["OBJCTRA"].split())
            idec = ":".join(hrs["OBJCTDEC"].split())
            ira, idec = d2hms(ira, idec, conv=1)
    except:
        logger.info('these are no parameters of object ra and dec in header of this fits')
    
     
    r,c=np.shape(imgs)
    icountmax=np.nanmax(imgs)
     
    telsID  = ttf[0]
    filtID  = ttf[2][-1]
    objID   = ttf[1]
    obsDate = date
    ximg    = c
    yimg    = r
    binfact = int(hrs['XBINNING'])
    pscale  = 0.297 * binfact # pixel scale in unit of arcsec
    saturate_value=hrs['SATURATE']
   
    refCatn  = ancdir + "GaiaStar_%s.ldac"%objID
    refCat2n = ancdir + "GaiaStar_%s.cat"%objID
    if not (os.path.exists(refCatn)):
        refSed_gaia(rootpath,ttfname,date,ira,idec)
    refCat2  = Table.read(refCat2n, format="ascii")
    raRef    = refCat2["ra"]
    decRef   = refCat2["dec"]
    magRef   = refCat2["mag%s"%filtID.upper()]
    nrefStar = len(raRef)

    scipath=rootpath+'reception/'+str(date)+'/sci/'+ttfname+'/'
    tscidir = scipath 
    mkdir(tscidir)
  

    # 3S configuration
    sexParam1  = confdir + "default.param"
    sexParam2  = confdir + "default.paramNew"
    sexConf    = confdir + "default.sex"
    swarpConf  = confdir + "default.swarp"
    scampConf  = confdir + "default.scamp"
    psfexConf  = confdir + "default.psfex"

    sexComd1   = "sex %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s "
    sexComd2   = "-DETECT_MINAREA 3 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 "
    sexComd3   = "-SATUR_LEVEL %8.1f -CHECKIMAGE_TYPE NONE "
    sexComd4   = "-PSF_NAME %s "
    sexComd    = sexComd1 + sexComd2 + sexComd3
    sexComdPSF = sexComd + sexComd4

    scampComd1 = "scamp %s -c %s -ASTREFCAT_NAME %s -MOSAIC_TYPE LOOSE "
    scampComd2 = "-FWHM_THRESHOLDS 2,20 -SN_THRESHOLDS 20,1000 "
    scampComd  = scampComd1 + scampComd2

    swarpComd1 = "swarp @%s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s "
    swarpComd2 = "-CENTER_TYPE ALL -PIXEL_SCALE %3.1f -RESAMPLE_DIR %s "
    swarpComd  = swarpComd1 + swarpComd2

    psfexComd  = "psfex %s -c %s -PSF_SIZE %d,%d -PSF_DIR %s "
    
    
    t0 = time.time()
    nobj = len(lpathsubbkg)
    #fileguide record
    # hdr=fits.getheader(iobjFrameNew)
    # saturate_value=hdr['SATURATE']
    # isexComd = sexComd%(iobjFrameNew,confdir,ildacn,saturate_value)
    # os.system(isexComd)
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
     
    hdrlist=[]
    for i in range(nobj):
        logger.info('####### the preocessed fits is: ', lpathsubbkg[i])
        #对每一个frame run sextractor
        iobjFrameNew = lpathsubbkg[i]#[:-3] + "sub.fit"
        iobjFrame_o0=iobjFrameNew[:-14]+date+'_subbkg_o0.fits'
        ipsfFrameNew = iobjFrame_o0[:-4] + "psf"
        try:
            ihdr,iimgMat= readfits(iobjFrameNew)
        except:
            logger.info('the fits file is broken',iobjFrameNew)
            continue

        binfact=int(ihdr['XBINNING'])
        pscale  = 0.297 * binfact # pixel scale in unit of arcsec
        
        try:
            ira  = ":".join(ihdr["OBJCTRA"].split())
            idec = ":".join(ihdr["OBJCTDEC"].split())
            ira, idec = d2hms(ira, idec, conv=1)
        except:
            logger.info('these are no parameters of object ra and dec in header of this fits')
            continue
 
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
        ihdr["CD1_1"]   = -pscale/3600.0
        ihdr["CD1_2"]   = 0.0
        ihdr["CD2_1"]   = 0.0
        ihdr["CD2_2"]   = -pscale/3600.0
        fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True)
 
        # run sextractor to extract objects
        ildacn = iobjFrame_o0[:-4] + "ldac"
        saturate_value=ihdr['SATURATE']
        isexComd = sexComd%(iobjFrame_o0,sexConf,sexParam1,ildacn,saturate_value)
        os.system(isexComd)

        # construct PSF model
        ipsfComd = psfexComd%(ildacn, psfexConf, 21, 21, tscidir)
        os.system(ipsfComd)
        os.system("mv *pdf %s"%(figdir))
        os.system("mv *fits %s"%tscidir)
        try:
            ildac = fits.getdata(ildacn,ext=2)
            iximg, iyimg = ildac["X_IMAGE"], ildac["Y_IMAGE"]
            ixMin, ixMax = iximg.min(), iximg.max()
            iyMin, iyMax = iyimg.min(), iyimg.max()
        except:
            logger.info('the fits file is broken',iobjFrame_o0)
            continue
        
        # refine the WCS
        iflag, imag = ildac["FLAGS"], ildac["MAG_AUTO"]
        igid = iflag==0
        iximgG, iyimgG, imagG = iximg[igid], iyimg[igid], imag[igid]
        isciTMPn = iobjFrame_o0[:-5] + "_sci.cat"
        isciTMP  = open(isciTMPn,"w")
        for k in range(len(imagG)):
            isciTMP.write("%4d %9.3f %9.3f %6.2f\n"%(k+1,iximgG[k],iyimgG[k],imagG[k]))
        isciTMP.close()

        # pixel coordinates of the reference stars
        
        #iimgMat, ihdr = fits.getdata(iobjFrameNew,header=True)
        iwcs     = wcs.WCS(ihdr)
        irefTMPn = iobjFrame_o0[:-5] + "_ref.cat"
        irefTMP  = open(irefTMPn,"w")
        for k in range(nrefStar):
            kimgpos = iwcs.all_world2pix(raRef[k],decRef[k],1)
            kxpos   = kimgpos[0].flatten()[0]
            kypos   = kimgpos[1].flatten()[0]
            if kxpos<ixMin or kxpos>ixMax: continue
            if kypos<iyMin or kypos>iyMax: continue
            irefTMP.write("%4d %9.3f %9.3f %6.2f\n"%(k+1,kxpos,kypos,magRef[k]))
        irefTMP.close()
        
        logger.info("   Run Scamp to refine the astrometry")
        try:
            try:
                iscampComd = scampComd%(ildacn,scampConf,refCatn)
                os.system(iscampComd)    
                import warnings
            except warnings as e:
                logger.info('the scamp1 is failed,we will try match process')
                match_pos(iobjFrame_o0,ancdir,objID)
                ldac_gen(iobjFrame_o0)
                iscampComd = scampComd%(ildacn,scampConf,refCatn)
        except:
            logger.info('astrometry is failed!#####'+ iobjFrame_o0)
            continue
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
        logger.info("    Generate final SExtractor catalog")
        isexComd = sexComdPSF%(iobjFrame_o0,sexConf,sexParam2,ildacn,saturate_value,ipsfFrameNew)
        os.system(isexComd)    

        ## some diagnosis charts and files
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
        gid        = (iflag==0) & (isnr>50.0) 
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
        logger.info("    Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))
        ifigname = figdir + "astrometry_%s.pdf"%(ildacn.split("/")[-1][:-5])
        plotFun.astrometryPlot(idra, iddec, sigma=(iraSig,idecSig), figout=ifigname)

        # 3) determine the zeropoint
        izpfig = figdir + "zpoint_%s.pdf"%(ildacn.split("/")[-1][:-5])
        ipsffig = figdir + "psfPhot_%s.pdf"%(ildacn.split("/")[-1][:-5])
        izpArr = magRef[idRef] - imagPSF[ids]
        imask  = sigma_clip(izpArr, sigma=3.0, maxiters=3, masked=True)
        izpArr = izpArr[~imask.mask]
        izp    = np.median(izpArr)
        izpErr = mad(izpArr)
        plotFun.zpPlot(izpArr, imagPSF[ids][~imask.mask], figout=izpfig)
        logger.info("    Initial zeropoint is %6.3f mag determined by %4d stars"%(izp, len(izpArr)))
        
        bkgNoise = np.sqrt(np.pi)*np.median(ifwhm) * (np.min([ihdr["SKYRMS1"],ihdr["SKYRMS2"]]))
        plotFun.mag2Err(ifluxPSF[psf_id],ifluxErrPSF[psf_id],ifluxAUTO[psf_id],ifluxErrAUTO[psf_id],bkgNoise=bkgNoise,zpoint=izp,gain=ihdr["GAIN"],figout=ipsffig)
        
        # 4) show the Whisker plot
        iwskfig     = figdir + "ewhiskerEll_%s.pdf"%(ildacn.split("/")[-1][:-5])
        iwskfwhmfig = figdir + "ewhiskerFWHM_%s.pdf"%(ildacn.split("/")[-1][:-5])
        plotFun.ewhisker(iximg, iyimg, iaimg, ibimg, itheta, figout=iwskfig)
        #plotFun.ewhiskerEllipse(iximg, iyimg, iaimg, ibimg, itheta, scale=30.0, figout=iwskfig)
        plotFun.ewhiskerFWHM(iximg,iyimg,ifwhm,itheta,scale=20,figout=iwskfwhmfig)  
    
        # 5) update image header
        ihdr["ZPMAG"] = izp
        ihdr["ZPERR"] = izpErr
        fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True) 
         

        # 6）split ccd and update the zpmaag and zperr
        
        dateobs,exptime=ihdr["DATE-OBS"],ihdr["EXPTIME"]
        ora=ihdr['OBJCTRA']
        odec=ihdr['OBJCTDEC']
        ldacdata=fits.open(ildacn)
        ldacdata.writeto(iobjFrame_o0[:-5] + '_'+date+"_sexcat.fits",overwrite=True) 

        sub_osubimg1=iobjFrame_o0[:-6]+'1.fits'
        sub_osubimg2=iobjFrame_o0[:-6]+'2.fits'
        if os.path.exists(sub_osubimg1): 
            os.remove(sub_osubimg1)
        os.system("imcopy %s[%d:%d,%d:%d] %s"%(iobjFrame_o0,1,int(c/2),1,r,sub_osubimg1))
        if os.path.exists(sub_osubimg2): 
            os.remove(sub_osubimg2)
        os.system("imcopy %s[%d:%d,%d:%d] %s"%(iobjFrame_o0,int(c/2)+1,c,1,r,sub_osubimg2))

        subhdr1=fits.getheader(sub_osubimg1)
        subhdr2=fits.getheader(sub_osubimg2)
        subhdr1['GAIN']=subhdr1['GAIN1']
        subhdr1['RNOISE']= subhdr1['RNOISE1']
        subhdr1['SATURATE']= subhdr1['SATURAT1']
        subhdr1['SKYBKG']= subhdr1['SKYBKG1']
        subhdr1['SKYRMS']= subhdr1['SKYRMS1']


        subhdr2['GAIN']=subhdr2['GAIN2']
        subhdr2['RNOISE']= subhdr2['RNOISE2']
        subhdr2['SATURATE']= subhdr2['SATURAT2']
        subhdr2['SKYBKG']= subhdr2['SKYBKG2']
        subhdr2['SKYRMS']= subhdr2['SKYRMS2']

        del subhdr1['GAIN1']
        del subhdr1['GAIN2']
        del subhdr2['GAIN1']
        del subhdr2['GAIN2']

        del subhdr1['RNOISE1']
        del subhdr1['RNOISE2']
        del subhdr2['RNOISE1']
        del subhdr2['RNOISE2']
    

        del subhdr1['SATURAT1']
        del subhdr1['SATURAT2']
        del subhdr2['SATURAT1']
        del subhdr2['SATURAT2']
        
        del subhdr1['SKYBKG1']
        del subhdr1['SKYBKG2']
        del subhdr2['SKYBKG1']
        del subhdr2['SKYBKG2']


        del subhdr1['SKYRMS1']
        del subhdr1['SKYRMS2']
        del subhdr2['SKYRMS1']
        del subhdr2['SKYRMS2']
        
        submat1=fits.getdata(sub_osubimg1)
        submat2=fits.getdata(sub_osubimg2)
        fits.writeto(sub_osubimg1,submat1,subhdr1,overwrite=True) 
        fits.writeto(sub_osubimg2,submat2,subhdr2,overwrite=True) 


        sci_osubimg1=iobjFrameNew[:-15]+'_sciimg_o1.fits'
        sci_osubimg2=iobjFrameNew[:-15]+'_sciimg_o2.fits'
        scimat1=fits.getdata(sci_osubimg1)
        scimat2=fits.getdata(sci_osubimg2)
        fits.writeto(sci_osubimg1,scimat1,subhdr1,overwrite=True)
        fits.writeto(sci_osubimg2,scimat2,subhdr2,overwrite=True)

 
        bkg_osubimg1=iobjFrameNew[:-15]+'_bkg_o1.fits'
        bkg_osubimg2=iobjFrameNew[:-15]+'_bkg_o2.fits'
        bkgmat1=fits.getdata(bkg_osubimg1)
        bkgmat2=fits.getdata(bkg_osubimg2)
        fits.writeto(bkg_osubimg1,bkgmat1,subhdr1,overwrite=True)
        fits.writeto(bkg_osubimg2,bkgmat2,subhdr2,overwrite=True)



        hdrlist.append([dateobs,exptime])
        ldacList.write(ildacn+"\n")
        subImgList.write(iobjFrame_o0+"\n")

        #os.system('rm %s'%(iobjFrameNew))
        #os.system('rm %s'%(iobjFrameNew[-4]+'ldac'))
    
    ldacList.close()
    subImgList.close()
    hdrlistname=scipath+ttfname + "_time_o0.npy"
    save_npy(hdrlistname,hdrlist)
    logger.info("^_^ ALL astrometry and photometry processes are ok!")
 

def pro(tid_target_filter,date):
    rootpath=get_current_dir()
    photometry(rootpath,tid_target_filter,date)

    logger.info('the astrometric process is done')

#tid_target_filter = sys.argv[1]
#date = sys.argv[2] 
###ccdindex= int(sys.argv[3])
#pro(tid_target_filter,date)



