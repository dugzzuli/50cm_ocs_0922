# obtain the gaia catalog for a given sky region
import os, sys
from collections import Counter
from scipy.spatial import cKDTree as ckdt
from astropy.stats import sigma_clip
from astropy.io import fits
import numpy as np
import os, sys
from scipy.optimize import leastsq
import time
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import math
from loguru import logger
import time
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from execute import yphotutils
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

def filtTransform(magG, magBP, magRP):
    """
    filter transformation between Gaia and Y50 BVRI
    
    Reference:
    https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/
    chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html

    Parameters:
    magG, 1-d array
        Gaia G-band magnitudes
    magBP, 1-d array
        Gaia BP-band magnitudes
    magRP, 1-d array
        Gaia RP-band magnitudes
    """
    dBR = magBP - magRP
    if dBR<-0.5 or dBR>2.75: return np.nan, np.nan, np.nan, np.nan

    magV = magG - (-0.017600 - 0.00686*dBR - 0.17320*dBR**2)
    magR = magG - (-0.003226 + 0.38330*dBR - 0.13450*dBR**2)
    magI = magG - (+0.020850 + 0.74190*dBR - 0.09631*dBR**2)
    
    dGV = magG-magV
    if dGV<-1.4 or dGV>-0.03: return np.nan, magV, magR, magI
    
    # B-band
    a, b, c, d = -0.001768, -0.229700, -0.023850, -0.029070-dGV
    p = (b*c)/(6.0*a*a) - d/(2*a) - (b*b*b)/(27.0*a*a*a)
    q = p*p + (c/(3.0*a)-(b*b)/(9.0*a*a))**3
    p, q = complex(p), complex(q)
    ppq = p + np.sqrt(q)
    pmq = p - np.sqrt(q)
    u = ppq**(1/3)
    v = pmq**(1/3)
    w = complex(-0.5, 0.5*np.sqrt(3))
    dBV1 = u + v - b/(3.0*a)
    dBV2 = w*u + w**2*v - b/(3.0*a)
    dBV3 = w**2*u + w*v - b/(3.0*a)
    dBVR = [dBV1, dBV2, dBV3]    

    dBV  = [i.real for i in dBVR if i.real>0][0]
    magB = dBV + magV

    return magB, magV, magR, magI

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

def crossmatch(ra1, dec1, ra2, dec2, aperture=1.5):
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


def RefCatDownload(refsavepath,targetname,fra,fdec):
    # refsavepath = '/Users/cxl/50cm/refcata/'
    refpath = '%sGaiaDR3_%s.txt'%(refsavepath,targetname)
    if os.path.exists(refpath) != True:
        print(fra,fdec)
        ra,dec=d2hms(fra,fdec,conv=0)
        print('@@@@@@@@@@@@@@')
        print(ra,dec)
        rah,ram,ras = ra.split(':')
        decd,decm,decs = dec.split(':')
        
        print(rah,ram,ras)
        print(refpath)
        url = 'curl -o %s \
            "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=gaia_dr3_source&spatial=box&radunits=arcsec&objstr=%sh+%sm+%ss+%sd+%sm+%ss&size=4000&outfmt=1&selcols=ra,dec,ref_epoch,phot_g_mean_mag"'%(refpath,rah,ram,ras,decd,decm,decs)
        os.system(url)
    return refpath


def ScampRefCat(refpath,blimmag=0,dlimmag=25):
    scamprefpath = refpath.split('.')[0]+'.ldac'
    scamprefregpath = refpath.split('.')[0]+'.reg'
    #if os.path.exists(scamprefpath) != True or os.path.exists(scamprefregpath) != True:
    gaiacata0 = Table.read(refpath,format='ipac')
    gaiacata = gaiacata0[(blimmag<gaiacata0['phot_g_mean_mag']) & (gaiacata0['phot_g_mean_mag']<dlimmag)]
    gaiara = np.array(gaiacata['ra'])
    gaiadec = np.array(gaiacata['dec'])
    gaiamag = np.array(gaiacata['phot_g_mean_mag'])
    gaiamagerr = np.ones((len(gaiara))) * 0.010857
    obsdate = np.array(gaiacata['ref_epoch'])
    flags = np.ones((len(gaiara))) * 0
    aerr = np.ones((len(gaiara))) * 0
    berr = np.ones((len(gaiara))) * 0
    thetaerr = np.ones((len(gaiara))) * 0
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    col01 = fits.Column(name="OBJECT_POS",format="J",array=[1])
    col02 = fits.Column(name="OBJECT_COUNT",format="J",array=[len(gaiamag)])
    col03 = fits.Column(name="CHANNEL_NR",format="J",array=[0])
    col04 = fits.Column(name="CHANNEL_NAME",format="16A",array=["I"])
    cols = fits.ColDefs([col01,col02,col03,col04])
    tbhdu1 = fits.BinTableHDU.from_columns(cols)
    tbhdu1.header["EXTNAME"] = "LDAC_IMHEAD"
    tbhdu1.header["AUTHOR"] = "LDAC by DZ"
    col1 = fits.Column(name='X_WORLD', format='1D',unit="deg", array=gaiara)
    col2 = fits.Column(name='Y_WORLD', format='1D', unit="deg", array=gaiadec)
    col3 = fits.Column(name='ERRA_WORLD', format='1E', unit="deg", array=aerr)
    col4 = fits.Column(name='ERRB_WORLD', format='1E', unit="deg", array=berr)
    col5 = fits.Column(name='ERRTHETA_WORLD', format='1E', unit="deg", array=thetaerr)
    col6 = fits.Column(name='MAG', format='E', array=gaiamag)
    col7 = fits.Column(name='MAGERR', format='1E', array=gaiamagerr)
    col8 = fits.Column(name='FLAGS', format='1I', array=flags)
    col9 = fits.Column(name='OBSDATE', format='1D', array=obsdate)
    #cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9])
    cols = fits.ColDefs([col1,col2,col3,col4,col6,col7,col8,col9])
    tbhdu2 = fits.BinTableHDU.from_columns(cols)
    tbhdu2.header["EXTNAME"] = "LDAC_OBJECTS"
    thdulist = fits.HDUList([prihdu,tbhdu1,tbhdu2])
    thdulist.writeto(scamprefpath,overwrite=True)
    reglist = open(scamprefregpath,'w')
    reglist.write('# Region file for DS9'+"\n")
    reglist.write('global color=red font=\'helvetica 10 normal\' select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'+"\n")
    reglist.write('fk5'+"\n")
    for i in np.arange(len(gaiara)):
        regionline='circle(%s, %s, 3.00")'%(gaiara[i],gaiadec[i])
        reglist.write(regionline+"\n")
    return scamprefpath


def ldac_scamp(x,y,mag,a_err=None,b_err=None,theta_err=None,mag_err=None,
               obsdate=None,flag=None,outcat="scamp.ldac"):
    """
    create LDAC for Scamp use
    """
    if type(a_err)==type(None): a_err = mag*0.0
    if type(b_err)==type(None): b_err = mag*0.0
    if type(theta_err)==type(None): theta_err = mag*0.0
    if type(mag_err)==type(None): mag_err = mag*0.0+1.0857/100.0
    if type(obsdate)==type(None): obsdate = mag*0.0 + 2015
    if type(flag)==type(None): flag=mag*0.0

    # first construct an empty fits table
    prihdr = fits.Header()
    prihdr["AUTHOR"]="DZ LIU"
    prihdu = fits.PrimaryHDU(header=prihdr)

    # First extension
    col01 = fits.Column(name="OBJECT_POS",format="J",array=[1])
    col02 = fits.Column(name="OBJECT_COUNT",format="J",array=[len(mag)])
    col03 = fits.Column(name="CHANNEL_NR",format="J",array=[0])
    col04 = fits.Column(name="CHANNEL_NAME",format="16A",array=["I"])
    cols = fits.ColDefs([col01,col02,col03,col04])
    tbhdu1 = fits.BinTableHDU.from_columns(cols)
    tbhdu1.header["EXTNAME"] = "LDAC_IMHEAD"
    tbhdu1.header["AUTHOR"] = "LDAC by DZ"

    # Second extension
    col1 = fits.Column(name='X_WORLD', format='1D',unit="deg", array=x)
    col2 = fits.Column(name='Y_WORLD', format='1D', unit="deg", array=y)
    col3 = fits.Column(name='ERRA_WORLD', format='1E', unit="deg", array=a_err)
    col4 = fits.Column(name='ERRB_WORLD', format='1E', unit="deg", array=b_err)
    col5 = fits.Column(name='ERRTHETA_WORLD', format='1E', unit="deg", array=theta_err)
    col6 = fits.Column(name='MAG', format='E', array=mag)
    col7 = fits.Column(name='MAGERR', format='1E', array=mag_err)
    col8 = fits.Column(name='FLAGS', format='1I', array=flag)
    col9 = fits.Column(name='OBSDATE', format='1D', array=obsdate)

    #cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9])
    cols = fits.ColDefs([col1,col2,col3,col4,col6,col7,col8,col9])
    tbhdu2 = fits.BinTableHDU.from_columns(cols)
    tbhdu2.header["EXTNAME"] = "LDAC_OBJECTS"

    thdulist = fits.HDUList([prihdu,tbhdu1,tbhdu2])
    thdulist.writeto(outcat,overwrite=True)
    return

def refSed_gaia(rootpath,ttfname,date,refSky_ra,refSky_dec,maglim_min=10.0,maglim_max=21.0):
    upath=rootpath+'reception/'+str(date)+'/'
    scipath = upath+'sci/'+ttfname+'/'
    
    #basic setting
    #dirname, filename = os.path.split(os.path.abspath(__file__))
    imgdir  = rootpath + "images/"
    confdir = rootpath + "config/"
    ancdir  = rootpath + "ancData/"
    figdir  = rootpath + "figures/"
  
    t0 = time.time()  

    tid,objID,filterid=ttfname.split('_')
    from config import config

    refCatDir =config["refcat"]["gaia"] #'/home/50cm/data2/GAIAEDR3New/'   
    # directory of the reference gaia catalogs

    # target coordinate in degree
    raCen, decCen = float(refSky_ra),float(refSky_dec)
    magGLim       = [maglim_min,maglim_max]# magnitude limit of the Gaia G-band 
    catSkyn =  "GaiaStar_%s_%s.cat"%(objID,date)
    #catSkyn       = "GaiaStar_%s_%s.cat"%(objID,date)
    width     = 1
    ##############################################
    #coslim=math.cos(decCen/180*math.pi)
    #logger.info('coslim=',coslim)
 
    # ra1, ra2   = raCen-0.5*width,  raCen+0.5*width
    # #dec1, dec2 = decCen-0.5*width, decCen+0.5*width
    # dec1, dec2 = decCen-0.5*width/coslim, decCen+0.5*width/coslim
    #ra1, ra2   = raCen-0.5*width,  raCen+0.5*width
    # #dec1, dec2 = decCen-0.5*width, decCen+0.5*width
    #dec1, dec2 = decCen-(0.5*width*coslim), decCen+ (0.5*width*coslim)
    # if(decCen<0.5*width):
    #  dec1, dec2 = 0, decCen+0.5*width/coslim

    #ra1,ra2,dec1,dec2=14.600,   15.600, 43.447,   44.841
    logger.info(f'decCen= {decCen}')
    logger.info(type(decCen))

    coslim = math.cos(decCen*(math.pi/180))
    logger.info(f'coslim={coslim}')
    ra1, ra2   = raCen-(0.5*width/coslim),  raCen+(0.5*width/coslim)
    dec1, dec2 = decCen-(0.5*width), decCen+ (0.5*width)
 
    tarSky     = (ra1, dec1, ra2, dec2)

    logger.info("^_^ Extract GAIA catalog in target sky: %s"%objID)
    logger.info("^_^ Sky coverage: ra=[%8.3f, %8.3f], cosdec=[%8.3f, %8.3f], dec=[%8.3f, %8.3f]"%(ra1,ra2,dec1,dec2,decCen-0.5*width,decCen+0.5*width))
    logger.info("^_^ Magnitude limit: Gmag=[%7.3f, %7.3f]"%(magGLim[0],magGLim[1]))
    logger.info(" ")

    # catalog list
    rcatDir=refCatDir
    catListn = rcatDir + "catalog.info"
    catList  = open(catListn,"r").read().splitlines()[1:]
    ncat     = len(catList)
    logger.info("^_^ Total %s GAIA catalogs"%ncat)

    catSky   = open(ancdir+catSkyn,"w")
    title    = "#id ra dec epoch mag magBP magRP magB magV magR magI\n"
    fmt      = "%10d %15.8f %15.8f %8.1f %8.3f %8.3f %8.3f "
    catSky.write(title)
    nstar    = 0
    badstar  = 0
    obsTime= f"{date[:4]}-{date[4:6]}-{date[6:]}"
    for i in range(ncat):
        icatInf = catList[i].split()
        icatn   = icatInf[0]
        iraMin, iraMax   = float(icatInf[1]), float(icatInf[2])
        idecMin, idecMax = float(icatInf[3]), float(icatInf[4])
        iSky  = (iraMin, idecMin, iraMax, idecMax)
        iarea = overlapRect(tarSky, iSky)
        if iarea is None: continue
        logger.info("    Catalog (%d/%d) %s covers the target sky"%(i+1,ncat,icatn))
        # read the catalog
        icat  = Table.read(rcatDir+icatn, format="fits")
        inobj = len(icat)
        logger.info("    Total %s stars in the catalog"%inobj)
        ira, idec, imag = icat["ra"], icat["dec"], icat["phot_g_mean_mag"]
        imagBP, imagRP  = icat["phot_bp_mean_mag"], icat["phot_rp_mean_mag"]
        iepoch          = icat["ref_epoch"]
         
        ipmra           = icat["pmra"]
        ipmdec          = icat["pmdec"]
        idup            = icat["duplicated_source"]
        for j in range(inobj):
            if (math.isnan(imag[j]) or type(imag[j]) != np.float64):
                logger.info('the nan value is: ',imag[j])
                continue
            if type(imagBP[j]) != np.float64 or math.isnan(imagBP[j]) : continue
            if type(imagRP[j]) != np.float64 or math.isnan(imagRP[j]): continue
            if imagBP[j]>30.0 or imagRP[j]>30.0: continue
            if imag[j]<magGLim[0] or imag[j]>magGLim[1]: continue
            jpos = (ira[j], idec[j])
            gid  = pointRect(jpos, tarSky)
            if gid:
                nstar += 1
                icrd = SkyCoord(ra=jpos[0]*u.deg, dec=jpos[1]*u.deg, frame="icrs",
                            pm_ra_cosdec=ipmra[j]*u.mas/u.yr, pm_dec=ipmdec[j]*u.mas/u.yr,
                            obstime=Time(iepoch[j], format="jyear"))
                icrdx = icrd.apply_space_motion(new_obstime=Time(obsTime))
                jrax  = icrdx.ra.degree
                jdecx = icrdx.dec.degree
                jrline = fmt%(nstar,jrax,jdecx,iepoch[j],imag[j],imagBP[j],imagRP[j])
                #jrline = fmt%(nstar,ira[j],idec[j],iepoch[j],imag[j],imagBP[j],imagRP[j])
                if "nan" in jrline: badstar+=1; continue
                # add BVRI magnitudes
                #imagB, imagV, imagR, imagI = filtTransform(imag[j], imagBP[j], imagRP[j])
                #jaline = "%8.3f %8.3f %8.3f %8.3f\n"%(imagB, imagV, imagR, imagI)
                jaline = "%8.3f %8.3f %8.3f %8.3f\n"%(0, 0, 0, 0)
                jline  = jrline + jaline
                catSky.write(jline)

    logger.info("^_^ Total %d stars returned"%(nstar-badstar))
    catSky.close()
        
    # covert the ascii catalog to ldac_fits and also give ds9 region file
    ldacName = ancdir + catSkyn[:-3] + "ldac"
    regName  = ancdir + catSkyn[:-3] + "reg"

    pcat       = Table.read(ancdir+catSkyn,format="ascii")
    ra, dec    = pcat["ra"],    pcat["dec"]
    epoch, mag = pcat["epoch"], pcat["mag"]

    ldac_scamp(ra,dec,mag,a_err=None,b_err=None,theta_err=None,mag_err=None,obsdate=epoch,flag=None,outcat=ldacName)
    wds9reg(ra,dec,flag=None,radius=3.0,unit="arcsec",color="green",outfile=regName)

    t1 = time.time()
    dt = t1 - t0
    logger.info("^_^ Total %.2f min used"%(dt/60.0))

def refSed_gaia_filter(rootpath,ttfname,date,refSky_ra,refSky_dec,maglim_min=8.0,maglim_max=21.0):
    ancdir  = rootpath + "ancData/"
    t0 = time.time()  
    temp_arr=ttfname.split('_')
    if len(temp_arr)==4:
        tid,objID,filterid=temp_arr[0],temp_arr[1]+"_"+temp_arr[2],temp_arr[3] #yphotutils.ttfname_header(ttfname)
    else:
        tid,objID,filterid=temp_arr
    # try:
    #     tid,objID,filterid=ttfname.split('_')
    # except:
    #     tid,obj1,obj2,filterid=ttfname.split('_')
    #     objID = obj1+'_'+obj2
    from config import config
    refCatDir = config["refcat"]["gaia"]
    logger.info(f'the date is {date}')
    obsTime= f"{date[:4]}-{date[4:6]}-{date[6:]}"
     
    logger.info(f'the date is {obsTime}')
    # directory of the reference gaia catalogs
    # target coordinate in degree
    raCen, decCen = float(refSky_ra),float(refSky_dec)
     
    magGLim  = [maglim_min,maglim_max]# magnitude limit of the Gaia G-band 
    catSkyn       = "GaiaStar_%s_%s_%s.cat"%(objID,filterid,date)
    width     = 1.5
    ##############################################
    #coslim=math.cos(decCen/180*math.pi)
    logger.info(type(decCen))
    logger.info(f'decCen= {decCen}')
 
    coslim = math.cos(decCen*(math.pi/180))
    logger.info(f'coslim= {coslim}')
    ra1, ra2   = raCen-(0.5*width/coslim),  raCen+(0.5*width/coslim)
    dec1, dec2 = decCen-(0.5*width), decCen+ (0.5*width)
    tarSky     = (ra1, dec1, ra2, dec2)
    logger.info("^_^ Extract GAIA catalog in target sky: %s"%objID)
     
    # catalog list
    rcatDir=refCatDir
    catListn = rcatDir + "catalog.info"
    catList  = open(catListn,"r").read().splitlines()[1:]
    ncat     = len(catList)
    logger.info("^_^ Total %s GAIA catalogs"%ncat)
    catSky   = open(ancdir+catSkyn,"w")
    title    = "#id ra dec epoch mag magBP magRP magB magV magR magI\n"
    fmt      = "%10d %15.8f %15.8f %8.1f %8.3f %8.3f %8.3f "
    catSky.write(title)
    nstar    = 0
    badstar  = 0
    for i in range(ncat):
        icatInf = catList[i].split()
        icatn   = icatInf[0]
        iraMin, iraMax   = float(icatInf[1]), float(icatInf[2])
        idecMin, idecMax = float(icatInf[3]), float(icatInf[4])
        iSky  = (iraMin, idecMin, iraMax, idecMax)
        iarea = overlapRect(tarSky, iSky)
        if iarea is None: continue
        logger.info("    Catalog (%d/%d) %s covers the target sky"%(i+1,ncat,icatn))
        # read the catalog
        icat  = Table.read(rcatDir+icatn, format="fits")
        inobj = len(icat)
        logger.info("    Total %s stars in the catalog"%inobj)
        ira, idec, imag = icat["ra"], icat["dec"], icat["phot_g_mean_mag"]
        imagBP, imagRP  = icat["phot_bp_mean_mag"], icat["phot_rp_mean_mag"]
        iepoch          = icat["ref_epoch"]
         
        ipmra           = icat["pmra"]
        ipmdec          = icat["pmdec"]
        idup            = icat["duplicated_source"]
        for j in range(inobj):
            if (math.isnan(imag[j]) or type(imag[j]) != np.float64):
                logger.info('the nan value is: ',imag[j])
                continue
            if type(imagBP[j]) != np.float64 or math.isnan(imagBP[j]) : continue
            if type(imagRP[j]) != np.float64 or math.isnan(imagRP[j]): continue
            if imagBP[j]>30.0 or imagRP[j]>30.0: continue
            if imag[j]<magGLim[0] or imag[j]>magGLim[1]: continue
            jpos = (ira[j], idec[j])
            gid  = pointRect(jpos, tarSky)
            if gid:
                nstar += 1
                icrd = SkyCoord(ra=jpos[0]*u.deg, dec=jpos[1]*u.deg, frame="icrs",
                            pm_ra_cosdec=ipmra[j]*u.mas/u.yr, pm_dec=ipmdec[j]*u.mas/u.yr,
                            obstime=Time(iepoch[j], format="jyear"))
                icrdx = icrd.apply_space_motion(new_obstime=Time(obsTime))
                jrax  = icrdx.ra.degree
                jdecx = icrdx.dec.degree
                jrline = fmt%(nstar,jrax,jdecx,iepoch[j],imag[j],imagBP[j],imagRP[j])
                #jrline = fmt%(nstar,ira[j],idec[j],iepoch[j],imag[j],imagBP[j],imagRP[j])
                if "nan" in jrline: badstar+=1; continue
                # add BVRI magnitudes
                #imagB, imagV, imagR, imagI = filtTransform(imag[j], imagBP[j], imagRP[j])
                #jaline = "%8.3f %8.3f %8.3f %8.3f\n"%(imagB, imagV, imagR, imagI)
                jaline = "%8.3f %8.3f %8.3f %8.3f\n"%(0, 0, 0, 0)
                jline  = jrline + jaline
                catSky.write(jline)
    logger.info("^_^ Total %d stars returned"%(nstar-badstar))
    logger.info(ancdir+catSkyn)
    catSky.close()
    # covert the ascii catalog to ldac_fits and also give ds9 region file
    ldacName = ancdir + catSkyn[:-3] + "ldac"
    regName  = ancdir + catSkyn[:-3] + "reg"
    pcat       = Table.read(ancdir+catSkyn,format="ascii")
    ra, dec    = pcat["ra"],    pcat["dec"]
    epoch, mag = pcat["epoch"], pcat["mag"]
    ldac_scamp(ra,dec,mag,a_err=None,b_err=None,theta_err=None,mag_err=None,obsdate=epoch,flag=None,outcat=ldacName)
    wds9reg(ra,dec,flag=None,radius=3.0,unit="arcsec",color="green",outfile=regName)
    t1 = time.time()
    dt = t1 - t0
    logger.info("^_^ Total %.2f min used"%(dt/60.0))

  
def refSed_gaia_filter_pm(rootpath,ttfname,date,refSky_ra,refSky_dec,maglim_min=8.0,maglim_max=21.0):
    ancdir  = rootpath + "ancData/"
    t0 = time.time()  
    try:
        tid,objID,filterid=ttfname.split('_')
    except:
        tid,obj1,obj2,filterid=ttfname.split('_')
        objID = obj1+'_'+obj2
    from config import config
    refCatDir = config["refcat"]["gaia"]
    logger.info(f'the date is {date}')
    obsTime= f"{date[:4]}-{date[4:6]}-{date[6:]}"
     
    logger.info(f'the date is {obsTime}')
    # directory of the reference gaia catalogs
    # target coordinate in degree
    raCen, decCen = refSky_ra,refSky_dec
     
    magGLim  = [maglim_min,maglim_max]# magnitude limit of the Gaia G-band 
        #logger.info("^_^ Magnitude limit: Gmag=[%7.3f, %7.3f]"%(magGLim[0],magGLim[1]))
    magGLim  = [maglim_min,maglim_max]# magnitude limit of the Gaia G-band 
    catSkyn       = "GaiaStar_%s_%s_%s.cat"%(objID,filterid,date)
    width     = 1.5
    ##############################################
    #coslim=math.cos(decCen/180*math.pi)
    coslim = math.cos(decCen*(math.pi/180))
    logger.info(f'coslim= {coslim}')
    ra1, ra2   = raCen-(0.5*width/coslim),  raCen+(0.5*width/coslim)
    dec1, dec2 = decCen-(0.5*width), decCen+ (0.5*width)
    tarSky     = (ra1, dec1, ra2, dec2)
    logger.info("^_^ Extract GAIA catalog in target sky: %s"%objID)
     
    # catalog list
    rcatDir=refCatDir
    catListn = rcatDir + "catalog.info"
    catList  = open(catListn,"r").read().splitlines()[1:]
    ncat     = len(catList)
    logger.info("^_^ Total %s GAIA catalogs"%ncat)
    catSky   = open(ancdir+catSkyn,"w")
    title    = "#id ra dec epoch mag magBP magRP magB magV magR magI\n"
    fmt      = "%10d %15.8f %15.8f %8.1f %8.3f %8.3f %8.3f "
    catSky.write(title)
    nstar    = 0
    badstar  = 0
    for i in range(ncat):
        icatInf = catList[i].split()
        icatn   = icatInf[0]
        iraMin, iraMax   = float(icatInf[1]), float(icatInf[2])
        idecMin, idecMax = float(icatInf[3]), float(icatInf[4])
        iSky  = (iraMin, idecMin, iraMax, idecMax)
        iarea = overlapRect(tarSky, iSky)
        if iarea is None: continue
        logger.info("    Catalog (%d/%d) %s covers the target sky"%(i+1,ncat,icatn))
        # read the catalog
        #########################################################pm change################
        icat  = Table.read(rcatDir+icatn, format="fits")
        inobj = len(icat)
        logger.info("    Total %s stars in the catalog"%inobj)
        ira, idec, imag = icat["ra_J2023"], icat["dec_J2023"], icat["phot_g_mean_mag"]
        imagBP, imagRP  = icat["phot_bp_mean_mag"], icat["phot_rp_mean_mag"]
        iepoch          = icat["ref_epoch"]
        ipmra           = icat["pmra"]
        ipmra_err       = icat["pmra_error"]
        ipmdec          = icat["pmdec"]
        ipmdec_err       = icat["pmdec_error"]
        idup            = icat["duplicated_source"]
        istar =icat['classprob_dsc_combmod_star']
        irpmra =ipmra_err/ipmra
        irpmdec =ipmdec_err/ipmdec
        icat = icat[np.where((ipmra!=999.0)&(imagBP>30.0)&(imagRP>30.0) & (imag<magGLim[0]) &(imag>magGLim[1])&(irpmra<1)&(irpmdec<1))&(istar>0.5)]
        for j in range(inobj):
            # if (math.isnan(imag[j]) or type(imag[j]) != np.float64):
            #     logger.info('the nan value is: ',imag[j])
            #     continue
            # if type(imagBP[j]) != np.float64 or math.isnan(imagBP[j]) : continue
            # if type(imagRP[j]) != np.float64 or math.isnan(imagRP[j]): continue
            # if imagBP[j]>30.0 or imagRP[j]>30.0: continue
            # if imag[j]<magGLim[0] or imag[j]>magGLim[1]: continue
            
            
            
            jpos = (ira[j], idec[j])
            gid  = pointRect(jpos, tarSky)
            
            if gid:
                nstar += 1
                jrline = fmt%(nstar,ira[j],idec[j],iepoch[j],imag[j],imagBP[j],imagRP[j])
                if "nan" in jrline: badstar+=1; continue
            
                # add BVRI magnitudes
                #imagB, imagV, imagR, imagI = filtTransform(imag[j], imagBP[j], imagRP[j])
                #jaline = "%8.3f %8.3f %8.3f %8.3f\n"%(imagB, imagV, imagR, imagI)
                jaline = "%8.3f %8.3f %8.3f %8.3f\n"%(0, 0, 0, 0)
                jline  = jrline + jaline
                catSky.write(jline)
                
            # if gid:
            #     nstar += 1
            #     icrd = SkyCoord(ra=jpos[0]*u.deg, dec=jpos[1]*u.deg, frame="icrs",
            #                 pm_ra_cosdec=ipmra[j]*u.mas/u.yr, pm_dec=ipmdec[j]*u.mas/u.yr,
            #                 obstime=Time(iepoch[j], format="jyear"))
            #     icrdx = icrd.apply_space_motion(new_obstime=Time(obsTime))
            #     jrax  = icrdx.ra.degree
            #     jdecx = icrdx.dec.degree
            #     jrline = fmt%(nstar,jrax,jdecx,iepoch[j],imag[j],imagBP[j],imagRP[j])
            #     #jrline = fmt%(nstar,ira[j],idec[j],iepoch[j],imag[j],imagBP[j],imagRP[j])
            #     if "nan" in jrline: badstar+=1; continue
                 
            #     # add BVRI magnitudes
            #     #imagB, imagV, imagR, imagI = filtTransform(imag[j], imagBP[j], imagRP[j])
            #     #jaline = "%8.3f %8.3f %8.3f %8.3f\n"%(imagB, imagV, imagR, imagI)
            #     jaline = "%8.3f %8.3f %8.3f %8.3f\n"%(0, 0, 0, 0)
            #     jline  = jrline + jaline
            #     catSky.write(jline)
    logger.info("^_^ Total %d stars returned"%(nstar-badstar))
    catSky.close()
    # covert the ascii catalog to ldac_fits and also give ds9 region file
    ldacName = ancdir + catSkyn[:-3] + "ldac"
    regName  = ancdir + catSkyn[:-3] + "reg"
    pcat       = Table.read(ancdir+catSkyn,format="ascii")
    ra, dec    = pcat["ra"],    pcat["dec"]
    epoch, mag = pcat["epoch"], pcat["mag"]
    ldac_scamp(ra,dec,mag,a_err=None,b_err=None,theta_err=None,mag_err=None,obsdate=epoch,flag=None,outcat=ldacName)
    wds9reg(ra,dec,flag=None,radius=3.0,unit="arcsec",color="green",outfile=regName)
    t1 = time.time()
    dt = t1 - t0
    logger.info("^_^ Total %.2f min used"%(dt/60.0))
