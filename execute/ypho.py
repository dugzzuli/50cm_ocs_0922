import time
from execute.flatdiff import mkdir 
from loguru import logger
from collections import Counter
from scipy.spatial import cKDTree as ckdt
from astropy.stats import sigma_clip
from astropy.io import fits
import numpy as np
import os, sys
from astropy import wcs
from astropy import coordinates as coord, units as u
import datetime
from astropy.time import Time
from astropy.table import Table
from lib.LogInstance import Logger
from lib.phot.ybias import get_current_dir
from execute.plotmag import timeplus
import math

from lib.phot.yrefsed import refSed_gaia
from lib.phot import yrefsed
from execute.yphocom import _great_circle_distance, _spherical_to_cartesian, d2hms, mad #, save_npy
import lib.phot.PlotFun as plotFun
import warnings

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
    ira='0'
    idec='0'
    for i in range(0,len(lpathsubbkg)):
        try:
                hrs=fits.getheader(lpathsubbkg[i])
                
                if(len(hrs["OBJCTRA"])>1):
                    ira  = ":".join(hrs["OBJCTRA"].split())
                    idec = ":".join(hrs["OBJCTDEC"].split())
                    imgs,hrs=fits.getdata(lpathsubbkg[i],header=True)
                break
        except:
                logger.info('these are no parameters of object ra and dec in header of this fits')
                continue
              
    if (ira=='0' or idec=='0'):
        logger.info('these are no parameters of object ra and dec in header of this fits')
        return None
    else:
        ira, idec = d2hms(ira, idec, conv=1)
         
    
     
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
  

    # 3S configuration
    sexParam1  = confdir + "default.param"
    sexParam2  = confdir + "default.paramNew"
    sexConf    = confdir + "default.sex"
    swarpConf  = confdir + "default.swarp"
    scampConf  = confdir + "default.scamp"
    scampConf2  = confdir + "second.scamp"
    psfexConf  = confdir + "default.psfex"
    
    
 
    sexComd1   = "sex %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s "
    sexComd2   = "-DETECT_MINAREA 5 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 "
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
            iimgMat,ihdr= fits.getdata(iobjFrameNew,header=True)
        except:
            logger.info('the fits file is broken',iobjFrameNew)
            continue

        binfact=int(ihdr['XBINNING'])
        pscale  = 0.297 * binfact # pixel scale in unit of arcsec
        exptime=ihdr['EXPTIME']
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
        ihdr["CD1_1"]   = pscale/3600.0
        ihdr["CD1_2"]   = 0.0
        ihdr["CD2_1"]   = 0.0
        ihdr["CD2_2"]   = pscale/3600.0
        fits.writeto(iobjFrame_o0,iimgMat,ihdr,overwrite=True)

        # run sextractor to extract objects
        ildacn = iobjFrame_o0[:-4] + "ldac"
        saturate_value=ihdr['SATURATE']
        isexComd = sexComd%(iobjFrame_o0,sexConf,sexParam1,ildacn,saturate_value)
        os.system(isexComd)
        # construct PSF model
        ipsfComd = psfexComd%(ildacn, psfexConf, 21, 21, tscidir)
        os.system(ipsfComd)
        try:
            os.system("mv *pdf %s"%(figdir))
        
        except:
            logger.info('these are no pdf file in this directory')
            logger.error(str(traceback.format_exc()))
        try:
            
            os.system("mv *fits %s"%tscidir)
        except:
            logger.info('these are no fits file in this directory')
            logger.error(str(traceback.format_exc()))
    
        logger.info("   Run Scamp to refine the astrometry")
     
        
        ildac=fits.getdata(ildacn,ext=2)
            
        if(filtID =='v' or filtID =='u'):
            ixcenloc=np.where((ildac["X_IMAGE"]>c/2+10) | (ildac["X_IMAGE"]<c/2-10))
            #ixyloc=np.where((iximg<c-(1000/binfact))&(iximg>1000/binfact))#&(iyimg<6100)&(iyimg>30))

            ildac =ildac[ixcenloc]
        hdul=fits.open(ildacn)###
        hdul[1].data=ildac
        hdul.close()
            
        #ret=subprocess.getoutput
        iscampComd = scampComd%(ildacn,scampConf,refCatn)
        os.system(iscampComd)
        ihead   = ildacn[:-4] + "head"
        ihead1 = ildacn[:-4] + "ahead"
            
        os.system("mv %s %s"%(ihead,ihead1))
        iscampComd2 = scampComd%(ildacn,scampConf2,refCatn)
        os.system(iscampComd2)

        os.system("mv %s %s"%(ihead,ihead1))
        iscampComd3 = scampComd%(ildacn,scampConf,refCatn)
        os.system(iscampComd3)

        os.system("mv %s %s"%(ihead,ihead1))
        iscampComd4 = scampComd%(ildacn,scampConf2,refCatn)
        os.system(iscampComd4)

 
 
        ihead   = ildacn[:-4] + "head"
        iimgMat, ihdr = fits.getdata(iobjFrame_o0,header=True)
        if(os.path.exists(ihead)==False):
            logger.error(f"os.path.exists(ihead)==False: {os.path.exists(ihead)==False}")
            logger.info(f"scamp can not run a new ihead")

            return None,None
        
        with open(ihead, "r") as iheadList_fits:
        
            iheadList=iheadList_fits.read().splitlines()[3:-1]
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
        try:
            os.system("mv *pdf %s"%(figdir))
        
        except:
            logger.info('these are no pdf file in this directory')
            logger.error(str(traceback.format_exc()))
         
 
        
       # generate the final catalog
        logger.info("############ Generate final SExtractor catalog")
        isexComd = sexComdPSF%(iobjFrame_o0,sexConf,sexParam2,ildacn,saturate_value,ipsfFrameNew)
        os.system(isexComd)
        logger.info("############ Generate ok")


        ## some diagnosis charts and files
        ildac        = fits.getdata(ildacn,ext=2)
        iflag        = ildac["FLAGS"]
        isnr         = ildac["SNR_WIN"]
        iximg        = ildac["X_IMAGE"]
        iyimg        = ildac["Y_IMAGE"]
        ira          = ildac["ALPHA_J2000"]
        idec         = ildac["DELTA_J2000"]
       
        ifwhm        = ildac["FWHM_IMAGE"]
        iaimg        = ildac["AWIN_IMAGE"]
        ibimg        = ildac["BWIN_IMAGE"]
        itheta       = ildac["THETAWIN_IMAGE"]
        iflux_radius= ildac["FLUX_RADIUS"]
        #logger.info(np.array(ildac["MAG_AUTO"])[1]) 
        hdul=fits.open(ildacn)
        hdul[2].data=ildac
        #logger.info(np.array(hdul[1].data[0][0][5]),math.log10(exptime/10)*2.5,exptime)
        ildac["MAG_AUTO"]=ildac["MAG_AUTO"]+math.log10(exptime)*2.5
        ildac["MAG_PSF"] = ildac["MAG_PSF"]+math.log10(exptime)*2.5
        ildac["MAG_APER"]=ildac["MAG_APER"]+math.log10(exptime)*2.5
       
      
        ildac["FLUX_AUTO"] =ildac["FLUX_AUTO"]/float(exptime)
        ildac["FLUX_PSF"] =ildac["FLUX_PSF"]/float(exptime)
        ildac["FLUX_APER"] =ildac["FLUX_APER"]/float(exptime)
        
        ##fy 20230217#####
        len_ldac = len(ildac["FLUX_PSF"])
        for i in range(0,len_ldac):
            try:
                if(ildac["MAG_AUTO"][i]==0):
                    ildac["MAG_AUTO"][i] =  ildac["MAG_AUTO"][i]
                else:
                    ildac["MAG_AUTO"][i] = 1.0857 * (ildac["FLUXERR_AUTO"][i]/ildac["FLUX_AUTO"][i])
            except Exception as e:    
                ildac["MAG_AUTO"][i] =  ildac["MAG_AUTO"][i]
                
            try:
                ildac["MAG_APER"][i] = 1.0857 * (ildac["FLUXERR_APER"][i]/ildac["FLUX_APER"][i])
            except Exception as e:    
                ildac["MAG_APER"][i]=ildac["MAG_APER"][i]
                
            try:
                if(ildac["MAG_APER"][i]==0):
                    ildac["MAG_PSF"][i]=  ildac["MAG_PSF"][i]
                else:
                    ildac["MAG_PSF"][i]= 1.0857 * (ildac["FLUXERR_PSF"][i]/ildac["FLUX_PSF"][i])
            except Exception as e:
                ildac["MAG_PSF"][i]=  ildac["MAG_PSF"][i]
        ##fy 20230217#####
         
        hdul.writeto(ildacn,overwrite=True)
        hdul.close()
        
        ildac        = fits.getdata(ildacn,ext=2)
        imagAUTO     = ildac["MAG_AUTO"] 
        imagPSF      = ildac["MAG_PSF"] 
        imagErrPSF   = ildac["MAGERR_PSF"]
        ifluxAUTO    = ildac["FLUX_AUTO"]
        ifluxErrAUTO = ildac["FLUXERR_AUTO"]
        ifluxPSF     = ildac["FLUX_PSF"]
        ifluxErrPSF  = ildac["FLUXERR_PSF"]
        ifluxAPER     = ildac["FLUX_APER"]
        ifluxErrAPER  = ildac["FLUXERR_APER"]

         
        # 1) region file
        iregn  = iobjFrame_o0[:-4] + "reg"
        wds9reg(ira,idec,flag=None,radius=8.0,unit="arcsec",color="green",outfile=iregn)
        psf_id    = iflux_radius>1.5
        # select high-snr stars
        gid        = (iflag==0) & (isnr>20.0)
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
         
        logger.info("Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))
  
        # idra  = 3600.0*(ira[ids]-raRef[idRef])
        # iddec = 3600.0*(idec[ids]-decRef[idRef])
        # iraSig, idecSig = mad(idra), mad(iddec)
        # logger.info("    Astrometric scatter (%s stars): sigma_ra=%.3f, sigma_dec=%.3f"%(len(idRef),iraSig,idecSig))
        # ifigname = figdir + "astrometry_%s.pdf"%(ildacn.split("/")[-1][:-5])
        # plotFun.astrometryPlot(idra, iddec, sigma=(iraSig,idecSig), figout=ifigname)

        # 3) determine the zeropoint
        izpfig = figdir + "zpoint_%s.pdf"%(ildacn.split("/")[-1][:-5])
        ipsffig = figdir + "psfPhot_%s.pdf"%(ildacn.split("/")[-1][:-5])
        izpArr = magRef[idRef] - imagPSF[ids]
        imask  = sigma_clip(izpArr, sigma=3.0, maxiters=3, masked=True)
        izpArr = izpArr[~imask.mask]
        #izp    = np.median(izpArr)
        izp    = 0
        izpErr = mad(izpArr)
        #plotFun.zpPlot(izpArr, imagPSF[ids][~imask.mask], figout=izpfig)
        logger.info("    Initial zeropoint is %6.3f mag determined by %4d stars"%(izp, len(izpArr)))

        bkgNoise = np.sqrt(np.pi)*np.median(ifwhm) * (np.min([ihdr["SKYRMS1"],ihdr["SKYRMS2"]]))
        #plotFun.mag2Err(ifluxPSF[psf_id],ifluxErrPSF[psf_id],ifluxAUTO[psf_id],ifluxErrAUTO[psf_id],bkgNoise=bkgNoise,zpoint=izp,gain=ihdr["GAIN"],figout=ipsffig)

        # 4) show the Whisker plot
        iwskfig     = figdir + "ewhiskerEll_%s.pdf"%(ildacn.split("/")[-1][:-5])
        iwskfwhmfig = figdir + "ewhiskerFWHM_%s.pdf"%(ildacn.split("/")[-1][:-5])
        #plotFun.ewhisker(iximg, iyimg, iaimg, ibimg, itheta, figout=iwskfig)
        #plotFun.ewhiskerEllipse(iximg, iyimg, iaimg, ibimg, itheta, scale=30.0, figout=iwskfig)
        #plotFun.ewhiskerFWHM(iximg,iyimg,ifwhm,itheta,scale=20,figout=iwskfwhmfig)

     
     



        #hdrlist.append([dateobs,exptime])
        ldacList.write(ildacn+"\n")
        subImgList.write(iobjFrame_o0+"\n")

        #os.system('rm %s'%(iobjFrameNew))
        #os.system('rm %s'%(iobjFrameNew[-4]+'ldac'))
    
    ldacList.close()
    subImgList.close()
    #hdrlistname=scipath+ttfname + "_time_o0.npy"
    #save_npy(hdrlistname,hdrlist)
    logger.info("^_^ ALL astrometry and photometry processes are ok!")
 

def pro(tid_target_filter,date):
    rootpath=get_current_dir()
    photometry(rootpath,tid_target_filter,date)

    logger.info('the astrometric process is done')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20221216', help='输入处理日期')
    parser.add_argument('--tid_target_filter', type=str,default="y50a_ZTFJ2130+4420_mg", help='输入处理日期')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_ypho").get_logger
    pro(args.tid_target_filter,args.date)
    
# tid_target_filter = sys.argv[1]
# date = sys.argv[2] 
# ##ccdindex= int(sys.argv[3])
# pro(tid_target_filter,date)



