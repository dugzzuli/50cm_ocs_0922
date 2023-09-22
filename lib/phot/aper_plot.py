import os
import sys
from astropy.io import fits
from loguru import logger
import matplotlib.pyplot as plt
import numpy as np
from loguru import logger
from lib.LogInstance import Logger
from lib.phot.ybias import get_current_dir


 
 
def getindex_location(tra,tdec,ralist,declist):
    tmpval = 9999
    tmp= 1.5

    for i in range(len(ralist)):
        d=np.sqrt((ralist[i]-tra)**2+(declist[i]-tdec)**2)
        #logger.info(ra[i],',',dec[i],',d=',d)
        if(d<3/3600):
            if ( d < tmpval):
                tmpval = d
                tmp  = i
    return tmp


def single_frame_aper(rootpath,ttfname,date,tra,tdec):
    #tra=322.7362855
    #tdec=44.346236
    
    #rra=322.73015
    #rdec=44.3311416

    #cra=322.6965166
    #cdec=44.3353833

    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'
    #ldacListn_b   = tscidir +'y50b_ZTF2130+4420_jb_o0_ldac.list'
    #ldacListn_v   = tscidir +'y50a_ZTF2130+4420_jv_o0_ldac.list'
    ldac_bname = tscidir +'y50b_ZTF2130+4420_jb/y50b_ZTF2130+4420_jb_001_subbkg_o0.ldac'
    ldac_vname = tscidir +'y50a_ZTF2130+4420_jv/y50a_ZTF2130+4420_jv_001_subbkg_o0.ldac'
    ildac_b = fits.getdata(ldac_bname,ext=2)
    ildac_v = fits.getdata(ldac_vname,ext=2)
    vra, vdec, vid, vpx, vpy, vamag,vapmag,vapmagerr,vflux_aper,vflux_apererr= ildac_v["ALPHA_J2000"], ildac_v["DELTA_J2000"],ildac_v["NUMBER"],ildac_v["X_IMAGE"],ildac_v["Y_IMAGE"],ildac_v["MAG_AUTO"],ildac_v["MAG_APER"],ildac_v["MAGERR_APER"],ildac_v["FLUX_APER"],ildac_v["FLUXERR_APER"]
    bra, bdec, bid, bpx, bpy, bamag,bapmag, bapmagerr,bflux_aper,bflux_apererr= ildac_b["ALPHA_J2000"], ildac_b["DELTA_J2000"],ildac_b["NUMBER"],ildac_b["X_IMAGE"],ildac_b["Y_IMAGE"],ildac_b["MAG_AUTO"],ildac_b["MAG_APER"],ildac_b["MAGERR_APER"],ildac_b["FLUX_APER"],ildac_b["FLUXERR_APER"] 
    logger.info(vapmag[1], type(vapmag[1]))
    indexb=getindex_location(tra,tdec,bra, bdec)
    indexv=getindex_location(tra,tdec,vra, vdec)
    logger.info(indexb,indexv)
    logger.info(bapmag[indexb],bapmagerr[indexb],vapmag[indexv],vapmagerr[indexv])
    return bapmag[indexb],bapmagerr[indexb],vapmag[indexv],vapmagerr[indexv],bflux_aper[indexb]/bflux_apererr[indexb],vflux_aper[indexv]/vflux_apererr[indexv]

def plotshow_combine(x,ty,terr,cy,cerr,ry,rerr,tffname,filename):
    plt.switch_backend('agg')
    #x=[3.0,3.5, 4.0,4.5, 5.0,5.5, 6.0,6.5, 7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,16.0,17.0,18.0,19.0,20.0]
    x=[4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60]
    

    tyn=[]
    cyn=[]
    ryn=[]
    for i in range(1,len(ty)):
        logger.info(ty[i]-ty[i-1])
        tyn.append(ty[i]-ty[i-1])
        cyn.append(cy[i]-cy[i-1])
        ryn.append(ry[i]-ry[i-1])
    tyn=np.array(tyn)
    cyn=np.array(cyn)
    ryn=np.array(ryn)
    #plt.errorbar(x, ty-np.max(ty), yerr=terr,fmt='o',ecolor='b',color='b',elinewidth=2,capsize=4)
    #plt.errorbar(x, cy-np.max(cy), yerr=cerr,fmt='o',ecolor='y',color='y',elinewidth=2,capsize=4)
    #plt.errorbar(x, ry-np.max(ry), yerr=rerr,fmt='o',ecolor='r',color='r',elinewidth=2,capsize=4)
    
    plt.scatter(x, tyn, color='b',s=6)
    plt.scatter(x, cyn, color='y',s=6)
    plt.scatter(x, ryn, color='r',s=6)


    # plt.scatter(x, ty-ry, c='b',s=3)
    # plt.scatter(x,cy-ry,c='c',s=3)
    #plt.xlim(0, 0.7)
    plt.xlabel('MJD-59334.3(days)')
    filterid=tffname[-1]
    #plt.ylabel('B-V__delta_mag(mag)')
    plt.ylabel('aper_mag(mag[+1]-mag)')
       
    #plt.clf()
    ax = plt.gca()
    #ax.invert_yaxis()
    plt.savefig(filename)
   



def multiple_frame_aper(rootpath,ttfname,date):
    tra=322.7362855
    tdec=44.346236

    rra=322.73015
    rdec=44.3311416

    cra=322.6965166
    cdec=44.3353833
    
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'
    #x=[3.0, 3.25,3.5,3.75,4.0,4.25,4.5,4.75, 5.0,5.25,5.5,5.75, 6.0,6.25,6.5,6.75,7.0,7.25,7.5,7.75, 8.0,8.25,8.5,8.75, 9.0, 9.25,9.5,9.75,10.0,10.25 ]
    x=[2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60]
    tbam,tbamerr,tvam,tvamerr,tbsn,tvsn=single_frame_aper(rootpath,ttfname,date,tra,tdec)
    rbam,rbamerr,rvam,rvamerr,rbsn,rvsn=single_frame_aper(rootpath,ttfname,date,rra,rdec)
    cbam,cbamerr,cvam,cvamerr,cbsn,cvsn=single_frame_aper(rootpath,ttfname,date,cra,cdec)
    figname1=tscidir+date+'_bc_aper_mag.pdf'
    plotshow_combine(x,tbam,tbamerr,cbam,cbamerr,rbam,rbamerr,ttfname,figname1)
    figname2=tscidir+date+'_vc_aper_mag.pdf'
    plotshow_combine(x,tvam,tvamerr,cvam,cvamerr,rvam,rvamerr,ttfname,figname2)
    plt.switch_backend('agg')
    plt.scatter(x,tbsn, c='b',s=4)
    plt.scatter(x,rbsn, c='r',s=4)
    plt.scatter(x,cbsn, c='y',s=4)
    # plt.scatter(x,cy-ry,c='c',s=3)
    #plt.xlim(0, 0.7)
    plt.xlabel('MJD-59334.3(days)')
    #plt.ylabel('B-V__delta_mag(mag)')
    plt.ylabel('b_(aper_flux/aper_flux_err)')
    #plt.clf()
    ax = plt.gca()
    plt.savefig(tscidir+date+'_b_aper_sn.pdf') 
    plt.close()
    plt.switch_backend('agg')
    plt.scatter(x,tbsn, c='b',s=4)
    plt.scatter(x,rbsn, c='r',s=4)
    plt.scatter(x,cbsn, c='y',s=4)
    # plt.scatter(x,cy-ry,c='c',s=3)
    #plt.xlim(0, 0.7)
    plt.xlabel('MJD-59334.3(days)')
    #plt.ylabel('B-V__delta_mag(mag)')
    plt.ylabel('v_(aper_flux/aper_flux_err)')
    #plt.clf()
    ax = plt.gca()
    plt.savefig(tscidir+date+'_v_aper_sn.pdf')


def pro(ttfname,date):

    rootpath=get_current_dir()
    multiple_frame_aper(rootpath,ttfname,date)

    logger.info('the astrometric process is done')


               

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20221207', help='输入处理日期')
    parser.add_argument('--ttfname', type=str, help='输入处理日期')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_yphost").get_logger
    pro(args.ttfname,args.date)
    
# ttfname = sys.argv[1]
# date = sys.argv[2]
# pro(ttfname,date)

