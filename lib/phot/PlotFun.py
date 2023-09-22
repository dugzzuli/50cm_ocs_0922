# functions to plot
from astropy.stats import sigma_clip
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Ellipse
import pylab as pl
import os, sys
from loguru import logger

 
import math 
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

  
def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)
  


def astrometryPlot(dra, ddec, sigma=(0.0,0.0), figout=None):
    """
    Compare the measured celestial coordinates with reference 
    """
    raSig, decSig = sigma
    sigMax = max([raSig, decSig])
    x0, x1 = -3.0*sigMax, +3.0*sigMax
    y0, y1 = -3.0*sigMax, +3.0*sigMax
    fig = pl.figure(figsize=(6.5,6.0))
    ax  = fig.add_axes([0.15,0.1,0.80,0.85])
    ax.scatter(dra,ddec,marker=".",s=1, c="k")
    ax.plot([x0,x1],[0.0,0.0],"--",color="k",linewidth=2)
    ax.plot([0.0,0.0],[y0,y1],"--",color="k",linewidth=2)
    ax.set_xlim(x0,x1)
    ax.set_ylim(y0,y1)
    ax.set_xlabel("$\Delta$Ra (arcsec) [$\sigma$=%.3f]"%raSig,fontsize=15)
    ax.set_ylabel("$\Delta$Dec (arcsec) [$\sigma$=%.3f]"%decSig,fontsize=15)
    #ax.text(x1-0.29,y1-0.05,"$\sigma_{RA}$=%.3f arcsec"%raSig,fontsize=12)
    #ax.text(x1-0.30,y1-0.10,"$\sigma_{Dec}$=%.3f arcsec"%decSig,fontsize=12)
    for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(15)
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(15)
    pl.savefig(figout)
    pl.clf()
    pl.close()

    return


def astrometryPlot_new(ra, dec, refra,refdec, dra, ddec, sigma=(0.0,0.0), mu=(0.0,0.0),figout=None,mini_figout=None):
    
    if len(dra) ==0:
        logger.error("astrometryPlot fail")
        logger.info("astrometryPlot fail")
        return 
    if  len(ddec)==0:
        logger.error("astrometryPlot fail")
        logger.info("astrometryPlot fail")
        return 
    """
    Compare the measured celestial coordinates with reference 
    """
    # fig,axes = pl.subplots(2,1,figsize=(14,8))  
    # ax = axes.flatten()  
    pl.close()
    pl.style.use('seaborn-paper')
    
    
    # 子图中再绘制子图
    fig = pl.figure(constrained_layout=True,figsize=(12,6))
    
    gs0 = GridSpec(1, 1, figure=fig)#将figure切片为1行2列的两个子图
    
    gs00 = gridspec.GridSpecFromSubplotSpec(6, 2, subplot_spec=gs0[0])#将以上第一个子图gs0[0]再次切片为3行3列的9个axes
    #gs0[0]子图自由切片
    format_axes(fig)
    p1 = fig.add_subplot(gs00[:, :1])
    p2 = fig.add_subplot(gs00[:-3, 1])
    p3 = fig.add_subplot(gs00[-3:, -1])
    raSig, decSig = sigma
    raMu,decMu=mu
    sigMax = max([raSig, decSig])
    x0, x1 = -3.0*sigMax, +3.0*sigMax
    y0, y1 = -3.0*sigMax, +3.0*sigMax 
    logger.info(f'plotfun len(ra)={len(ra)},len(raRef)={len(refra)},len(idra)={len(dra)}')
    p1.scatter( refra,refdec,color='r',alpha = 0.6,linewidths=0.3,s=2)
    for j in range(0,len(refra)):
        
        p1.arrow(refra[j], refdec[j],dra[j]/30,ddec[j]/30,width = 0.001,ec='g',head_width=0.01)
        
    p1.scatter(np.median(refra),np.max(refdec)+0.05,s=4,color='r')
    p1.arrow(np.median(refra),np.max(refdec)+0.05,1/30,0.,width=0.001,ec='g',head_width=0.01)
    p1.text(np.median(refra)+9/180,np.max(refdec)+0.05,'1 arcsec')
    p1.invert_xaxis()
    # p1.plot([x0,x1],[0.0,0.0],"--",color="k",linewidth=2)
    # p1.plot([0.0,0.0],[y0,y1],"--",color="k",linewidth=2)
    # p1.set_xlim(x0,x1)
    # p1.set_ylim(y0,y1)
    # p1.set_xlabel("RA (deg) [$\sigma$=%.3f]"%raSig,fontsize=10)
    # p1.set_ylabel("DEC (deg) [$\sigma$=%.3f]"%decSig,fontsize=10)
    
    p1.set_xlabel("RA (deg)" %raSig,fontsize=10)
    p1.set_ylabel("DEC (deg)"%decSig,fontsize=10)
     
    # p2.set_xticks([])
    # p2.set_yticks([])
    # p3.set_xticks([])
    # p3.set_yticks([])
    #p2.set_xlim(x0,x1)
    #p2.set_ylim(y0,y1)
    bin_w = 0.005
    xy_margin =np.max([np.max(np.abs(x0)),np.max(np.fabs(y0))])
    lim = int(xy_margin/bin_w+1)*bin_w
    p2.set_xlim(-lim,lim)
    #p2.set_ylabel("persent",fontsize=10)
    p3.set_xlim(-lim,lim)
    #p3.set_xlabel("persent",fontsize=10)
    bins = np.arange(-lim, lim+bin_w,bin_w)
    p2.hist(dra,bins=bins)
    p3.hist(ddec, bins=bins)#,orientation='horizontal')
    #pl.suptitle("GridSpec Inside GridSpec",color='r')

    
    p2.set_xlabel("$\Delta$RA (arcsec) [$\sigma$=%.3f, $\mu$=%.3f ]"%(raSig,raMu),fontsize=10)
    p3.set_xlabel("$\Delta$DEC (arcsec) [$\sigma$=%.3f, $\mu$=%.3f ]"%(decSig,decMu),fontsize=10)
     
    #  p1.set_xlabel('RA(deg), $\mu$='+str(raMu)+', $\sigma$='+str(raSig),fontsize=10)
    # p1.set_ylabel('DEC(deg), $\mu$='+str(decMu)+', $\sigma$='+str(decSig),fontsize=10)
     
    pl.savefig(figout,dpi=400)
    pl.savefig(mini_figout,dpi=50)
    logger.info(figout)
    if not os.path.exists(figout):
        logger.info(figout + 'can not generateXXXXXXXXXXXX')
    pl.clf()
    pl.close()
    return

 

def ewhisker(ximg, yimg, aimg, bimg, theta, figout="ewhisker.pdf"):
    """
    Plot Whisker map of PSF ellipticity

    Parameters:
    catn: string
        LDAC catalog name
    figout: string
        output PSF Whisker map
    """
    ## plot the Whisker plot
    ell  = (aimg-bimg)/(aimg+bimg)
    mask  = sigma_clip(ell, sigma=3.0, maxiters=3, masked=True)
    ell = ell[~mask.mask]
    ximg, yimg, theta = ximg[~mask.mask], yimg[~mask.mask], theta[~mask.mask]

    ell1 = np.abs(ell) * np.cos(2.0*theta*np.pi/180.0)
    ell2 = np.abs(ell) * np.sin(2.0*theta*np.pi/180.0)

    eRef = float("%5.3f"%np.abs(np.mean(ell)))
    pl.quiver(ximg,yimg,ell1,ell2,angles=theta,width=0.002,scale=0.8,
                          pivot="middle",headwidth=0, headlength=4.5)
    pl.quiver(ximg.max()-200,yimg.max()-300,eRef,0.0,angles=0.0,color="r",width=0.002,scale=0.8,
                          pivot="middle",headwidth=0, headlength=4.5)
    pl.text(ximg.max()-400,yimg.max()-200,"e=%.3f"%eRef, color="r", fontsize=12)
    pl.xlim([ximg.min(),ximg.max()])
    pl.ylim([yimg.min(),yimg.max()])
    pl.xlabel("X (pixel)")
    pl.ylabel("Y (pixel)")
    pl.savefig(figout)
    pl.clf()
    pl.close()
    return

def ewhiskerFWHM(ximg, yimg, fwhm, theta, scale=10, figout="ewhiskerFWHM.pdf"):
    """
    Plot Whisker map of FWHM

    Parameters:
    catn: string
        LDAC catalog name
    figout: string
        output PSF Whisker map
    """
    ## plot the Whisker plot
    mask  = sigma_clip(fwhm, sigma=3.0, maxiters=3, masked=True)
    ximg, yimg, fwhm, theta = ximg[~mask.mask], yimg[~mask.mask], fwhm[~mask.mask], theta[~mask.mask]
    
    fwhmF = 0.5*fwhm*scale
    nobj  = len(fwhm)
    x1   = ximg + fwhmF*np.cos(theta*np.pi/180.0)
    x2   = ximg - fwhmF*np.cos(theta*np.pi/180.0)
    y1   = yimg + fwhmF*np.sin(theta*np.pi/180.0)
    y2   = yimg - fwhmF*np.sin(theta*np.pi/180.0)
    for i in range(nobj):
        pl.plot([x1[i],x2[i]], [y1[i],y2[i]],"k-",linewidth=1.5)

    eFWHM = 5.0 * scale
    pl.plot([ximg.max()-400,ximg.max()-400+eFWHM],[yimg.max()-300,yimg.max()-300],"r-",linewidth=3.0)
    pl.text(ximg.max()-500, yimg.max()-200,"FWHM=5.0", color="r", fontsize=12)
    pl.xlim([ximg.min(),ximg.max()])
    pl.ylim([yimg.min(),yimg.max()])
    pl.xlabel("X (pixel)")
    pl.ylabel("Y (pixel)")
    pl.savefig(figout)
    pl.clf()
    pl.close()
    return

def ewhiskerEllipse(ximg, yimg, aimg, bimg, theta, scale=10.0, figout="ewhiskerEllipse.pdf"):
    """
    show the ellipse
    """ 
    fig, ax = pl.subplots(subplot_kw={'aspect': 'equal'})
    nobj = len(ximg)
    for i in range(nobj):
        iell = Ellipse(xy=np.array([ximg[i],yimg[i]]), width=aimg[i]*2*scale, height=bimg[i]*2*scale, angle=theta[i])
        ax.add_artist(iell)
        iell.set_clip_box(ax.bbox)
        iell.set_alpha(np.random.rand())
        iell.set_facecolor(np.random.rand(3))

    ax.set_xlim(ximg.min(),ximg.max())
    ax.set_ylim(yimg.min(),yimg.max())

    pl.savefig(figout)
    pl.clf()
    pl.close()
    return

def zpPlot(zpArray, magArray, figout="zpCom.pdf"):
    """
    Show the zeropoint variation as a function of magnitude
    """
    zpMed    = np.median(zpArray)
    zpStd    = np.std(zpArray)
    magArray = magArray+zpMed
    xlim     = [magArray.min()-0.5, magArray.max()+0.5]
    pl.scatter(magArray, zpArray, color="black", marker="o", s=6)
    pl.plot(xlim, [zpMed, zpMed],"r-",linewidth=2.0)
    pl.plot(xlim, [zpMed-zpStd, zpMed-zpStd],"r--",linewidth=1.5)
    pl.plot(xlim, [zpMed+zpStd, zpMed+zpStd],"r--",linewidth=1.5)
    pl.xlim(xlim)
    pl.xlabel("PSF_MAG", fontsize=18)
    pl.ylabel("ZPoint", fontsize=18)
    pl.title("ZP = %.3f $\pm$ %.3f (#%d stars)"%(zpMed,zpStd,len(zpArray)),fontsize=15)
    pl.savefig(figout)
    pl.clf()
    pl.close()
    return

def mag2Err1(flux, fluxErr, fluxCom, fluxErrCom, bkgNoise=100.0, zpoint=0.0, gain=1.0, figout="magErr.pdf"):
    gid      = flux>0.0
    flux, fluxErr = flux[gid], fluxErr[gid]
    mag      = -2.5*np.log10(flux) + zpoint
    magErr   = 1.0857 * (fluxErr/flux)

    magCom      = -2.5*np.log10(fluxCom) + zpoint
    magErrCom   = 1.0857 * (fluxErrCom/fluxCom)    

    magBin   = np.linspace(mag.min(), mag.max(),200) - zpoint
    fluxX    = 10**(-0.4*(magBin))
    pshotErr = 1.0857/(np.sqrt(fluxX*gain))
    bkgErr   = 1.0857*bkgNoise/fluxX

    pl.plot(magCom, magErrCom,"ro",  markersize=2.0, alpha=0.5,label="AUTO Mag")
    pl.plot(mag,    magErr,   "ko",  markersize=2.0, label="PSF Mag")
    pl.plot(magBin+zpoint, pshotErr, "m--", linewidth=2.0, label="Shot Noise")
    pl.plot(magBin+zpoint, bkgErr,   "b--", linewidth=2.0, label="Background Noise")
    pl.xlabel("PSF_MAG", fontsize=18)
    pl.ylabel("PSF_MAGERR", fontsize=18)
    pl.legend(frameon=False, loc="upper left", fontsize=12)
    pl.xlim([mag.min()-0.5,mag.max()+0.5])
    pl.savefig(figout)
    pl.clf()
    pl.close()
    return


def mag2Err(flux, fluxErr, fluxCom, fluxErrCom, bkgNoise=100.0, zpoint=0.0, gain=1.0, figout="magErr.pdf"):
    # 202309081705 dugking
    logger.info(f"{flux, fluxErr, fluxCom, fluxErrCom, bkgNoise, zpoint, gain, figout}")
    gid      = flux>0.0
    flux, fluxErr = flux[gid], fluxErr[gid]
    mag      = -2.5*np.log10(flux) + zpoint
    magErr   = 1.0857 * (fluxErr/flux)
    magCom      = -2.5*np.log10(fluxCom) + zpoint
    magErrCom   = 1.0857 * (fluxErrCom/fluxCom)
    # 202309081705 dugking
    logger.info(f"mag:{mag}")
    magBin   = np.linspace(mag.min(), mag.max(),200) - zpoint
    fluxX    = 10**(-0.4*(magBin))
    pshotErr = 1.0857/(np.sqrt(fluxX*gain))
    bkgErr   = 1.0857*bkgNoise/fluxX
    pl.plot(magCom, magErrCom,"ro",  markersize=2.0, alpha=0.5,label="AUTO Mag")
    pl.plot(mag,    magErr,   "ko",  markersize=2.0, label="PSF Mag")
    pl.plot(magBin+zpoint, pshotErr, "m--", linewidth=2.0, label="Shot Noise")
    pl.plot(magBin+zpoint, bkgErr,   "b--", linewidth=2.0, label="Background Noise")
    pl.xlabel("PSF_MAG", fontsize=18)
    pl.ylabel("PSF_MAGERR", fontsize=18)
    pl.legend(frameon=False, loc="upper left", fontsize=12)
    pl.xlim([mag.min()-0.5,mag.max()+0.5])
    pl.savefig(figout)
    pl.clf()
    pl.close()
    return



