import os
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy import coordinates as coord, units as u
import datetime
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
# from PyAstronomy import pyas
import glob
import lib.phot.cal_airmass as cal_airmass
from loguru import logger

from lib.phot.ybias import get_current_dir




def plotshow(x,y,ye,ylabel,title,filename):
    plt.switch_backend('agg')
    plt.errorbar(x, y, yerr=ye,fmt='o',ecolor='r',color='b',elinewidth=2,capsize=4)
    #plt.xlim(0, 0.7)
    plt.xlabel('MJD-2459334.3(days)')
    plt.ylabel(ylabel)
    plt.title(title)
    #plt.clf()
    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(filename)


#def TransLDAC_time(ldacfilename):
#    data=fits.getdata(ldacfilename,ext=1)
#    obstime=np.array(data[0][0])[5]
#    exptime=np.array(data[0][0])[7]
#    obstime=obstime.split('/')[0].split('=')[1].replace('\'',"").strip().replace('T',' ')#.strip('\'')
#    exptime=exptime.split('/')[0].split('=')[1]
#    obstime=datetime.datetime.strptime(obstime,'%Y-%m-%d %H:%M:%S')
    # 计算偏移量
#    offset = datetime.timedelta(seconds=float(exptime)/2)
    # 获取修改后的时间并格式化
#    timeplus = obstime + offset
#    timeplus=str(timeplus).replace(' ','T')
#    t=Time(timeplus)
#    mjd=float(str(t.jd))
#    mjd=Time(mjd,format='mjd')
#    return obstime,float(exptime),timeplus,mjd


def airmass(tra,tdec,obstime):
    obstimeut=Time(obstime,format="datetime",scale='utc')
    gmg=coord.EarthLocation(100.03*u.deg,26.7089*u.deg,height=3200*u.m)#100.8667*u.deg,26.7089 东经100 °01′51″，北纬26 °42′32″，海拔3200米
    objectstar=coord.SkyCoord(tra*u.degree,tdec*u.degree,frame='icrs')
    altaz_obj =objectstar.transform_to(AltAz(obstime=obstimeut,location=gmg))
    #altaz_obj =objectstar.transform_to(AltAz(obstime=MJDTrans(mjdtime,30),location=gmg))
    airmass=altaz_obj.secz
    return float(airmass)


def ams_cal(filename,ira,idec):
    res=fits.open(filename)
    data = res[1].data
    date_obs = data.field(0)[0,5][11:30]
    expt     = float(data.field(0)[0,6][12:].split(' ')[0])
    amss,jd_mid,tbd_mid = cal_airmass.cam(date_obs,expt,ira,idec)
    return amss


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




def TransLDAC_time(ldacfilename):
    data=fits.getdata(ldacfilename,ext=1)
    
    obstime=np.array(data[0][0])[5]
    
    exptime=np.array(data[0][0])[7]
    
    #logger.info(data[0][0])
    obstime=obstime.split('/')[0].split('=')[1].replace('\'',"").strip().replace('T',' ')#.strip('\'')
    exptime=exptime.split('/')[0].split('=')[1]
    obstime=datetime.datetime.strptime(obstime,'%Y-%m-%d %H:%M:%S')
    # 计算偏移量
    offset = datetime.timedelta(seconds=float(exptime)/2)
    # 获取修改后的时间并格式化
    timeplus = obstime + offset
    logger.info(timeplus)
    #timeplus=str(timeplus).replace(' ','T')
    t=Time(timeplus)
    tjd=Time(np.array(str(t.jd)),format='jd')
    tb=tjd.tdb
    return timeplus,str(tb),str(tjd)


def timeformat(stime):
    a=str(stime[:-9])
    b=str(stime[11:])
    st = a+' '+b
    st=datetime.datetime.strptime(st,'%Y-%m-%d %H:%M:%S')
    return st
    


def delta_time(t0,ti):
    #t0=timelist[0]
    #delta_timelist=[]
    t0=timeformat(t0)
    dt=timeformat(ti)-t0
    ddt=str(dt.total_seconds())
    
    #logger.info('t0,ti and the ddt is:', t0, timeformat(ti),ddt)
    dddt=(float(ddt)/60.0/39.8876)%1   
    return dddt


def phasetime(timelist):
    t0=timelist[0][0]
    logger.info(timelist[0][0])
    deltalist=[]
    for i in range(0,len(timelist)):
        dt=delta_time(t0,timelist[i][0])
        deltalist.append(dt)
    return deltalist





def MJDTrans(objtime,exptime):
    #logger.info('-----------------------',objtime)
    t=Time(objtime)
    tjd=Time(np.array(str(t.jd)),format='jd')
    mjd1=float(str(t.jd))
    mjd2=mjd1+(exptime/24/60/60)/2
    return mjd2



def getindex_location(tra,tdec,ralist,declist):
    for i in range(len(ralist)):
        d=np.sqrt((ralist[i]-tra)**2+(declist[i]-tdec)**2)
        #logger.info(ra[i],',',dec[i],',d=',d)
        if(d<8/3600):
            return i
    #return 0
            





def airms_plot(rootpath,ttfname,date,mmag,tra,tdec):
    #tra=69.0428297
    #tdec= 19.87394052
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'+ttfname+'/'
    #tra,tdec=d2hms('21:30:57','44:20:21',conv=1)
    logger.info(tra,tdec)
    fblist = glob.glob(tscidir +'*sexcat.fits')
    ldacList=np.sort(fblist)
    logger.info(fblist)

    amx=[]
    magy_psf=[]
    magerry_psf=[]

    for i in range(len(ldacList)):
        ildac = fits.getdata(ldacList[i],ext=2)
        logger.info(ldacList[i][:-21]+'.fits')
        ztmag=fits.getheader(ldacList[i][:-21]+'.fits')['ZPMAG']
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        mag_aper,magerr_aper=ildac["MAG_APER"],ildac["MAGERR_APER"]
        mag_psf,magerr_psf=ildac["MAG_PSF"],ildac["MAGERR_PSF"]
        index=getindex_location(tra,tdec,ra, dec)
        try:
           utctime,tdbtime,jd=TransLDAC_time(ldacList[i])
           #logger.info('utctime,tdbtime,jd=',str(utctime),tdbtime,jd)
        except:
           logger.info('the header is wrong')
           continue
        #data_array.append([utctime,tdbtime,psfmag,psfarr])
        logger.info('#########################'+ldacList[i])
        #mjd1=MJDTrans(time,float(exptime)/2)
        amx.append(airmass(tra,tdec,utctime))
        logger.info('utctime,mag,airmass=',str(utctime),str(mag_psf[index]),str(airmass(tra,tdec,utctime)))
        magy_psf.append(mag_psf[index])#+ztmag)
        magerry_psf.append(magerr_psf[index])
        #logger.info(tdbtime,mag_psf[index]+ztmag,magerr_psf[index])

    plt.switch_backend('agg')
    plt.errorbar(amx, magy_psf, yerr=magerry_psf,fmt='o',ecolor='r',color='r',elinewidth=1,capsize=2,ms=3)
    #plt.plot(mjdx, magy_psf)
    plt.xlabel('airmass')
    plt.ylabel('M-MI(mag)')
    plt.title(ttfname+'_M='+str(mmag))
    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(rootpath+date+'_'+ttfname+'_airmass_mag.pdf')
    logger.info('ok')



def mag_plot_2(rootpath,ttfname1,ttfname2,date,tra,tdec):
    tra=69.0428297
    tdec= 19.87394052
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir1 = upath+'sci/'+ttfname1+'/'
    #tra,tdec=d2hms('21:30:57','44:20:21',conv=1)
    logger.info(tra,tdec)
    fblist = glob.glob(tscidir1 +'*sexcat.fits')
    ldacList=np.sort(fblist)
    logger.info(fblist)
    mjdx=[]
    magy_psf=[]
    magerry_psf=[]
    airmy=[]
    for i in range(len(ldacList)):
        ildac = fits.getdata(ldacList[i],ext=2)
        logger.info(ldacList[i][:-21]+'.fits')
        ztmag=fits.getheader(ldacList[i][:-21]+'.fits')['ZPMAG']
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        mag_auto, magerr_auto=ildac["MAG_AUTO"],ildac["MAGERR_AUTO"]
        mag_aper,magerr_aper=ildac["MAG_APER"],ildac["MAGERR_APER"]
        mag_psf,magerr_psf=ildac["MAG_PSF"],ildac["MAGERR_PSF"]
        index=getindex_location(tra,tdec,ra, dec)
        try:
           utctime,tdbtime,jd=TransLDAC_time(ldacList[i])
        except:
           logger.info('the header is wrong')
           continue
        #data_array.append([utctime,tdbtime,psfmag,psfarr])
        #logger.info('#########################'+ldacList[i])
        #mjd1=MJDTrans(time,float(exptime)/2)
        #airmy.append(airmass(tra,tdec,utctime))
        mjdx.append(float(tdbtime))
        magy_psf.append(mag_psf[index]-0.1)#+ztmag)
        magerry_psf.append(magerr_psf[index])
        #logger.info(tdbtime,mag_psf[index]+ztmag,magerr_psf[index])
    

    tscidir2 = upath+'sci/'+ttfname2+'/'
    fblist2 = glob.glob(tscidir2 +'*sexcat.fits')
    ldacList2=np.sort(fblist2)
    logger.info(fblist2)
    mjdx2=[]
    magy_psf2=[]
    magerry_psf2=[]
    airmy2=[]
    for i in range(len(ldacList2)):
        ildac2 = fits.getdata(ldacList2[i],ext=2)
        logger.info(ldacList2[i][:-21]+'.fits')
        ztmag2=fits.getheader(ldacList2[i][:-21]+'.fits')['ZPMAG']
        ra, dec   = ildac2["ALPHA_J2000"], ildac2["DELTA_J2000"]
        mag_auto, magerr_auto=ildac2["MAG_AUTO"],ildac2["MAGERR_AUTO"]
        mag_aper,magerr_aper=ildac2["MAG_APER"],ildac2["MAGERR_APER"]
        mag_psf,magerr_psf=ildac2["MAG_PSF"],ildac2["MAGERR_PSF"]
        index2=getindex_location(tra,tdec,ra, dec)
        try:
           utctime,tdbtime,jd=TransLDAC_time(ldacList2[i])
        except:
           logger.info('the header is wrong')
           continue
        #data_array.append([utctime,tdbtime,psfmag,psfarr])
        logger.info('#########################'+ldacList2[i])
        #mjd1=MJDTrans(time,float(exptime)/2)
        mjdx2.append(float(tdbtime))
        magy_psf2.append(mag_psf[index2])#+ztmag)
        magerry_psf2.append(magerr_psf[index2])
        #logger.info(tdbtime,mag_psf2[index]+ztmag2,magerr_psf2[index])
    plt.switch_backend('agg')

    plt.errorbar(mjdx, magy_psf, yerr=magerry_psf,fmt='o',ecolor='r',color='r',elinewidth=1,capsize=2,ms=3)
    plt.errorbar(mjdx2, magy_psf2, yerr=magerry_psf2,fmt='o',ecolor='b',color='b',elinewidth=1,capsize=2,ms=3)
   
    #plt.plot(mjdx, magy_psf)
    plt.xlabel('tdb')
    plt.ylabel('PSF_APER(mag)')
    plt.legend(['g-0.1','r'])
    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(rootpath+date+'_'+ttfname1+'_'+ttfname2+'_time_mag.pdf')
    logger.info('ok')








def mag_plot(rootpath,ttfname,date,tra,tdec):
    #tra=322.7362855
    #tdec=44.346236
    #cra=322.696516
    #cdec=44.335383
    #tra=69.0428297 
    #tdec= 19.87394052
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'+ttfname+'/'
    #tra,tdec=d2hms('21:30:57','44:20:21',conv=1)
    logger.info(tra,tdec)
    fblist = glob.glob(tscidir +'*sexcat.fits')
    ldacList=np.sort(fblist)
    logger.info(fblist)
    
    mjdx=[]
   

    magy_psf=[]
    magerry_psf=[]

    for i in range(len(ldacList)):
        ildac = fits.getdata(ldacList[i],ext=2)
        logger.info(ldacList[i][:-21]+'.fits')
        ztmag=fits.getheader(ldacList[i][:-21]+'.fits')['ZPMAG']
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        mag_auto, magerr_auto=ildac["MAG_AUTO"],ildac["MAGERR_AUTO"]
        mag_aper,magerr_aper=ildac["MAG_APER"],ildac["MAGERR_APER"]
        mag_psf,magerr_psf=ildac["MAG_PSF"],ildac["MAGERR_PSF"]
        index=getindex_location(tra,tdec,ra, dec)
           
          
        try:
           utctime,tdbtime,jd=TransLDAC_time(ldacList[i])
        except:
           logger.info('the header is wrong')
           continue
        #data_array.append([utctime,tdbtime,psfmag,psfarr])
        logger.info('#########################'+ldacList[i])
        #mjd1=MJDTrans(time,float(exptime)/2)
        mjdx.append(float(tdbtime))
        magy_psf.append(mag_psf[index])#+ztmag)
        magerry_psf.append(magerr_psf[index])
        logger.info(tdbtime,mag_psf[index],magerr_psf[index])
  
    plt.switch_backend('agg')

    plt.errorbar(mjdx, magy_psf, yerr=magerry_psf,fmt='o',ecolor='r',color='r',elinewidth=1,capsize=2,ms=3)
    #plt.plot(mjdx, magy_psf)
    plt.xlabel('tdb')
    plt.ylabel('PSF_APER(mag)')
    plt.title(ttfname)
    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(rootpath+date+'_'+ttfname+'_time_apermag.pdf')
    logger.info('ok')





def phase_plot(rootpath,ttfname_b,ttfname_v,date,tra,tdec):
    #tra=322.7362855
    #tdec=44.346236
    #cra=322.696516
    #cdec=44.335383
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir_b = upath+'sci/'+ttfname_b+'/'
    tscidir_v = upath+'sci/'+ttfname_v+'/'
    #tra,tdec=d2hms('21:30:57','44:20:21',conv=1)
    logger.info(tra,tdec)
    ldacListn_b   = tscidir_b +ttfname_b+"_o0_ldac.list"
    ldacListn_v   = tscidir_v +ttfname_v+"_o0_ldac.list"
    ldacList_b = open(ldacListn_b,"r").read().splitlines()
    ldacList_v = open(ldacListn_v,"r").read().splitlines()
    jlist_b=np.load(tscidir_b+ttfname_b+'_time_o0.npy')
    jlist_v=np.load(tscidir_v+ttfname_v+'_time_o0.npy')
    t0b=jlist_b[0][0]
    t0v=jlist_v[0][0]
    plist_b=[]#phasetime(jlist_b)
    plist_v=[]#phasetime(jlist_v)
    lenfile=min(len(ldacListn_b),len(ldacListn_v))
    b_magy_psf=[]
    b_magerry_psf=[]
    v_magy_psf=[]
    v_magerry_psf=[]
    for i in range(0,lenfile):
        dt_b=delta_time(t0b,jlist_b[i][0])
        dt_v=delta_time(t0v,jlist_v[i][0])
        plist_b.append(dt_b)
        plist_v.append(dt_v)

        ildac_b= fits.getdata(ldacList_b[i],ext=2)
        ztmag_b=fits.getheader(ldacList_b[i][:-14]+date+'_subbkg_o0.fits')['ZPMAG']
        rab, decb   = ildac_b["ALPHA_J2000"], ildac_b["DELTA_J2000"]
        mag_psf_b,magerr_psf_b=ildac_b["MAG_PSF"],ildac_b["MAGERR_PSF"]
        index_b=getindex_location(tra,tdec,rab, decb)
        
        b_magy_psf.append(mag_psf_b[index_b]+ztmag_b)
        b_magerry_psf.append(magerr_psf_b[index_b])
        
        ildac_v= fits.getdata(ldacList_v[i],ext=2)
        ztmag_v=fits.getheader(ldacList_v[i][:-14]+date+'_subbkg_o0.fits')['ZPMAG']
        rav, decv   = ildac_v["ALPHA_J2000"], ildac_v["DELTA_J2000"]
        mag_psf_v,magerr_psf_v=ildac_v["MAG_PSF"],ildac_v["MAGERR_PSF"]
        index_v=getindex_location(tra,tdec,rav, decv)
        
        v_magy_psf.append(mag_psf_v[index_v]+ztmag_v)
        v_magerry_psf.append(magerr_psf_v[index_v])
    

    for i in range(0,lenfile):
         plist_b.append(plist_b[i]+1)
         plist_v.append(plist_v[i]+1)      
    v_magy_psf=list(v_magy_psf)+(list(v_magy_psf))
    v_magerry_psf=list(v_magerry_psf)+(list(v_magerry_psf))
    b_magy_psf=list(b_magy_psf)+(list(b_magy_psf))
    b_magerry_psf=list(b_magerry_psf)+(list(b_magerry_psf))

    plt.switch_backend('agg')
    logger.info('plist_b:',np.shape(plist_b),' b_magy_psf:',np.shape(b_magy_psf))
    plt.errorbar(np.array(plist_b)-0.21, b_magy_psf, yerr=b_magerry_psf,fmt='o',ecolor='r',color='r',elinewidth=1,capsize=2,ms=3)
    plt.errorbar(np.array(plist_v)-0.21, np.array(v_magy_psf)+0.4, yerr=v_magerry_psf,fmt='o',ecolor='black',color='black',elinewidth=1,capsize=2,ms=3)
    plt.legend(['B','V'])
    plt.xlabel('phase')
    plt.ylabel('PSF_MAG(mag)')
    #plt.title(title)
    #plt.clf()
    #plt.xlim(-0.75,1)
    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(rootpath+date+'_'+date+'_bvphase_mag.pdf')
    logger.info('ok')
 


            
def plot(rootpath,ttfname,date,tra,tdec):
    tra=322.7362855
    tdec=44.346236
    #cra=322.696516
    #cdec=44.335383
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'+ttfname+'/'
    #tra,tdec=d2hms('21:30:57','44:20:21',conv=1)
    logger.info(tra,tdec)
    ldacListn   = tscidir +ttfname+"_o0_ldac.list"
    logger.info(ldacListn)
    ldacList = open(ldacListn,"r").read().splitlines()
    jlist=np.load(tscidir+ttfname+'_time_o0.npy')
    #mjd1=float(str(jd)[1:])
    #mjd2=mjd1+(exptime/24/60/60)/2
    #mag, magerr=ldacList["MAG_AUTO"],ldacList["MAGERR_AUTO"]
    #3.0,3.5, 4.0,4.5, 5.0,5.5, 6.0,6.5, 7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,16.0,17.0,18.0,19.0,20.0
    mjdx=[]
    magy=[]
    magerry=[]
    magy_auto=[]
    magerry_auto=[]
    magy_psf=[]
    magerry_psf=[]
    magy_aper=[]
    magerry_aper=[]
    for i in range(len(ldacList)):
        ildac = fits.getdata(ldacList[i],ext=2)
        ztmag=fits.getheader(ldacList[i][:-14]+date+'_subbkg_o0.fits')['ZPMAG']
        logger.info('ztmag=',ztmag)
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        mag_auto, magerr_auto=ildac["MAG_AUTO"],ildac["MAGERR_AUTO"]
        mag_aper,magerr_aper=ildac["MAG_APER"],ildac["MAGERR_APER"]
        mag_psf,magerr_psf=ildac["MAG_PSF"],ildac["MAGERR_PSF"]
        index=getindex_location(tra,tdec,ra, dec)
        time,exptime=jlist[i]
        mjd1=MJDTrans(time,float(exptime)/2)
        mjdx.append(mjd1)
        magy_auto.append(mag_auto[index]+ztmag)
        magerry_auto.append(magerr_auto[index])
        magy_psf.append(mag_psf[index]+ztmag)
        magerry_psf.append(magerr_psf[index])
        magy_aper.append(mag_aper[index][5]+ztmag)
        magerry_aper.append(magerr_aper[index][5])

        #logger.info(mjd2)
    logger.info('len of mjdx,magy,magerry:',len(mjdx),len(magy),len(magerry))
    logger.info('shape of mjdx, magy,magerry:',np.shape(mjdx),np.shape(magy),np.shape(magerry))
    figname_auto=tscidir+date+'_'+ttfname+'_auto_mag.pdf'
    figname_psf=tscidir+date+'_'+ttfname+'_psf_mag.pdf'
    figname_aper=tscidir+date+'_'+ttfname+'_aper_mag.pdf'
    plotshow(mjdx,np.array(magy_auto),np.array(magerry_auto),'AUTO_MAG(mag)',ttfname,figname_auto)
    plotshow(mjdx,np.array(magy_psf),np.array(magerry_psf),'PSF_MAG(mag)',ttfname,figname_psf)
    plotshow(mjdx,np.array(magy_aper),np.array(magerry_aper),'APER_MAG(mag)',ttfname,figname_aper)
    #plotshow(x,y,ye,ylabel,title,filename)
    #listh=[mjdx,magy,magerry]
    #np.save(tscidir+ttfname+'_mag.npy',listh)

def datacal_gate(rootpath,ttfname,date,tra,tdec,ccdindex):
    tra=322.7362855
    tdec=44.346236
    upath=rootpath+'reception/'+str(date)+'/'
    #tscidir = upath+'sci/fileguide/'
    tscidir =upath+'sci/'+ttfname+'/'
    #logger.info(tscidir)
    #tra,tdec=d2hms('21:30:57','44:20:21',conv=1)
     
    ldacListn   = tscidir +ttfname+"_o0_ldac.list"
    ldacList = open(ldacListn,"r").read().splitlines()
    #logger.info('2+++++++++++++++'+ldacListn)
    jlist=np.load(tscidir+ttfname+'_time_o0.npy')
    #mjd1=float(str(jd)[1:])
    #mjd2=mjd1+(exptime/24/60/60)/2
    #mag, magerr=ldacList["MAG_AUTO"],ldacList["MAGERR_AUTO"]
    mjdx=[]
    magy=[]
    magerry=[]
    #logger.info(jlist)
    for i in range(len(ldacList)):
        #ldacname=ldacList[i][:-15]+'_subbkg_o'+str(ccdindex)+'.ldac'
        ldacname=ldacList[i][:-15]+'_subbkg_o0.ldac'
        #ldacname=ldacList[i][:-15]+'_subbkg_o0_20210429_sexcat.fits'
        ildac = fits.getdata(ldacname,ext=2)
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        mag, magerr=ildac["MAG_AUTO"],ildac["MAGERR_AUTO"]
        index=getindex_location(tra,tdec,ra, dec)
        time,exptime=jlist[i]
        mjd1=MJDTrans(time,float(exptime)/2)
        mjdx.append(mjd1)
        magy.append(mag[index])
        logger.info('the'+str(i)+'th source index:'+str(index))
        magerry.append(magerr[index])
    return mjdx,magy,magerry

def plot_mag_combine(rootpath,ttfname,date,tra,tdec,ccdindex):
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'+ttfname+'/'
    mjdx,magy,magerry= datacal_gate(rootpath,ttfname,date,tra,tdec,ccdindex)
    logger.info('len of mjdx,magy,magerry:',len(mjdx),len(magy),len(magerry))
    #logger.info('shape of mjdx, magy,magerry:',np.shape(np.array(mjdx)),np.shape(np.array(magy)),np.shape(np.array(magerry)))
    #logger.info('shape of mjdx, magy,magerry:',np.shape(np.array(magy)))
    logger.info(magy)
    figname=tscidir+ttfname+'_mag.pdf'
    plotshow(mjdx,np.array(magy),np.array(magerry),ttfname,figname)
    listh=[mjdx,magy,magerry]
    np.save(tscidir+ttfname+'_mag.npy',listh)


def datacal_combine(rootpath,ttfname,date,tra,tdec):
    tra=322.7362855
    tdec=44.346236
    
    rra=322.73015
    rdec=44.3311416

    cra=322.6965166
    cdec=44.3353833

    upath=rootpath+'reception/'+str(date)+'/'
    #tscidir = upath+'sci/fileguide/'
    tscidir =upath+'sci/'+ttfname+'/'
    #logger.info(tscidir)
    #tra,tdec=d2hms('21:30:57','44:20:21',conv=1)

    ldacListn   = tscidir +ttfname+"_o0_ldac.list"
    ldacList = open(ldacListn,"r").read().splitlines()
    #logger.info('2+++++++++++++++'+ldacListn)
    jlist=np.load(tscidir+ttfname+'_time_o0.npy')
    #mjd1=float(str(jd)[1:])
    #mjd2=mjd1+(exptime/24/60/60)/2
    #mag, magerr=ldacList["MAG_AUTO"],ldacList["MAGERR_AUTO"]
    mjdx=[]
    magy_r=[]
    magy_c=[]
    magy_t=[]
    magerry=[]
    #logger.info(jlist)
    for i in range(len(ldacList)):
        #ldacname=ldacList[i][:-15]+'_subbkg_o'+str(ccdindex)+'.ldac'
        ldacname=ldacList[i][:-15]+'_subbkg_o0.ldac'
        #ldacname=ldacList[i][:-15]+'_subbkg_o0_20210429_sexcat.fits'
        ildac = fits.getdata(ldacname,ext=2)
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        mag, magerr=ildac["MAG_AUTO"],ildac["MAGERR_AUTO"]
        index_t=getindex_location(tra,tdec,ra, dec)
        index_c=getindex_location(cra,cdec,ra, dec)
        index_r=getindex_location(rra,rdec,ra, dec)
        time,exptime=jlist[i]
        mjd1=MJDTrans(time,float(exptime)/2)
        mjdx.append(mjd1)
        magy_r.append(mag[index_r])
        magy_c.append(mag[index_c])
        magy_t.append(mag[index_t])
        logger.info('the'+str(i)+'th source index_r='+str(index_r)+' index_t:'+str(index_t)+' index_c:'+str(index_c))
       
        #magerry.append(magerr[index])
    return mjdx,magy_t, magy_r, magy_c

def plotshow_combine(x,ty,ry,cy,tffname,filename):
    plt.switch_backend('agg')
    plt.scatter(x, ty-ry, c='b')
    plt.scatter(x,cy-ry,c='c')
    #plt.xlim(0, 0.7)
    plt.xlabel('MJD-59334.3(days)')
    filterid=tffname[-1]
    plt.ylabel(filterid+'_delta_mag(mag)')


    #plt.clf()
    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(filename)


def plot_combine(rootpath,ttfname,date,tra,tdec):
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'+ttfname+'/'
    mjdx,magy_t, magy_r, magy_c= datacal_combine(rootpath,ttfname,date,tra,tdec)
    logger.info('len of mjdx,magy,magerry:',len(mjdx),len(magy_t),len(magy_c),len(magy_r))
    #logger.info('shape of mjdx, magy,magerry:',np.shape(np.array(mjdx)),np.shape(np.array(magy)),np.shape(np.array(magerry)))
    #logger.info('shape of mjdx, magy,magerry:',np.shape(np.array(magy)))
    
    figname=tscidir+ttfname+'_diff_mag.pdf'
    plotshow_combine(mjdx,np.array(magy_t),np.array(magy_r),np.array(magy_c),ttfname,figname)
    #listh=[mjdx,magy,magerry]
    #np.save(tscidir+ttfname+'diff_mag.npy',listh)



def getindex_loc(tra,tdec,ralist,declist):
    #tmpval = 999
    #tmp= 1.5
    for i in range(len(ralist)):
        d=np.sqrt((ralist[i]-tra)**2+(declist[i]-tdec)**2)
        #logger.info(ra[i],',',dec[i],',d=',d)
        if(d<8/3600):
            #if ( d < tmpval):
            #    tmpval = d
            #    tmp  = i
            return i ,ralist[i],declist[i]
       
    logger.info('can not find the location')
          

def ams_plot(rootpath,ttfname,date,tra,tdec):
    upath=rootpath+'reception/'+str(date)+'/'
    tscidir = upath+'sci/'+ttfname+'/'
    #tra,tdec=d2hms('21:30:57','44:20:21',conv=1)
    #ldacListn   = tscidir +ttfname+"_o0_ldac.list"
    #ldacList = open(ldacListn,"r").read().splitlines()
    #logger.info(tscidir)
    fblist = glob.glob(tscidir +'*_sexcat.fits')
    ldacList=np.sort(fblist)
    #logger.info(ldacList)
    data_array=[]
    amsx=[]
    magy=[]
    magyerr=[]
    logger.info(str(len(ldacList))+'################################## mags of target in this file!')
    for i in range(len(ldacList)):
        ildac = fits.getdata(ldacList[i],ext=2)
        logger.info(ldacList[i][:-21]+'.fits')
        ztmag=fits.getheader(ldacList[i][:-21]+'.fits')['ZPMAG']
        ra, dec   = ildac["ALPHA_J2000"], ildac["DELTA_J2000"]
        mag_auto, magerr_auto=ildac["MAG_AUTO"],ildac["MAGERR_AUTO"]
        mag_aper,magerr_aper=ildac["MAG_APER"],ildac["MAGERR_APER"]
        mag_psf,magerr_psf=ildac["MAG_PSF"],ildac["MAGERR_PSF"]
        index,ira,idec=getindex_loc(tra,tdec,ra, dec)
        ams=ams_cal(ldacList[i],ira,idec)
        psfmag=mag_psf[index]#+ztmag
        #logger.info('index',index)
        logger.info('psfmag',mag_psf[index])
        psfarr=magerr_psf[index]
        obstime,utctime,tdbtime=TransLDAC_time(ldacList[i])
        data_array.append([utctime,tdbtime,ira,idec,ams,psfmag,psfarr])
        #if(psfmag>-13):
        
#  continue
        amsx.append(ams)
        magy.append(psfmag)
        magyerr.append(psfarr)
    plt.switch_backend('agg')
    plt.errorbar(amsx, magy, yerr=magyerr,fmt='o',ecolor='r',color='r',elinewidth=1,capsize=2,ms=3)
    #plt.plot(mjdx, magy_psf)
    plt.xlabel('airmass')
    plt.ylabel('MI(psf_mag)')
    plt.title(ttfname)
    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig(rootpath+date+'_'+ttfname+'_airmass_mag.pdf')
    logger.info('ok')


def pro(ttfname,date):
    #tra=322.7362855
    #tdec=44.346236

    #tra=69.0428297
    #tdec=19.87394052
    rootpath=get_current_dir()
    #ttfname_b,ttfname_v='y50b_ZTF2130+4420_jb','y50a_ZTF2130+4420_jv'
    #phase_plot(rootpath,ttfname_b,ttfname_v,date,tra,tdec)
    #plot(rootpath,ttfname,date,ccdindex,tra,tdec)
    #plot_combine(rootpath,ttfname,date,tra,tdec,ccdindex)
    #plot(rootpath,ttfname,date,tra,tdec)
    #airms_plot(rootpath,ttfname,date)
    #mag_plot(rootpath,ttfname,date)
    ttfname1='y50a_J043610.27+195226.1_mg'
    ttfname2='y50b_J043610.27+195226.1_mr'
    #mag_plot_2(rootpath,ttfname1,ttfname2,date,tra,tdec)
    #hzg1 g=14.5719 r=14.9510 g191 g=11.6792 r=12.0664 g193 g=15.6357 r=15.5846
    #mmag=14.9510
    mmag=0
    #hz21#tra,tdec=183.483333, 32.941667
    #g191#tra,tdec=76.375000,52.831667i
    #g193#tra,tdec=118.362500,52.493056

    tra,tdec=118.3636,52.4920
    #tra,tdec=183.48443,32.942044
    #tra,tdec=76.3775,52.831089
    ams_plot(rootpath,ttfname,date,tra,tdec)
    mag_plot(rootpath,ttfname,date,tra,tdec)
    logger.info('the astrometric process is done')



# ttfname = sys.argv[1]
# date = sys.argv[2]
# pro(ttfname,date)
