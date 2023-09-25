import traceback
import os
import sys
import glob
import numpy as np
from numpy.lib.polynomial import _roots_dispatcher
import pandas as pd
from astropy.io import fits
from config import config
from execute.flatdiff import mkdir
from loguru import logger
from lib.LogInstance import Logger
from lib.phot.ycalcom1 import savenpy, splitfilename
from lib.phot.ybias import get_current_dir

# from numba import jit
# @jit #(nopython=True)

def datamove(spath,tpath):#here spath must be a filepath not a dir
    # adding exception handling
    #if(1==1):
    if(os.path.exists(tpath)):
        logger.info('the tpath is: {}'.format(tpath))
        logger.info('the files you want to move may exists')
        sfile=len(glob.glob(spath+'*.fit*'))
        tfile=len(glob.glob(tpath+'*.fit*'))
        if(tfile-sfile>0 or tfile==0):
            logger.info('the files you want to move may not complete in your target directory,we will move them again')
            
            for fit_file in glob.glob(spath+'*.fit*'):
            #logger.info(fit_file)
                try:
                #shutil.copy(fit_file, tpath)
                    os.system("ln -sf %s %s"%(fit_file, tpath)) 
                except:
                    logger.info("Unexpected error:", sys.exc_info())
            for wcs_file in glob.glob(spath+'*.wcs*'):
            #logger.info(fit_file)
                try:
                #shutil.copy(fit_file, tpath)
                    os.system("ln -sf %s %s"%(wcs_file, tpath)) 
                except:
                    logger.info("Unexpected error:", sys.exc_info())
        logger.info('the files are moved successfully')

        
    else:
        mkdir(tpath)
        for fit_file in glob.glob(spath+'*.fit*'):
            #logger.info(fit_file)
            try:
                #shutil.copy(fit_file, tpath)
                os.system("ln -sf %s %s"%(fit_file, tpath)) 
            except:
                logger.info("Unexpected error:{}".format( sys.exc_info()))
        
        logger.info('the files are moved successfully')



def single_datamove(spath,tpath,filename):
    # adding exception handling
    #if(1==1):
    if(os.path.exists(tpath) and os.path.exists(spath+filename)):
        fit_file = spath+filename
        logger.info('the tpath is: {}'.format(tpath))
        logger.info('the files you want to move may exists')
          
        try:
            os.system("ln -sf %s %s"%(fit_file, tpath)) 
        except:
            logger.info("Unexpected error:", sys.exc_info())
        logger.info('the files are moved successfully')

        
    else:
        mkdir(tpath)
         
        try:
            #shutil.copy(fit_file, tpath)
            os.system("ln -sf %s %s"%(fit_file, tpath)) 
        except:
            logger.info("Unexpected error:{}".format( sys.exc_info()))
        
        logger.info('the files are moved successfully')



def ttfname_header(filename):
    hdr=fits.getheader(filename)
    #OBJECT  = 'y50b_HZ44_' 
    tt=hdr['OBJECT']
    
    tid,target = tt.split('_')[0],tt.split('_')[1]
    filterid='m'+hdr['FILTER']
    
    logger.info('++++++++++++++++++++++++++++++++++++++++')

    logger.info(tid+'+'+target+'+'+filterid)
    logger.info('++++++++++++++++++++++++++++++++++++++++')
    return tid,target,filterid
 


def dataclassify(rawdir,fileguide_path):
    fblist = glob.glob(rawdir+'*bias*.fit*')
    if(fblist==False):
       fblist = glob.glob(rawdir+'*BIAS*.fit*')
    fflist = glob.glob(rawdir+'*flat*.fit*')
    if(fflist==False):
        fflist = glob.glob(rawdir+'*FLAT*.fit*')
    
    fdlist = glob.glob(rawdir+'*dark*.fit*')
    if(fdlist==False):
        fflist = glob.glob(rawdir+'*DARK*.fit*')
    flist= glob.glob(rawdir+'*.fit*')
    scilist=list(set(flist).difference(set(fblist)).difference(set(fflist)).difference(set(fdlist)))
    #sci images classification by target name 
    s_objset=set()
    for j in range(0,len(scilist)):
        filepath,tempfilename = os.path.split(scilist[j])
        try:
            tid,objectname,filterid,file_id=splitfilename(tempfilename)
            filterid='m'+filterid[-1]
            hdr=fits.getheader(scilist[j])
        except:
            logger.error(tempfilename+ ' can not be splited\n'+str(traceback.format_exc()))
            
            try:
                tid,objectname,filterid=ttfname_header(filepath)
            except Exception as e:
                logger.error(str(traceback.format_exc()))
                logger.error('the header information can not extract tid, target and filter')
                continue

        s_objset.add((tid,filterid,objectname,hdr['XBINNING']))
    # for var in np.sort(scilist):
    #     logger.info(var,end='\n')
    objlist=list(s_objset)
    #将文件分类并存放
    mkdir(fileguide_path)
    for i in range(0,len(objlist)):
        tid,filterid,target,bins=objlist[i]
        filterid=filterid
        taglist = glob.glob(rawdir+tid+'_'+target+'*'+filterid+'*.fit*')
        # logger.info('###########################################################')
        # for var in np.sort(taglist):
        #     logger.info(var,end='\n')
        ttfname=tid+'_'+target+'_'+filterid
        savenpy(fileguide_path,'sci_'+ttfname+'.npy',taglist)
        logger.info(ttfname+'  count: '+str(len(taglist))+ '  bins='+str(bins) )
    return objlist



def darkclassify(rawdir,fileguide_path):
    fdlist = glob.glob(rawdir+'*dark*.fit*')
    dlistbin = []
    for j in range(0,len(fdlist)):
        filepath,tempfilename = os.path.split(fdlist[j])
        try:
           tid,objectname,file_id=tempfilename.split('_')
           hdr=fits.getheader(fdlist[j])
           dlistbin.append(list([fdlist[j],str(tid), hdr['XBINNING']]))
        except:
           logger.info('the dark file: ' +tempfilename+' is can not be splited, please check' )
           continue
    # for var in np.sort(scilist):
    dlistbin=np.array(dlistbin)
    try:
       darkdf=pd.DataFrame({'filename':list(dlistbin[:,0]),'tid':list(dlistbin[:,1]),'bins':list(dlistbin[:,2])})
       darkgroup = darkdf.groupby(['tid','bins'])
       for key,values in darkgroup:
           tdname=str(key[0])+'_dark_bin'+str(key[1])
           logger.info(tdname+'  count: {}'.format(len(list(values['filename']))))
           savenpy(fileguide_path,tdname+'.npy',list(values['filename']))
    except:
         logger.info('the dark  file is not exists or filename is not qualify,please check!')
    return dlistbin



def biasclassify(rawdir,fileguide_path):
    fblist = glob.glob(rawdir+'*bias*.fit*')
    fflist = glob.glob(rawdir+'*flat*.fit*')
    fdlist = glob.glob(rawdir+'*dark.fit*')
    blistbin = []
    for j in range(0,len(fblist)):
        filepath,tempfilename = os.path.split(fblist[j])
        try:
           tid,objectname,file_id=tempfilename.split('_')
           hdr=fits.getheader(fblist[j])
           blistbin.append(list([fblist[j],str(tid), hdr['XBINNING']]))
        except:
           logger.info('the bias file: ' +tempfilename+' is can not be splited, please check' )
           continue
    # for var in np.sort(scilist):
    blistbin=np.array(blistbin)
    try:
       biasdf=pd.DataFrame({'filename':list(blistbin[:,0]),'tid':list(blistbin[:,1]),'bins':list(blistbin[:,2])})
       biasgroup = biasdf.groupby(['tid','bins'])
       for key,values in biasgroup:
           tbname=str(key[0])+'_bias_bin'+str(key[1])
           logger.info(tbname+'  count: {}'.format(len(list(values['filename']))))
           savenpy(fileguide_path,tbname+'.npy',list(values['filename']))
    except:
         logger.info('the bias  file is not exists or filename is not qualify,please check!')
    return blistbin



def flatclassify(rawdir,fileguide_path):
    fflist = glob.glob(rawdir+'*flat*.fit')
    flistbin = []
    for j in range(0,len(fflist)):
        filepath,tempfilename = os.path.split(fflist[j])
        try:
            tid,objectname,filterid,file_id=tempfilename.split('_')
            hdr=fits.getheader(fflist[j])
            flistbin.append(list([fflist[j],str(tid),str(filterid), hdr['XBINNING']]))
        except:
            logger.info('the flat file: ' +tempfilename+' is can not be splited, please check' )
            continue
    # for var in np.sort(scilist):
    flistbin=np.array(flistbin)
    try:
        flatdf=pd.DataFrame({'filename':list(flistbin[:,0]),'tid':list(flistbin[:,1]),'filterid':list(flistbin[:,2]),'bins':list(flistbin[:,3])})
        flatgroup = flatdf.groupby(['tid','filterid','bins'])
        for key,values in flatgroup:
        #logger.info("the tid and bins：",key,)
        #logger.info(values,'\n')
            tbname=str(key[0])+'_flat_'+str(key[1])+'_bin'+str(key[2])
            logger.info(tbname+'  count: {}'.format(len(list(values['filename']))))
            savenpy(fileguide_path,tbname+'.npy',list(values['filename']))
    except:
         logger.info('the flat file is not exists or filename is not qualify,please check!')

    return flistbin     



def pro(date):
    
    loggerloguru = Logger(date, '', "_ycep").get_logger
    #filedate = sys.argv[2]
    rootdir=get_current_dir()
    loggerloguru.info("[_yrcep] "+'#######################'+str(rootdir))
    rawpath=rootdir+'reception/'+str(date)+'/raw/'
    fileguide_path=rawpath+'fileguide/'
    #sourcedir='/data2/50cm/raw/'+str(date)+'/'
    sourcedir=config["data"]["rawdata"]+str(date)+'/'
    
    if(os.path.exists(sourcedir)):
       datamove(sourcedir,rawpath)
       countall = len(glob.glob(rawpath+'*.fit*'))
        
       dataclassify(rawpath,fileguide_path)
       flatclassify(rawpath,fileguide_path)
       biasclassify(rawpath,fileguide_path)
       darkclassify(rawpath,fileguide_path)
       loggerloguru.info("[_yrcep] "+'the data '+ str(countall)+ ' is collected')
    else:
       loggerloguru.info("[_yrcep] "+'there is no raw files this day!please check!!')



if __name__ == "__main__":
    import argparse,time
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20221207', help='输入处理日期')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_ycom").get_logger
    starttime=time.time()
    pro(args.date)
    endtime=time.time()
    loggerloguru.info("[_yprecom1] "+"time spending:{}s".format(endtime-starttime))