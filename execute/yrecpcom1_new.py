import traceback
import os
import sys
import glob
import numpy as np
import pandas as pd
from astropy.io import fits
from config import config
from execute.flatdiff import mkdir
from loguru import logger
from lib.LogInstance import Logger
from lib.phot.ycalcom1 import savenpy, splitfilename
from lib.phot.ybias import get_current_dir


def datamove(spath,tpath):#here spath must be a filepath not a dir

    if(os.path.exists(tpath)):
        logger.info('the tpath is: {}'.format(tpath))
        logger.info('the files you want to move may exists')
        if(1==1):
            for fit_file in glob.glob(spath+'*.fit*'):
                # print(fit_file)
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
    

def ttfname_header(filename):
    hdr=fits.getheader(filename)
    #OBJECT  = 'y50b_HZ44_' 
    tt=hdr['OBJECT']
    try:
        tid,target = tt.split('_')[0],tt.split('_')[1]
        filterid='m'+hdr['FILTER'].lower()

        # logger.info('++++++++++++++++++++++++++++++++++++++++')

        # logger.info(tid+'+'+target+'+'+filterid)
        # logger.info('++++++++++++++++++++++++++++++++++++++++')
        return tid,target,filterid
    except:
        return None
    

def ttfname_header_new(filename):
    hdr=fits.getheader(filename)
    #OBJECT  = 'y50b_HZ44_' 
    target=hdr['OBJECT']#.lower()
    tid = hdr['TELEID'].lower()
     
    filterid='m'+hdr['FILTER'].lower()

    # logger.info('++++++++++++++++++++++++++++++++++++++++')

    # logger.info(tid+'+'+target+'+'+filterid)
    # logger.info('++++++++++++++++++++++++++++++++++++++++')
    return tid,target,filterid

# from numba import jit
def dataclassify_new1(rawdir,fileguide_path):
    fblist = glob.glob(rawdir+'*BIAS*.fit*')
    if(len(fblist)==0):
       fblist = glob.glob(rawdir+'*bias*.fit*')
    fflist = glob.glob(rawdir+'Flat*.fit*')
    if(len(fflist)==0):
        fflist = glob.glob(rawdir+'*flat*.fit*')
    
    fdlist = glob.glob(rawdir+'*DARK*.fit*')
    othertest = glob.glob(rawdir+'*test*.fit*')
    otherfocus = glob.glob(rawdir+'*focus*.fit*')
    if(len(fdlist)==0):
        fdlist = glob.glob(rawdir+'*dark*.fit*')
    flist= glob.glob(rawdir+'*.fit*')
    scilist=list(set(flist).difference(set(fblist)).difference(set(fflist)).difference(set(fdlist)).difference(set(othertest)).difference(set(otherfocus)))
    # logger.info('^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    # logger.info(scilist)
    # logger.info('^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    #sci images classification by target name 
    s_objset=set()
    for j in range(0,len(scilist)):
        filepath,tempfilename = os.path.split(scilist[j])
        try:
            #logger.info(scilist[j])
            tid,objectname,filterid=ttfname_header_new(scilist[j])
            hdr=fits.getheader(scilist[j])
        except Exception as e:
            #logger.error(str(traceback.format_exc()))
            logger.error('the header information can not extract tid, target and filter')
            try:
                logger.info(scilist[j])
                filepath,tempfilename = os.path.split(scilist[j])
                sname=tempfilename.split('_')
                #tid,objectname,filterid,file_id=splitfilename(tempfilename)
                tid =sname[0]
                objectname = sname[1]
                filterid = sname[2]
                hdr=fits.getheader(scilist[j])
            except Exception as e:
                #logger.error(str(traceback.format_exc()))
                logger.error('the filename can not splited ')
                continue

        s_objset.add((tid,filterid,objectname,hdr['XBINNING']))
    # for var in np.sort(scilist):
    #     logger.info(var,end='\n')
    objlist=list(s_objset)
    #将文件分类并存放
    mkdir(fileguide_path)
    for i in range(0,len(objlist)):
        tid,filterid,target,bins=objlist[i]
        if(target=='FLAT'):
            continue
        filterid=filterid
        #taglist = glob.glob(rawdir+tid+'_'+target+'*'+filterid+'*.fit*')
        up_tid=tid.upper()
        taglist = glob.glob(rawdir+tid+'_'+target+'_'+filterid+'*.fit*')
        
        # taglist1= glob.glob(rawdir+up_tid+'_'+target+'_'+filterid+'*.fit*')
        # taglist.extend(taglist1)
        # logger.info('####################################################')
        # logger.info(rawdir+tid+'_'+target+'*'+filterid)
        # logger.info('####################################################')
        # for var in np.sort(taglist):
        #     logger.info(var,end='\n')
        ttfname=tid+'_'+target+'_'+filterid
        
        if(len(taglist)>0):
            savenpy(fileguide_path,'sci_'+ttfname+'.npy',taglist)
            logger.info(ttfname+'  count: '+str(len(taglist))+ '  bins='+str(bins) )
        # logger.info('###########################################################')
        # for var in np.sort(taglist):
        #     logger.info(var,end='\n')
         
    return objlist

def dataclassify_new(rawdir,fileguide_path):
    #IMAGETYP= 'Bias Field' /        Type of image
    #OBJECT  = 'BIAS    ',OBJECT  = 'FLAT    '
    #TELEID  = 'Y50A    '
    fblist = glob.glob(rawdir+'*BIAS*.fit*')
    fb1= glob.glob(rawdir+'*bias*.fit*')
    fblist.extend(fb1)
    #logger.info(fblist)
    fflist = glob.glob(rawdir+'*FLAT*.fit*')
    if(len(fflist)==0):
        fflist = glob.glob(rawdir+'*Flat*.fit*')
    #logger.info(fflist)
    fdlist = glob.glob(rawdir+'*DARK*.fit*')
    #logger.info(fdlist)
    flist= glob.glob(rawdir+'*.fit*')
    scilist=list(set(flist).difference(set(fblist)).difference(set(fflist)).difference(set(fdlist)))
    # logger.info('%%%%%%%%%%%%%%%%%%%%%')
    # logger.info(scilist)
    #sci images classification by target name 
    s_objset=set()
    for j in range(0,len(scilist)):
        
        #filepath,tempfilename = os.path.split(scilist[j])
        try:
            logger.info(scilist[j])
            tid,objectname,filterid=ttfname_header_new(scilist[j])
        except Exception as e:
            #logger.error(str(traceback.format_exc()))
            logger.error('the header information can not extract tid, target and filter')
            try:
                logger.info(scilist[j])
                filepath,tempfilename = os.path.split(scilist[j])
                sname=tempfilename.split()
                #tid,objectname,filterid,file_id=splitfilename(tempfilename)
                tid =sname[0]
                objectname = sname[1]
                filterid = sname[2]
                hdr=fits.getheader(scilist[j])
            except Exception as e:
                #logger.error(str(traceback.format_exc()))
                logger.error('the filename can not splited ')
                continue
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



def darkclassify_new(rawdir,fileguide_path):
    fdlist = glob.glob(rawdir+'*DARK*.fit*')
    dlistbin = []
    for j in range(0,len(fdlist)):
        filepath,tempfilename = os.path.split(fdlist[j])
        try:
          
            hdr=fits.getheader(fdlist[j])
                
            tid = hdr['TELEID'].lower()
            dlistbin.append(list([fdlist[j],str(tid), hdr['XBINNING']]))
        #    exp =int(hdr['EXPOSURE'])
        #    dlistbin.append(list([fdlist[j],str(tid), hdr['XBINNING'],str(exp)]))
        except:
            logger.info('the dark file: ' +tempfilename+' is can not be classified, please check' )
            try:
                tid,objectname,file_id=tempfilename.split('_')
            except:
                #  logger.error(str(traceback.format_exc()))
                logger.info('the dark file: ' +tempfilename+' is can not be splited, please check' )
                continue
    # for var in np.sort(scilist):
    dlistbin=np.array(dlistbin)
    # try:
    #    darkdf=pd.DataFrame({'filename':list(dlistbin[:,0]),'tid':list(dlistbin[:,1]),'bins':list(dlistbin[:,2]),'exp':list(dlistbin[:,3])})
    #    darkgroup = darkdf.groupby(['tid','bins','exp'])
    #    for key,values in darkgroup:
    #        tdname=str(key[0])+'_dark_bin'+str(key[1])+'_'+str(key[2])+'s'
    #        logger.info(tdname+'  count: {}'.format(len(list(values['filename']))))
    #        savenpy(fileguide_path,tdname+'.npy',list(values['filename']))
    # except:
    #      logger.info('the dark  file is not exists or filename is not qualify,please check!')
    # return dlistbin
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

def biasclassify_new(rawdir,fileguide_path):
    fblist = glob.glob(rawdir+'*BIAS*.fit*')
 
    blistbin = []
    for j in range(0,len(fblist)):
        filepath,tempfilename = os.path.split(fblist[j])
        try:
            #tid,objectname,file_id=tempfilename.split('_')
            hdr=fits.getheader(fblist[j])
            tid = hdr['TELEID'].lower()
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

def flatclassify_new(rawdir,fileguide_path):
    fflist = glob.glob(rawdir+'*FLAT*.fit*')
    if(len(fflist)==0):
        fflist = glob.glob(rawdir+'*Flat*.fit*')
    flistbin = []
    for j in range(0,len(fflist)):
        filepath,tempfilename = os.path.split(fflist[j])
        try:
            #tid,objectname,filterid,file_id=tempfilename.split('_')
            hdr=fits.getheader(fflist[j])
            tid = hdr['TELEID'].lower()
            filterid = 'm'+hdr['FILTER'].lower()
            flistbin.append(list([fflist[j],str(tid),str(filterid), hdr['XBINNING']]))
        except:
            logger.info('the flat file: ' +tempfilename+' can not be found in header, please check' )
            try:
                tid,objectname,filterid,file_id=tempfilename.split('_')
                flistbin.append(list([fflist[j],str(tid),str(filterid), hdr['XBINNING']]))
            except:
                logger.info('the flat file: ' +tempfilename+'  can not be splited, please check' )
                continue
            #continue    
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


# @jit #(nopython=True)
 
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



def biasclassify(rawdir,fileguide_path):
    fblist = glob.glob(rawdir+'*bias*.fit*')
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
    fflist = glob.glob(rawdir+'*flat*.fit*')
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
    
    loggerloguru = Logger(date, '', "_yrecp").get_logger
    #filedate = sys.argv[2]
    rootdir=get_current_dir()
    loggerloguru.info("[_yrcep] "+'#######################'+str(date))
    rawpath=rootdir+'reception/'+str(date)+'/raw/'
    fileguide_path=rawpath+'fileguide/'
    #sourcedir='/data2/50cm/raw/'+str(date)+'/'
    if(int(date)<20230328):
        sourcedir1=config["data"]["rawdata"]+str(date)+'/'
        sourcedir2=config["data"]["rawdata"]+str(date)+'_OCStest/50A/'
        sourcedir3=config["data"]["rawdata"]+str(date)+'_OCStest/50B/'
    else:
        sourcedir1=config["data"]["rawdata"]+str(date)+'/'
        sourcedir2=config["data"]["rawdata"]+str(date)+'_OCS/50A/'
        sourcedir3=config["data"]["rawdata"]+str(date)+'_OCS/50B/'
    if(int(date)>20230213):
        if(os.path.exists(sourcedir1)):
            datamove(sourcedir1,rawpath)
            #dataclassify_new1(rawpath,fileguide_path)
        datamove(sourcedir2,rawpath)
        datamove(sourcedir3,rawpath)
        countall = len(glob.glob(rawpath+'*.fit*'))
        dataclassify_new1(rawpath,fileguide_path)
        flatclassify_new(rawpath,fileguide_path)
        biasclassify_new(rawpath,fileguide_path)
        darkclassify_new(rawpath,fileguide_path)
        loggerloguru.info("[_yrecp] "+'the data  {}'.format(countall) + 'images  is collected')
    elif(os.path.exists(sourcedir1) and int(date)<20230214):
        datamove(sourcedir1,rawpath)
        if(os.path.exists(sourcedir2)):
            datamove(sourcedir2,rawpath)
            datamove(sourcedir3,rawpath)
        countall = len(glob.glob(rawpath+'*.fit*'))
        dataclassify(rawpath,fileguide_path)
        flatclassify(rawpath,fileguide_path)
        biasclassify(rawpath,fileguide_path)
        darkclassify(rawpath,fileguide_path)
    else:
        loggerloguru.info("[_yrecp] "+'there is no raw files this day!please check!!')
if __name__ == "__main__":
    import argparse,time
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20230225', help='输入处理日期')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_ycom").get_logger
    starttime=time.time()
    pro(args.date)
    endtime=time.time()
    loggerloguru.info("[_yrecp] "+"time spending:{}s".format(endtime-starttime))
