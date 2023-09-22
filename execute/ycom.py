import argparse
import time
from config import config

from lib.LogInstance import Logger
import traceback
import matplotlib
import glob
from tqdm import tqdm
matplotlib.use('Agg')
import os
import yprecom1 as ypre
import yphostcom1 as ypho
import yrecpcom1 as yr

def compro(date):
    yr.pro(date)
    filepath=config["data"]["rootpath"]+'/reception/'+date+'/raw/fileguide/'
    
    scilist = glob.glob(filepath + 'sci*')
    loggerloguru.info(f"scilist:{len(scilist)}")
    for item in tqdm(scilist):
        loggerloguru.info("[_ycom] "+item + ' is being processed#############')
        filepath, tempfilename = os.path.split(item)
        sci, tid, objectname, filterid = tempfilename.split('_')
        filterid, lastn = filterid.split('.')
        loggerloguru.info(objectname)
        tid_target_filter = tid + '_' + objectname + '_' + filterid

        try:
            ypre.pro_combine(tid_target_filter, date)
            
            ypho.pro(tid_target_filter, date)
            
        except Exception as e:

            loggerloguru.info(
                '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            loggerloguru.error(traceback.format_exc())
            
            loggerloguru.info(
                '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            continue
if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20221207', help='输入处理日期')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_ycom").get_logger
    
    starttime=time.time()
    compro(args.date)
    endtime=time.time()
    loggerloguru.info("[_ycom] "+"time spending:{}s".format(endtime-starttime))


