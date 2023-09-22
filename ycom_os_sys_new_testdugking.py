import argparse
import gc
from multiprocessing import JoinableQueue, Process
import time
from config import config

from lib.LogInstance import Logger
import traceback
import matplotlib
import glob
from tqdm import tqdm
matplotlib.use('Agg')
import os
from execute import yprecom1_new as ypre
from execute import yphostcom1_new as ypho
from execute import yrecpcom1_new as yr
from lib.gchelp import __clear_env, elapsed_time
#仅仅用于测试


@elapsed_time
def compro_consumer(date):
    
    yr.pro(date)
    filepath=config["data"]["rootpath"]+'/reception/'+date+'/raw/fileguide/'
    scilist = glob.glob(filepath + 'sci*')
    loggerloguru.info(f"scilist:{len(scilist)}")
    for item in tqdm(scilist):
        loggerloguru.info("[_ycom] "+item + ' is being processed#############')
        filepath, tempfilename = os.path.split(item)
        _, tid, objectname, filterid = tempfilename.split('_')
        filterid, lastn = filterid.split('.')
        loggerloguru.info(objectname)
        tid_target_filter = tid + '_' + objectname + '_' + filterid

        try:
            # __clear_env()
            
            ypre.pro_combine(tid_target_filter, date)
            
            ypho.pro(tid_target_filter, date)
            import yapercor as apr
            apr.pro(tid_target_filter, date)
            pass
                
        except Exception as e:

            loggerloguru.info(
                '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            loggerloguru.error(traceback.format_exc())
            
            loggerloguru.info(
                '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        finally:
            
            
            for k in list(locals().keys()):
                # if locals[k] is np.nan:
                try:
                    del locals[k]
                except:
                    continue
            gc.collect()
            
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20230207', help='输入处理日期')
    parser.add_argument('--maxsize', type=int, default=30, help='最大队列')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_ycom").get_logger

    starttime=time.time()
    compro_consumer(args.date)
    endtime=time.time()
    loggerloguru.info("[_ycom] "+"time spending:{}s".format(endtime-starttime))


