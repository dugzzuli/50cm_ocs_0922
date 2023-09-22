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
# from execute import yprecom1 as ypre
# from execute import yphostcom1 as ypho
# from execute import yrecpcom1 as yr
from execute import yprecom0_new as ypre
from execute import yphost0 as ypho
from execute import yrecpcom1_new as yr

from lib.gchelp import __clear_env

def compro_pro(date,q):
    #yr.pro(date)
    filepath=config["data"]["rootpath"]+'/reception/'+date+'/raw/fileguide/'
    #scilist = glob.glob(filepath + 'sci_'+'*NGC5457*')
    scilist = glob.glob(filepath + 'sci_'+'*12377*')
    print(scilist)
    
    loggerloguru.info(f"scilist:{len(scilist)}")
    for item in tqdm(scilist):
        
        q.put(item)
    q.join()
        

def compro_consumer(date,q):    
    time.sleep(1)
    while True:
        item=q.get()
        loggerloguru.info("[_ycom] "+item + ' is being processed#############')
        filepath, tempfilename = os.path.split(item)
        _, tid, objectname, filterid = tempfilename.split('_')
        filterid, lastn = filterid.split('.')
        loggerloguru.info(objectname)
        tid_target_filter = tid + '_' + objectname + '_' + filterid
        ypho.pro(tid_target_filter, date)
        try:
            # __clear_env()
            
            ypre.pro_combine(tid_target_filter, date)
            
            ypho.pro(tid_target_filter, date)
            
            import yapercor as apr
            
            apr.pro(tid_target_filter, date)
            
            
                
        except Exception as e:

            loggerloguru.info(
                '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            loggerloguru.error(traceback.format_exc())
            
            loggerloguru.info(
                '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        finally:
            q.task_done()
            
            for k in list(locals().keys()):
                # if locals[k] is np.nan:
                try:
                    del locals[k]
                except:
                    continue
            gc.collect()
            
            
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20221102', help='输入处理日期')
    parser.add_argument('--maxsize', type=int, default=2, help='最大队列')
    args = parser.parse_args()  # 参数解析
    print(args.date)
    loggerloguru = Logger(args.date, '', "_ycom").get_logger
    
    q = JoinableQueue()
    starttime=time.time()
    
    maxsize=args.maxsize
    q = JoinableQueue(maxsize=maxsize)
   
    # 创建 多进程 生产
    p1 = Process(target=compro_pro, args=(args.date,q))
    

    # 创建多进程 消费
    task_p_list=[]
    for task_con in range(maxsize):
        temp_task=Process(target=compro_consumer, args=(args.date, q))
        temp_task.daemon = True
        task_p_list.append(temp_task)
        

    # 全部启动
    p_l = [p1]
    for p in p_l:
        p.start()
        
    for c in (task_p_list):
        c.start()
        
    # 等待生产进程结束
    p1.join()
    
    print("生产结束")
    # 因为生产者内部,一直等待队列空才会退出,所以生产结束意味着全部结束了

    endtime=time.time()
    loggerloguru.info("[_ycom] "+"time spending:{}s".format(endtime-starttime))


