import glob
import traceback
import os
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
 
# from execute import yprecom1 as ypre
# from execute import yphostcom1 as ypho
# from execute import yrecpcom1 as yr
from execute import yprecom0_new as ypre
from execute import yphostcom0 as ypho
from execute import yrecpcom1_new as yr
 
 
import argparse
import gc
from multiprocessing import JoinableQueue, Process
import time
 
import traceback
 
import glob
 
matplotlib.use('Agg')

import ycom_ttf
 
 

def mtarget(datelist):
 
     
    #sci_root=rootdir+'reception/'+str(date)+'/sci/'
    for i in range(len(datelist)):
        date=datelist[i]
        
        os.system("python yrecp.py --date=" + date )
         
        try:
                
                os.system("python ycom_ttf.py --date=" + date)
        except:
                continue
                

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()  # 创建parser
    # parser.add_argument('--date', type=str, default="20230517", help='输入处理日期')
    # #parser.add_argument('--filename', type=str, default="/home/mpilot/1m6/rawdata/20221107/mb_testdata_20221107/mb_sc_h04hz22_u_20221107185350_006.fits", help='输入处理日期')
    # parser.add_argument('--filename', type=str, default="/home/mpilot/1m6/rawdata/20230517/my_testdata_20230517/my_sc_tp00929_g_20230517185731_215.fits", help='输入处理日期')
   #ngc5457_date=['20230518','20230510','20230429','20230420','20230418','20230417','20230406','20230331','20230314','20230312','20230310','20230307']
   date_12377=['20221129','20221128','20221127','20221126','20221125','20221124']
    #ngc5457_date=['20230518','20230510','20230508','20230429','20230420','20230418','20230417']

   #mtarget(ngc5457_date)
   mtarget(date_12377)
    #args = parser.parse_args()  # 参数解析
    #loggerloguru = Logger(args.date, '', "_msingle").get_logger 
 
    #tx=forTimeGetAlt() #获取开始时间
     
 
