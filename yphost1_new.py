import numpy as np
import os, sys

from loguru import logger
from lib.LogInstance import Logger
from lib.phot.ybias import get_current_dir
from execute import yphostcom0 as ypho
from loguru import logger

def pro(tid_target_filter,date):
    rootpath=get_current_dir()
    ypho.photometry(rootpath,tid_target_filter,date)

    logger.info('the astrometric process is done')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20230313', help='输入处理日期')
    parser.add_argument('--tid_target_filter', type=str,default="y50b_NGC3799_mg", help='输入处理日期')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_yphost").get_logger
    pro(args.tid_target_filter,args.date)
    
    
