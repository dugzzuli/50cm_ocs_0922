import argparse
from locale import ALT_DIGITS
import os
from re import I
import shutil
import sys
import glob
import time
import numpy as np
from numpy.lib.polynomial import _roots_dispatcher
import pandas as pd
from astropy.io import fits
from lib.LogInstance import Logger
from execute import yrecpcom1_new as yr
from execute.flatdiff import mkdir
from loguru import logger
from lib.phot.ycalcom1 import savenpy
from config import config
 



# date = sys.argv[1]
# yr.pro(date)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20230313', help='输入处理日期')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_yrecp").get_logger
    
    starttime=time.time()
    yr.pro(args.date)
    endtime=time.time()
    logger.info("[_ycom] "+"time spending:{}s".format(endtime-starttime))
