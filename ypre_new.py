import argparse
import ccdproc
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
from astropy.nddata import CCDData
import astropy.units as u
from astropy.table import Table
#import matplotlib.pyplot as plt
import time as ts
import re
import sys
import os
import glob
from lib.LogInstance import Logger 
import lib.phot.ygn as ygn 
import lib.phot.ybias as yb#creat_folder, single_fits_pack, creat_hdu_header,gen_overscan
import lib.phot.ylflat  as yf
from numpy.core.overrides import array_function_from_dispatcher 
from photutils.segmentation import make_source_mask
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground,Background2D 
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
import lib.phot.ygn as ygn 
import lib.phot.ybias as yb#creat_folder, single_fits_pack, creat_hdu_header,gen_overscan
import lib.phot.ylflat  as yf
import lib.phot.ydark as yd
import logging
from execute import yprecom0_new as ypre
from loguru import logger
import time

# tid_target_filter = sys.argv[1]
# date = sys.argv[2]
# ypre.pro_combine(tid_target_filter,date)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20230313', help='输入处理日期')
    parser.add_argument('--tid_target_filter', type=str,default="y50b_NGC3799_mg", help='输入处理日期')
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_ycom").get_logger
    
    starttime=time.time()
    ypre.pro_combine(args.tid_target_filter,args.date)
    endtime=time.time()
    logger.info("[_ypre] "+"time spending:{}s".format(endtime-starttime))


