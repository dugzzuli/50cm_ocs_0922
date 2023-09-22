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
from execute import yrecpcom1 as yr
from execute.flatdiff import mkdir
from loguru import logger
from lib.phot.ycalcom1 import savenpy
from config import config
from loguru import logger

def lib_bak(libpath,calpath,date):
    #flat
    sfile= glob.glob(calpath+'*master_flat*')
    for i in range(0,len(sfile)):
        name= sfile[i].split('/')[-1]
        print(calpath+name)
        rename =name.split('.fits')[0]+'_'+date+'.fits'
        print(libpath+rename)
        os.system('cp %s %s' % (calpath+name,libpath+rename))
        print('ok')
        #sfile= glob.glob(calpath+'*master_flat*')
    gfile= glob.glob(calpath+'*gn*')
    for i in range(0,len(gfile)):
        name= gfile[i].split('/')[-1]
        print(calpath+name)
        rename =name.split('.npy')[0]+'_'+date+'.npy'
        print(libpath+rename)
        os.system('cp %s %s' % (calpath+name,libpath+rename))
        print('ok')
    bfile = glob.glob(calpath+'*master_bias*')
    for i in range(0,len(bfile)):
        name= bfile[i].split('/')[-1]
        print(calpath+name)
        rename =name.split('.fits')[0]+'_'+date+'.fits'
        print(libpath+rename)
        os.system('cp %s %s' % (calpath+name,libpath+rename))
        print('ok') 
    dfile = glob.glob(calpath+'*master_dark*')
    for i in range(0,len(dfile)):
        name= dfile[i].split('/')[-1]
        print(calpath+name)
        rename =name.split('.fits')[0]+'_'+date+'.fits'
        print(libpath+rename)
        os.system('cp %s %s' % (calpath+name,libpath+rename))
        print('ok') 
libpath='/home/50cm/50cm_scripts/50cm_pro/lib/cal3/'
calpath='/home/50cm/50cm_scripts/50cm_pro/reception/20230913/cal/'
#ate='20221102'
date='20230913'
lib_bak(libpath,calpath,date)
 
