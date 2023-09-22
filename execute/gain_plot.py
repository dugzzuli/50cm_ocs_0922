import os
import sys
from loguru import logger

import matplotlib.pyplot as plt
import numpy as np
import datetime
from PIL import Image

from loguru import logger
from lib.LogInstance import Logger
from lib.phot.ybias import get_current_dir



def load_gain(recption_path,day,during_days=-15):
    
    now=datetime.datetime.strptime(day,'%Y%m%d').date()
     
    #datelist=[]
    gain_50a_list=[]
    gain_50b_list=[]
    filelist_50a=['y50a_mg_bin1','y50b_mg_bin1']#,'y50a_mr_bin2','y50a_mi_bin2']
    #filelist_50b=['y50b_mg_bin2']#,'y50b_mr_bin2','y50b_mi_bin2']
    
    for item in filelist_50a:
        gn_name=item+'_gn.npy'
        try:
            gn=np.load(recption_path+str(day)+'/cal/'+gn_name,allow_pickle=True)
             
            gain_l=gn[0][0]
            gain_l_std=gn[0][2]
            gain_r=gn[1][0]
            gain_r_std=gn[1][2]
            if(max(gain_l_std,gain_r_std)>0.005):
                filelist_50a=['y50a_mr_bin2','y50b_mr_bin2']
                for item in filelist_50a:
                    gn_name=item+'_gn.npy'
                     
                    gn=np.load(recption_path+str(day)+'/cal/'+gn_name,allow_pickle=True)
             
                    gain_l=gn[0][0]
                    gain_l_std=gn[0][2]
                    gain_r=gn[1][0]
                    gain_r_std=gn[1][2]

        except:
            logger.info('there is no gain calculation this day!')
            continue
        gr_list=[]
        gl_list=[]
        datelist=[]
 
        for i in range(0,abs(during_days)):

            date = now + datetime.timedelta(days = during_days+i+1)
            #gn=[[a_gain,a_noise,a_gain_std,a_noise_std],[b_gain,b_noise,b_gain_std,b_noise_std],a_gainlist,a_noiselist,b_gainlist,b_noiselist]
             
            date1=str(date).replace('-','')
             
            logger.info(recption_path+str(date1)+'/cal/'+gn_name)
            
            if (os.path.exists(recption_path+str(date1)+'/cal/'+gn_name)):

                datelist.append(date)
                gn=np.load(recption_path+str(date1)+'/cal/'+gn_name,allow_pickle=True)
                #logger.info(date1+item+'_gain')
                #logger.info(gn[0])
                gl_list.append(np.array(gn[0]))
                gr_list.append(np.array(gn[1]))
        gl_list=np.array(gl_list)
        gr_list=np.array(gr_list)
        #datelist = [datetime.strptime(d, '%Y/%m/%d').date() for d in datelist]
         
        logger.info(datelist)

        plt.switch_backend('agg')
        plt.figure(dpi=300,figsize=(6,6))
        u1=np.max([gl_list[:,0].max(),gr_list[:,0].max()])
        u2=np.min([gl_list[:,0].min(),gr_list[:,0].min()])
        umean=(u1+u2)/2
        plt.ylim((u1+u2)/2-0.075,(u1+u2)/2+0.075)
        plt.errorbar(datelist, gl_list[:,0], yerr=gl_list[:,2],fmt='o',ecolor='r',color='r',elinewidth=2,capsize=4,label='gate_left')
        plt.errorbar(datelist, gr_list[:,0], yerr=gr_list[:,2],fmt='o',ecolor='b',color='b',elinewidth=2,capsize=4,label='gate_right')
        plt.legend()
        plt.xlabel('Date',fontsize=10)
        plt.legend()
        plt.ylabel(item+'_gain(e-/ADU)')

        plt.title(day+'_gain='+str(round(gain_l,3))+'+/-'+str(round(gain_l_std,3))+','+str(round(gain_r,3))+'+/-'+str(round(gain_r_std,3)))
        plt.xticks(rotation=15,fontsize=8)
        plt.savefig(recption_path+str(day)+'/cal/figure/'+day+item+'_gain.pdf')
        plt.clf() 
        plt.close()
    logger.info('the gain checking process is done')


def join(png1, png2, flag='horizontal'):
    """
    :param png1: path
    :param png2: path
    :param flag: horizontal or vertical
    :return:
    """
    img1, img2 = Image.open(png1), Image.open(png2)
    # 统一图片尺寸，可以自定义设置（宽，高）
    img1 = img1.resize((1500, 1000), Image.ANTIALIAS)
    img2 = img2.resize((1500, 1000), Image.ANTIALIAS)
    size1, size2 = img1.size, img2.size
    if flag == 'horizontal':
        joint = Image.new('RGB', (size1[0] + size2[0], size1[1]))
        loc1, loc2 = (0, 0), (size1[0], 0)
        joint.paste(img1, loc1)
        joint.paste(img2, loc2)
        joint.save('horizontal.png')
    elif flag == 'vertical':
        joint = Image.new('RGB', (size1[0], size1[1] + size2[1]))
        loc1, loc2 = (0, 0), (0, size1[1])
        joint.paste(img1, loc1)
        joint.paste(img2, loc2)
        joint.save('vertical.png')
 

 
def pro(date):
    
    rootpath=get_current_dir()
    recption_path=rootpath+'reception/'
    load_gain(recption_path,date)

   


               

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()  # 创建parser
    parser.add_argument('--date', type=str, default='20221207', help='输入处理日期')
    
    args = parser.parse_args()  # 参数解析
    loggerloguru = Logger(args.date, '', "_yphost").get_logger
    pro(args.date)
    # date= sys.argv[1]
    # pro(date)
 
