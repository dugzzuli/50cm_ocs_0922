#2023 0128 创立




# ycom.py 原始串行



# git提交
## 1.  git add .
## 2.  git commit -m "每天要修改的内容"
## 3.  git push origin main

# 配置文件路径

# datacenter
data:
  rawdata: '/data2/50cm/raw/'
  rootpath: '/data2/workspace/50cm_pro/'
refcat:
  gaia: '/data3/Other_surveys/GaiaEDR3FITS/GAIAEDR3New/'

# m1
data:
  rawdata: '/home/50cm/data2/raw/'
  rootpath: '/home/50cm/50cm_scripts/50cm_pro'
refcat:
  gaia: '/home/50cm/data2/GAIAEDR3New/'

# 注意
config.yaml 修改之后需要cp config.yaml ../config.yaml