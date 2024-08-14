import pandas as pd
import numpy as np
import os

df = pd.read_csv('chrom_hg19.sizes',sep = '\t',header = None)

ff = '/mnt/w/DSB_TAD_density/0.NHEK_hic2matrix/1-hic/GSE63525_NHEK_combined_30.hic'
di = 'NHEK' ##cell_type
for i in range(len(df)):
    chr = df.iloc[i, 0]
    leng = df.iloc[i, 1]
    os.system(f'mkdir -p 1.KR_matrix/{di.split("_")[0]}/{chr}')
    os.system(
        f'java -jar ./juicer_tools_1.22.01.jar dump norm KR {ff} {chr}:0:{leng} BP 10000 ./1.KR_matrix/{di.split("_")[0]}/{chr}/{chr}_10kb.KRnorm')
    os.system(
        f'java -jar ./juicer_tools_1.22.01.jar dump observed KR {ff} {chr}:0:{leng} {chr} BP 10000 ./1.KR_matrix/{di.split("_")[0]}/{chr}/{chr}_10kb.KRnormMatrix')

    os.system(
        f'java -jar ./juicer_tools_1.22.01.jar dump observed NONE {ff} {chr}:0:{leng} {chr} BP 10000 ./1.KR_matrix/{di.split("_")[0]}/{chr}/{chr}_10kb.RAWobserved')

    os.system(
        f'java -jar ./juicer_tools_1.22.01.jar dump oe KR {ff} {chr}:0:{leng} {chr} BP 10000 ./1.KR_matrix/{di.split("_")[0]}/{chr}/{chr}_10kb.KRoeMatrix')

##################
"""
Process NaN
"""
path = '1.KR_matrix'
dirs = os.listdir(path)
for di in dirs:
    for dii in os.listdir(path + '/' +di):
        f_kr = path + '/' + di + '/' +dii+'/'+ dii+'_10kb.KRnormMatrix'
        f_oe = path + '/' + di + '/' +dii+'/'+ dii+'_10kb.KRoeMatrix'

        df_kr = pd.read_csv(f_kr,sep = '\t',header = None)
        df_oe = pd.read_csv(f_oe,sep = '\t',header = None)
        df_kr.dropna(inplace = True)
        df_kr.to_csv(f"{path + '/' + di + '/' +dii+'/'+dii}_10kb_process.KRnormMatrix",sep = '\t',header = None,index = False)

        df_oe.dropna(inplace=True)
        df_oe.to_csv(f"{path + '/' + di + '/' +dii+'/'+dii}_10kb_process.KRoeMatrix",sep = '\t',header = None,index = False)

