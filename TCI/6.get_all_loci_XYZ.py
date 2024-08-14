import pandas as pd
import numpy as np
import os
import math

def mkdir(path):
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False


def savetxt(filename, x):
    np.savetxt(filename, x, delimiter='\t', fmt='%s')


def sort_list(list_in):
    list_out = sorted(list_in, key=lambda items: int(items.split('\t')[1]))
    list_out = sorted(list_out, key=lambda items: items.split('\t')[0])
    return list_out
mkdir('2.XYZ')
path = 'K:/2024_NC/HiC_normalize/KR_MCFS/NHEK/7.FLAMINGO'
dirs = os.listdir(path)
all_list = []
data = pd.DataFrame()
for i in dirs:
    if i == 'chrY':
        continue
    print(i)
    f = path + '/' + i + '/'+i
    res_list = []
    df = pd.read_csv(f,sep ='\t')
    df['chr'] = [i for j in range(len(df))]
    df['start'] = (df['frag_id']-1)*10000+1
    df['end'] = (df['frag_id'])*10000
    df = df.loc[:,['chr','start','end',df.columns[0],df.columns[1],df.columns[2],df.columns[3]]]
    df.to_csv('2.XYZ/'+i,sep = '\t',header = True,index = False)
    data = data.append(df)
data.to_csv('2.XYZ/'+'all_coordinate',sep = '\t',header = True,index =False)



