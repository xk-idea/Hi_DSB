from gcMapExplorer import lib as gmlib
import os
import sys
import numpy as np
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from math import *
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
import pandas as pd
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
#######################################################################
resolution = 10000
path = '/mnt/k/2024_NC/HiC_normalize/KR_MCFS/NHEK/1.KR_matrix/NHEK'
path_fragment = '/mnt/w/DSB_TAD_density/call_ISTAD/hic_KR_sparse2dense/HCT116/2.1each_chr_10_fragment'

dirs = os.listdir(path)
mkdir('2.KR_MCFS_normalized_dense_matrix')
mkdir('2.KR_norm_process')
log = {}
for di in dirs:
    print(di)
    # di = chr1
    out = f'2.log_out/{di}'
    mkdir(f'2.log_out/{di}')

    KR_out = '2.KR_norm_process/'+di
    mkdir('2.KR_norm_process/'+di)

    log[di] = []
    for di_ in os.listdir(path + '/' + di):
        # di_ = chr1_10kb.KRnorm
        if str(resolution // 1000)+'kb.RAWobserved' not in di_:
            continue

        f_kr = path + '/' + di + '/' + di_
        cooReader = gmlib.importer.CooMatrixHandler(f'{f_kr}')
        cooReader.save_ccmaps(f'{out}/{di}_{resolution // 1000 }kb_raw_from_text.ccmap', xlabels=f'{di}')
        cooReader.save_gcmap(f'{out}/{di}_{resolution // 1000 }kb_raw_from_text.gcmap', xlabels=f'{di}',coarseningMethod='sum', compression='lzf')

        ## normalize with KR
        nrom_KR_gcmap = f'{KR_out}/{di}_{resolution // 1000 }kb_raw_from_text_KR.gcmap'
        gmlib.normalizer.normalizeGCMapByKR(f'{out}/{di}_{resolution // 1000 }kb_raw_from_text.gcmap', nrom_KR_gcmap, tol=1e-10)

        ##normalize with MCFs
        os.system(
            f'gcMapExplorer normMCFS -i "{KR_out}/{di}_{resolution // 1000 }kb_raw_from_text_KR.gcmap" -fi gcmap -o "{out}/{di}_{resolution // 1000}kb_raw_from_text_MCFS.gcmap" -fo gcmap')

        filename = f'{out}/{di}_{resolution // 1000}kb_raw_from_text_MCFS.gcmap'
        gcmap = gmlib.gcmap.GCMAP(filename, mapName=f'{di}')
        print(len(gcmap.matrix[:]))


        df_fra = pd.read_csv(path_fragment + '/' + di , sep = '\t',header = None)
        diff = abs(len(df_fra)-len(gcmap.matrix[:]))
        my_res = pd.DataFrame(np.pad(gcmap.matrix[:], ((0, diff), (0, diff)), 'constant', constant_values=0))
        print(len(df_fra),len(my_res))
        my_res.to_csv(f'2.KR_MCFS_normalized_dense_matrix/{di}_dense', sep='\t', header=None, index=False)





