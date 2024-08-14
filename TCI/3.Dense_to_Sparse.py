import pandas as pd
import numpy as np
import os
def mkdir(path):
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

def savetxt(filename, x):
    np.savetxt(filename, x, delimiter='\t', fmt='%s')

def sort_list(list_in):  # use when using '\t' to seperate items
    list_out = sorted(list_in, key=lambda items: int(items.split('\t')[1]))
    list_out = sorted(list_out, key=lambda items: items.split('\t')[0])
    return list_out

##################################################################################
path = '2.KR_MCFS_normalized_dense_matrix'
dirs = os.listdir(path)
mkdir('3.dense2sparse')
for di in dirs:
    f = path + '/' + di
    with open(f,'r')as file:
        lines = file.readlines()
        nrows = len(lines)
        matrix = np.zeros((nrows, nrows))
        if nrows < 5:
            continue
        for i in range(len(lines)):
            line = lines[i].strip().split('\t')
            matrix[i, :] = line[:]

    with open(f'3.dense2sparse/{di.split("_")[0]}_sparse','w') as fout:
        for i in range(nrows):
            for j in range(nrows):
                if i > j :
                    continue
                if matrix[i,j] ==0:
                    continue
                fout.write(str(i+1)+'\t'+str(j+1)+'\t'+str(matrix[i, j]) + '\n')

