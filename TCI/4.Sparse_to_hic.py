import numpy as np
import os
import pandas as pd

def savetxt(filename, x):
    np.savetxt(filename, x, delimiter=' ', fmt='%s')

chrname = []
for i in range(1,23):
    chrname.append("chr"+str(i))
chrname.append("chrX")

resolution = 10000

file_path = "3.dense2sparse/"
dirs = os.listdir(file_path)

os.system('mkdir -p 4.sparse_2_hic')
for di in dirs:
    chr = di.split('_')[0]
    print(chr)
    os.system(f'mkdir -p 4.sparse_2_hic/{chr}')
    f = file_path + '/' + di
    os.system(f'hicConvertFormat -m {f} -bf /mnt/w/DSB_TAD_density/call_ISTAD/hic_KR_sparse2dense/NHEK/2.1each_chr_10_fragment/{chr} --inputFormat hicpro --outputFormat cool -o 4.sparse_2_hic/{chr}/{chr}_matrix.cool')
    os.system(f'hicConvertFormat -m 4.sparse_2_hic/{chr}/{chr}_matrix.cool -o 4.sparse_2_hic/{chr}/{chr}.ginteractions --inputFormat cool --outputFormat ginteractions')
    command = f"""awk -F "\\t" '{{print 0, $1, $2, 1, 1, $4, $5, 1, $7}}' 4.sparse_2_hic/{chr}/{chr}.ginteractions.tsv > 4.sparse_2_hic/{chr}/{chr}.ginteractions.tsv.short"""
    os.system(command)
    os.system(f'sort -k2,2d -k6,6d 4.sparse_2_hic/{chr}/{chr}.ginteractions.tsv.short > 4.sparse_2_hic/{chr}/{chr}.ginteractions.tsv.short.sorted')
    print('###########################generate .hic#########################')
    os.system(f'java -Xmx20g -jar juicer_tools_1.22.01.jar pre -r 10000,20000,50000,100000,250000,500000,1000000 4.sparse_2_hic/{chr}/{chr}.ginteractions.tsv.short.sorted 4.sparse_2_hic/{chr}/{chr}.hic hg19')


