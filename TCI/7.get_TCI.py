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

def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a - d, np.cross(b - d, c - d))) / 6


def concave_hull_volume(pts, lenCutoff):
    dt = Delaunay(pts)  # tetrahedrons formed up with 4 points
    tets = dt.points[dt.simplices]
    max_length = [0 for i in range(len(tets))]
    min_length = [0 for i in range(len(tets))]
    min2max = [0 for i in range(len(tets))]

    vol = 0
    count_valid_dts = 0
    count_total_dts = len(tets)

    total_vol = tetrahedron_volume(tets[:, 0], tets[:, 1], tets[:, 2], tets[:, 3])
    for i in range(len(tets)):

        max_length[i] = max(np.linalg.norm(tets[i, 0] - tets[i, 1]), np.linalg.norm(tets[i, 0] - tets[i, 2]),
                            np.linalg.norm(tets[i, 0] - tets[i, 3]), np.linalg.norm(tets[i, 1] - tets[i, 2]),
                            np.linalg.norm(tets[i, 1] - tets[i, 3]), np.linalg.norm(tets[i, 2] - tets[i, 3]))
        min_length[i] = min(np.linalg.norm(tets[i, 0] - tets[i, 1]), np.linalg.norm(tets[i, 0] - tets[i, 2]),
                            np.linalg.norm(tets[i, 0] - tets[i, 3]), np.linalg.norm(tets[i, 1] - tets[i, 2]),
                            np.linalg.norm(tets[i, 1] - tets[i, 3]), np.linalg.norm(tets[i, 2] - tets[i, 3]))
        min2max[i] = min_length[i] / float(max_length[i])

        if max_length[i] > lenCutoff:
            vol += 0
        else:
            vol += total_vol[i]  # tets[:,0-4] are the four points
            count_valid_dts += 1

    if count_total_dts == 0 or vol == 0 or len(min2max) == 0:
        print(pts)
        return 0, 0, 0, 0

    return vol, len(pts) / float(vol), float(sum(min2max)) / len(min2max), float(
        count_valid_dts) / count_total_dts  # point density


def get_density(pos):
    lenCutoff = 1000
    vol, density, min2max, validprop = concave_hull_volume(pos, lenCutoff)
    return vol, density

def simulate_random_walk_volume(total_loop,persistence_length, kuhn_segment_length, num_steps):
    temp = []
    for l in range(1,total_loop+1): 
        polymer_chain = np.zeros((num_steps, 3))
        for step in range(1, num_steps):
            step_length = kuhn_segment_length
            step_direction = np.random.normal(size=3)
            step_direction /= np.linalg.norm(step_direction)
            step_direction *= step_length
            step_direction = (1 - 1/persistence_length) * step_direction + \
                             (1/persistence_length) * polymer_chain[step-1, :]
            polymer_chain[step, :] = polymer_chain[step-1, :] + step_direction
        vol, density = get_density(polymer_chain)
        temp.append(vol)
    return max(temp)

persistence_length = 50 
kuhn_segment_length = 2 * persistence_length  # Kuhn segment 
total_loop = 100

df_coordinate = pd.read_csv('XYZ/all_coordinate',sep = '\t')
df_tad = pd.read_csv('0.NHEK_All_chr_TAD',sep = '\t',header = None)

dic_coordiante = {}
mkdir('3.each_TAD_coordinate')
for i in range(len(df_coordinate)):
    if df_coordinate.iloc[i,0] not in dic_coordiante:
        dic_coordiante[df_coordinate.iloc[i,0]] = []
        dic_coordiante[df_coordinate.iloc[i, 0]].append(df_coordinate.iloc[i,:].values.tolist())
    else:
        dic_coordiante[df_coordinate.iloc[i, 0]].append(df_coordinate.iloc[i, :].values.tolist())

dic_tad = {}
for i in range(len(df_tad)):
    if df_tad.iloc[i,0] not in dic_tad:
        dic_tad[df_tad.iloc[i,0]] = []
        dic_tad[df_tad.iloc[i, 0]].append(df_tad.iloc[i,:].values.tolist())
    else:
        dic_tad[df_tad.iloc[i, 0]].append(df_tad.iloc[i, :].values.tolist())

df_coordinate.index = df_coordinate.apply(lambda x:'_'.join(map(str,x[:3])),axis = 1)
for k in dic_tad:
    for v in dic_tad[k]:
        tad_start = v[1]
        tad_end = v[2]
        if ((tad_end-tad_start)/10000) < 5 : 
            continue
        tad_coordinate = []
        save_tad_coordinate = []
        for v_c in dic_coordiante[k]:
            coor_start = v_c[1]
            coor_end = v_c[2]
            if int(tad_start)+1 <= int(coor_start) and int(tad_end) >= int(coor_end):
                tad_coordinate.append(v_c[4:])
                save_tad_coordinate.append([v_c[i] for i in range(len(v_c)) if i not in [3]])
        savetxt(f'3.each_TAD_coordinate/{"_".join(map(str, v))}', save_tad_coordinate)

        if len(tad_coordinate) == 0 :
            v.append(str(0))
            continue
        tad_coordinate_array = np.array(tad_coordinate)
        if ((tad_end-tad_start)/10000)*0.8 > len(tad_coordinate_array):
            v.append(str(0))
            continue
        vol, density = get_density(tad_coordinate_array)
        v.append(str(vol))
        num_steps = (v[2] - v[1]) // 10000
        polymer_chain_vol = simulate_random_walk_volume(total_loop,persistence_length, kuhn_segment_length, num_steps)
        v.append(str(polymer_chain_vol))
        nor = float(vol)/float(polymer_chain_vol)
        v.append(str(float(vol)/float(polymer_chain_vol)))

        hull = ConvexHull(tad_coordinate_array)
        hull_vertices = tad_coordinate_array[hull.vertices]
        center = np.mean(hull_vertices, axis=0)
        point1 = center
        distance = []
        for ar in tad_coordinate_array:
            point2 = ar
            ojld = math.sqrt((point2[0] - point1[0]) ** 2 + (point2[1] - point1[1]) ** 2 + (point2[2] - point1[2]) ** 2)
            distance.append(ojld)
        average_distance = sum(distance)/len(tad_coordinate_array)
        v.append(str(average_distance))

res = []
res.append('chr'+'\t'+'start'+'\t'+'end'+'\t'+'TAD_volume'+'\t'+'polymer_volume'+'\t'+'TCI'+'\t'+'average_ojld_distance')
for key in dic_tad:
    for value in dic_tad[key]:
        res.append('\t'.join(map(str,value)))

savetxt('3.TAD_volume',res)






