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

def simulate_random_walk_volume(total_loop, persistence_length, kuhn_segment_length, num_steps, cache=None):
    if cache is not None and num_steps in cache:
        return cache[num_steps]

    temp = []
    for l in range(1, total_loop + 1): 

        polymer_chain = np.zeros((num_steps, 3))

        for step in range(1, num_steps):

            step_length = kuhn_segment_length
            step_direction = np.random.normal(size=3)
            step_direction /= np.linalg.norm(step_direction)
            step_direction *= step_length

            step_direction = (1 - 1 / persistence_length) * step_direction + \
                            (1 / persistence_length) * polymer_chain[step - 1, :]


            polymer_chain[step, :] = polymer_chain[step - 1, :] + step_direction
        vol, density = get_density(polymer_chain)
        temp.append(vol)

    max_vol = max(temp)
    if cache is not None:
        cache[num_steps] = max_vol

    return max_vol


persistence_length = 50  
kuhn_segment_length = 2 * persistence_length  # Kuhn segment 长度
total_loop = 10

path= '2.allDSB/3.run_FLAMINGO_top'
dirs = os.listdir(path)
res = []
res.append('chr' + '\t' + 'start' + '\t' + 'end' + '\t' + 'TAD_volume' + '\t' + 'polymer_volume' + '\t' + 'TCI' + '\t' + 'average_ojld_distance' )

polymer_volume_cache = {}
for di in dirs:
    if len(os.listdir(path + '/' + di)) !=2:
        continue
    f = path + '/' + di + '/coordinate'
    df = pd.read_csv(f,sep = '\t',header = None,skiprows=1)
    df_coor = df.iloc[:,1:].values.tolist()
    coor_array = np.array(df_coor)
    if len(coor_array) <= 4:
        continue

    vol, density = get_density(coor_array)

    start = int(di.split('_')[1])
    end = int(di.split('_')[2])
    num_steps = (end - start) // 10000
    if num_steps in polymer_volume_cache:
        polymer_chain_vol = polymer_volume_cache[num_steps]
    else:
        polymer_chain_vol = simulate_random_walk_volume(total_loop, persistence_length, kuhn_segment_length, num_steps,
                                                        cache=polymer_volume_cache)


    nor = float(vol) / float(polymer_chain_vol)



    hull = ConvexHull(coor_array)
    hull_vertices = coor_array[hull.vertices]
    center = np.mean(hull_vertices, axis=0)
    point1 = center
    distance = []
    for ar in coor_array:
        point2 = ar
        ojld = math.sqrt((point2[0] - point1[0]) ** 2 + (point2[1] - point1[1]) ** 2 + (point2[2] - point1[2]) ** 2)
        distance.append(ojld)
    average_distance = sum(distance) / len(coor_array)

    res.append(di.split('_')[0]+'\t'+di.split('_')[1]+'\t'+di.split('_')[2] + '\t' + str(vol) + '\t' + str(polymer_chain_vol) + '\t' + str(nor)+'\t'+str(average_distance))

savetxt('3.TAD_volume_sameTLsameRW', res)




