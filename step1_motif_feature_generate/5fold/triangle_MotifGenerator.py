# coding=utf-8
import copy
from msilib.schema import ListBox
import numpy as np
import sys
import getopt
import math
import utility

SEED = 0
FOLD = 0
T = utility.T
try:
    opts, args = getopt.getopt(sys.argv[1:], "s:f:", ['SEED=', 'FOLD='])
except getopt.GetoptError:
    print('exception')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-s', '--SEED'):
        SEED = int(arg)
    if opt in ('-f', '--FOLD'):
        FOLD = int(arg)
nm = 495  # number of miRNAs
nd = 383  # number of diseases
nc = 5430  # number of miRNA-disease associations
ConnectDate = np.loadtxt('../../data/known disease-miRNA association number ID.txt', dtype=int) - 1
A = np.zeros((nm, nd))
for i in range(len(ConnectDate)):
    A[ConnectDate[i, 0], ConnectDate[i, 1]] = 1
AA = copy.deepcopy(A)
train_id_all, test_id_all = utility.tradition_5fold_index(SEED, FOLD)
for i in range(len(test_id_all)):
    AA[test_id_all[i][0]][test_id_all[i][1]] = 0

DS1 = np.loadtxt('../../data/disease semantic similarity matrix 1.txt')
DS2 = np.loadtxt('../../data/disease semantic similarity matrix 2.txt')
DS = (DS1 + DS2) / 2
FS = np.loadtxt('../../data/miRNA functional similarity matrix.txt')


def Getgauss_miRNA(adjacentmatrix, nm):
    """
    MiRNA Gaussian interaction profile kernels similarity
    """
    KM = np.zeros((nm, nm))

    gamaa = 1
    sumnormm = 0
    for i in range(nm):
        normm = np.linalg.norm(adjacentmatrix[i]) ** 2
        sumnormm = sumnormm + normm
    gamam = gamaa / (sumnormm / nm)

    for i in range(nm):
        for j in range(nm):
            KM[i, j] = math.exp(-gamam * (np.linalg.norm(adjacentmatrix[i] - adjacentmatrix[j]) ** 2))
    return KM


def Getgauss_disease(adjacentmatrix, nd):
    """
    Disease Gaussian interaction profile kernels similarity
    """
    KD = np.zeros((nd, nd))
    gamaa = 1
    sumnormd = 0
    for i in range(nd):
        normd = np.linalg.norm(adjacentmatrix[:, i]) ** 2
        sumnormd = sumnormd + normd
    gamad = gamaa / (sumnormd / nd)

    for i in range(nd):
        for j in range(nd):
            KD[i, j] = math.exp(-(gamad * (np.linalg.norm(adjacentmatrix[:, i] - adjacentmatrix[:, j]) ** 2)))
    return KD

def __get_com_idx_arr(arr_1, arr_2):
    n = len(arr_1)
    arr = []

    for i in range(0, n):
        if arr_1[i] and arr_2[i]:
            arr.append(i)

    return arr

KM = Getgauss_miRNA(AA,nm)
KD = Getgauss_disease(AA,nd)
#integrating miRNA functional similarity and Gaussian interaction profile kernels similarity
FS_integration=np.zeros((nm,nm))
for i in range(nm):
    for j in range(nm):
        if FS[i,j] > 0:
            FS_integration[i,j] = FS[i,j]
        else:
            FS_integration[i,j] = KM[i,j]

#integrating disease semantic similarity and Gaussian interaction profile kernels similarity
DS_integration=np.zeros((nd,nd))
for i in range(nd):
    for j in range(nd):
        if DS[i,j] > 0:
            DS_integration[i,j] = DS[i,j]
        else:
            DS_integration[i,j] = KD[i,j]

for i in range(nm):
    FS_integration[i][i] = 0
for i in range(nd):
    DS_integration[i][i] = 0


i_related_k = []
for i in range(nm):
    i_m = []
    for m in range(nm):
        if m == i: continue
        if FS_integration[i, m] >= T:
            i_m.append(m)
    i_related_k.append(i_m)

j_related_l = []
for j in range(nd):
    j_d = []
    for d in range(nd):
        if d == j: continue
        if DS_integration[j, d] >= T:
            j_d.append(d)
    j_related_l.append(j_d)
    
m_m = []
for i in range(len(i_related_k)):
    for j in i_related_k[i]:
        if i <= j:
            m_m.append([i, j])
        else:
            m_m.append([j, i])

m_m = np.array(m_m)
m_m = np.unique(m_m, axis=0)


AA_N = np.array(AA)

m_m_d = []
for i, j in m_m:
    m_m_d.append(__get_com_idx_arr(AA_N[i], AA_N[j]))

motifs = []
for i in range(len(m_m)):
    for j in m_m_d[i]:
        motifs.append([m_m[i][0], m_m[i][1], j + nm])
        
d_d = []
for i in range(len(j_related_l)):
    for j in j_related_l[i]:
        if i <= j:
            d_d.append([i, j])
        else:
            d_d.append([j, i])
d_d = np.array(d_d)
d_d = np.unique(d_d, axis=0)
d_d_m = []
for i, j in d_d:
    i_column = AA_N[:, i]
    j_column = AA_N[:, j]
    d_d_m.append(__get_com_idx_arr(i_column, j_column))

for i in range(len(d_d)):
    for j in d_d_m[i]:
        motifs.append([d_d[i][0] + nm, d_d[i][1] + nm, j])

# triangle of 3miRNA and 3disease
m_m_m = set()
for i, j in m_m:
    temp = [k for k in i_related_k[i] if k in i_related_k[j]]
    for x in temp:
        mtriangle_list = [i, j, x]
        mtriangle_list.sort()
        m_m_m.add((mtriangle_list[0], mtriangle_list[1], mtriangle_list[2]))

d_d_d = set()
for i, j in d_d:
    temp = [k for k in j_related_l[i] if k in j_related_l[j]]
    for x in temp:
        dtriangle_list = [i, j, x]
        dtriangle_list.sort()
        d_d_d.add((dtriangle_list[0], dtriangle_list[1], dtriangle_list[2]))


for i in m_m_m:
    motifs.append([i[0], i[1], i[2]])
for i in d_d_d:
    motifs.append([i[0] + nm, i[1] + nm, i[2] + nm])

output_file = 'triangle_motif_SEED_{0}_FOLD_{1}'.format(SEED, FOLD)
            
output = open(output_file, 'w')
for line in motifs:
    line = [str(x) for x in line]
    line = '%s' % '\t'.join(line)
    output.write(line + "\n")
output.close()
