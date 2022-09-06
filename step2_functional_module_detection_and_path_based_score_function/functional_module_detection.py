# coding=utf-8
import time
import copy
import numpy as np
import pandas as pd
import math
from sklearn.cluster import KMeans
import utility


SEED = 0
nm = 495  # number of miRNAs
nd = 383  # number of diseases
nc = 5430  # number of miRNA-disease associations
nn = 495 * 383 - 5430  # number of unknown samples
ConnectDate = np.loadtxt('../data/known disease-miRNA association number ID.txt', dtype=int) - 1
A = np.zeros((nm, nd))
for i in range(nc):
    A[ConnectDate[i, 0], ConnectDate[i, 1]] = 1


T = utility.T
DS1 = np.loadtxt('../data/disease semantic similarity matrix 1.txt')
DS2 = np.loadtxt('../data/disease semantic similarity matrix 2.txt')
DS = (DS1+DS2)/2
FS = np.loadtxt('../data/miRNA functional similarity matrix.txt')

input_motif_list = ['triangle', 'triangle2miRNA', 'triangle2disease', 'mdm', 'mddd']

def func_(FOLD):
    print('now process: SEED:', SEED, ' FOLD:', FOLD)
    start = time.perf_counter()
    tran_matrix_file_template = '../data/5_fold_motifs/{0}_step2_SEED_' + str(SEED) + '_FOLD_' + str(FOLD)
    dic_file_template = '../data/5_fold_motifs/{0}_step1_SEED_' + str(SEED) + '_FOLD_' + str(FOLD) + '_dic'

    def initial_M_matrix(array):
        for i in range(len(array)):
            for j in range(len(array)):
                array[i][j] = float(array[i][j])

        tran_matrix = np.array(array)
        for i in range(tran_matrix.shape[0]):
            _sum = np.sum(tran_matrix[i])
            for j in range(tran_matrix.shape[1]):
                if tran_matrix[i, j] > 0:
                    tran_matrix[i, j] = tran_matrix[i, j] / _sum
        return tran_matrix
    
    
    def transferMto_M(M, dic):
        connects = np.argwhere(M > 0)
        _M = np.zeros((nm + nd, nm + nd))
        for i in connects:
            _M[dic[i[0]]][dic[i[1]]] = M[i[0]][i[1]]
        return _M
    
    
    def Getgauss_miRNA(adjacentmatrix,nm):
        """
        MiRNA Gaussian interaction profile kernels similarity
        """
        KM = np.zeros((nm,nm))
        gamaa=1
        sumnormm=0
        for i in range(nm):
            normm = np.linalg.norm(adjacentmatrix[i])**2
            sumnormm = sumnormm + normm  
        gamam = gamaa/(sumnormm/nm)
        for i in range(nm):
            for j in range(nm):
                KM[i,j]= math.exp (-gamam*(np.linalg.norm(adjacentmatrix[i]-adjacentmatrix[j])**2))
        return KM
           
    def Getgauss_disease(adjacentmatrix,nd):
        """
        Disease Gaussian interaction profile kernels similarity
        """
        KD = np.zeros((nd,nd))
        gamaa=1
        sumnormd=0
        for i in range(nd):
            normd = np.linalg.norm(adjacentmatrix[:,i])**2
            sumnormd = sumnormd + normd
        gamad=gamaa/(sumnormd/nd)
        for i in range(nd):
            for j in range(nd):
                KD[i,j]= math.exp(-(gamad*(np.linalg.norm(adjacentmatrix[:,i]-adjacentmatrix[:,j])**2)))
        return KD

    def judge_edge_in_cluster(node1, node2, clusters):
        res = False
        for i in clusters:
            if (node1 in i) and (node2 in i):
                res = True
                break
        return res

    AA = copy.deepcopy(A)
    train_id_all, test_id_all = utility.tradition_5fold_index(SEED, FOLD)
    for i in range(len(test_id_all)):
        AA[test_id_all[i][0]][test_id_all[i][1]] = 0
    
    KM = Getgauss_miRNA(AA,nm)
    KD = Getgauss_disease(AA,nd)
    
    # integrating miRNA functional similarity and Gaussian interaction profile kernels similarity
    FS_integration=np.zeros((nm,nm))
    for i in range(nm):
        for j in range(nm):
            if FS[i,j] > 0:    
                FS_integration[i,j] = FS[i,j]
            else:
                FS_integration[i,j] = KM[i,j]
      
    # integrating disease semantic similarity and Gaussian interaction profile kernels similarity
    DS_integration=np.zeros((nd,nd))
    for i in range(nd):
        for j in range(nd):
            if DS[i,j] > 0:   
                DS_integration[i,j] = DS[i,j]
            else:
                DS_integration[i,j] = KD[i,j]

    AM = np.zeros((nm+nd, nm+nd), dtype=float)
    
    for i in range(nm):
        for j in range(nm):
            if FS_integration[i,j] >= T:
                AM[i][j] = FS_integration[i,j]
                AM[j][i] = FS_integration[i,j]
    
    for i in range(nd):
        for j in range(nd):
            if DS_integration[i,j] >= T:
                AM[i+nm][j+nm] = DS_integration[i,j]
                AM[j+nm][i+nm] = DS_integration[i,j]
    
    for i in range(nm):
        for j in range(nd):
            if AA[i][j] != 0:
                AM[i][j+nm] = 1
                AM[j+nm][i] = 1
    
    dic_motif_to__M = dict()
    for i in input_motif_list:
        f = open(dic_file_template.format(i))
        dic_index_to_originNum = dict()
        for line in f:
            items = line.split(" ")
            dic_index_to_originNum[int(items[1]) - 1] = int(items[0])
        f.close()
    
        f = open(tran_matrix_file_template.format(i))
        array = [line.split("\t") for line in f]
        M = initial_M_matrix(array)
        f.close()
        dic_motif_to__M[i] = transferMto_M(M, dic_index_to_originNum)
    
    
    for i in range(AM.shape[0]):
        AM[i][i] = 0
    AM2 = np.zeros((nm + nd, nm + nd), dtype=float)
    for i in range(AM.shape[0]):
        _sum = np.sum(AM[i])
        for j in range(AM.shape[1]):
            if AM[i, j] > 0:
                AM2[i, j] = AM[i, j] / _sum
    
    AM_weight = np.zeros((nm+nd, nm+nd), dtype=float)
    count = np.zeros((nm+nd, nm+nd), dtype=float)
    for motif in input_motif_list:
        AMs = np.zeros((nm + nd, 2*(nm + nd)), dtype=float)
        for i in range(nm+nd):
            for j in range(nm+nd):
                AMs[i][j] = dic_motif_to__M[motif][i][j]
                AMs[i][j+nm+nd] = AM2[i][j]
        K = 4
        kmeans = KMeans(n_clusters=K,
                        init='k-means++',
                        n_init=100,
                        max_iter=5000,
                        tol=0.0001,
                        random_state=2020,
                        algorithm='auto')
        M = pd.DataFrame(AMs)
        kmeans.fit(M)
        results = kmeans.labels_
        clus = []
        for _ in range(K):
            clus.append([])
        for idx, c in enumerate(results):
            clus[c].append(idx)
        clus_sorted = sorted(clus, key=lambda i:len(i), reverse=True)
        for i in range(nm+nd):
            for j in range(nm+nd):
                if ((dic_motif_to__M[motif][i][j] > 0) or (AM2[i][j] > 0)) and judge_edge_in_cluster(i, j, clus_sorted):
                    AM_weight[i][j] += dic_motif_to__M[motif][i][j] + AM2[i][j]
                    count[i][j] += 1

    for i in range(nm+nd):
        for j in range(nm+nd):
            if count[i][j] > 1:
                AM_weight[i][j] = AM_weight[i][j] / count[i][j]
    
    np.savetxt('../data/AM_weight/AM_weight_SEED_{0}_FOLD_{1}'.format(SEED, FOLD), AM_weight)
    end = time.perf_counter()
    print('SEED:', SEED, ' FOLD:', FOLD, 'cost time:', end-start)
