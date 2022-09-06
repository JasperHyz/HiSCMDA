# coding=utf-8
import time
import copy
import numpy as np
import pandas as pd
import math
import os
import sys
import getopt
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

def func_(FOLD):
    def judge_node1_miRNA_node2_disease(node1, node2):
        if ((node1 < nm) and (node2 < nm)) or ((node1 >= nm) and (node2 >= nm)):
            return False
        elif node1 < nm and node2 >= nm:
            return 1
        else:
            return 2


    def l1_path_score(node1, node2, AM, score):
        if judge_node1_miRNA_node2_disease(node1, node2) > 0:
            if judge_node1_miRNA_node2_disease(node1, node2) == 1:
                p1 = AM[node1][node2]
                score[node1][node2 - nm] += p1
            else:
                p1 = AM[node1][node2]
                score[node2][node1 - nm] += p1


    def l2_path_score(node1, node2, node3, AM, score):
        if judge_node1_miRNA_node2_disease(node1, node3) > 0:
            if judge_node1_miRNA_node2_disease(node1, node3) == 1:
                p1 = AM[node1][node2]
                p2 = AM[node2][node3]
                score[node1][node3 - nm] += p1 * p2
            else:
                p1 = AM[node1][node2]
                p2 = AM[node2][node3]
                score[node3][node1 - nm] += p1 * p2


    def l3_path_score(node1, node2, node3, node4, AM, score):
        if judge_node1_miRNA_node2_disease(node1, node4) > 0:
            if judge_node1_miRNA_node2_disease(node1, node4) == 1:
                p1 = AM[node1][node2]
                p2 = AM[node2][node3]
                p3 = AM[node3][node4]
                score[node1][node4 - nm] += p1 * p2 * p3
            else:
                p1 = AM[node1][node2]
                p2 = AM[node2][node3]
                p3 = AM[node3][node4]
                score[node4][node1 - nm] += p1 * p2 * p3

    print('now process: SEED:', SEED, ' FOLD:', FOLD)
    start = time.perf_counter()
    score = np.zeros((nm, nd), dtype=float)
    AM_prob = np.loadtxt('../data/AM_weight/AM_weight_SEED_{0}_FOLD_{1}'.format(SEED, FOLD))

    for i in range(AM_prob.shape[0]):
        AM_prob[i][i] = 0

    AM_prob_related = []
    for i in range(nm+nd):
        temp = []
        for j in range(nm+nd):
            if AM_prob[i][j] != 0:
                temp.append(j)
        AM_prob_related.append(temp)

    for i in range(nm+nd):
        for j in AM_prob_related[i]:
            if i == j:
                continue
            l1_path_score(i, j, AM_prob, score)
            for l in AM_prob_related[j]:
                if l == j or l == i:
                    continue
                l2_path_score(i, j, l, AM_prob, score)
                for t in AM_prob_related[l]:
                    if t == i or t == j or t == l:
                        continue
                    l3_path_score(i, j, l, t, AM_prob, score)

    train_id_all, test_id_all = utility.tradition_5fold_index(SEED, FOLD)
    dataset_n = np.argwhere(A == 0)
    m_d_score = np.zeros((len(dataset_n) + len(test_id_all), 3),dtype=float)

    for i in range(len(test_id_all)):
        m_d_score[i][0] = test_id_all[i][0]
        m_d_score[i][1] = test_id_all[i][1]
        m_d_score[i][2] = score[test_id_all[i][0]][test_id_all[i][1]]

    for i in range(len(dataset_n)):
        m_d_score[i+len(test_id_all)][0] = dataset_n[i][0]
        m_d_score[i+len(test_id_all)][1] = dataset_n[i][1]
        m_d_score[i+len(test_id_all)][2] = score[dataset_n[i][0]][dataset_n[i][1]]
    df = pd.DataFrame(m_d_score)
    probability = np.zeros((nm, nd),dtype=float)
    df['rank'] = df[2].rank(ascending=0, method='max')
    n = df.shape[0]
    for i in range(0, n):
        probability[int(df.iloc[i][0])][int(df.iloc[i][1])] = 1 - (df.iloc[i]['rank'] - 1) / (n - 1)
    np.savetxt('../data/result/prob_SEED_{0}_FOLD_{1}'.format(SEED, FOLD), probability)
    np.savetxt('../data/result/score_SEED_{0}_FOLD_{1}'.format(SEED, FOLD), score)
    end = time.perf_counter()
    print('SEED:', SEED, ' FOLD:', FOLD, 'cost time:', end-start)



