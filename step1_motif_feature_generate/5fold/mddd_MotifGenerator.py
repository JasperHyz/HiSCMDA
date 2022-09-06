# coding=utf-8
import copy
import numpy as np
import sys
import getopt
import math
import utility


SEED = 0
FOLD = 0
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


m_related_d = []
for i in range(nm):
    temp = []
    for j in range(nd):
        if AA[i][j] == 1:
            temp.append(j)
    m_related_d.append(temp)


motifs = []
for i in range(nm):
    for d1 in range(len(m_related_d[i]) - 2):
        for d2 in range(d1 + 1, len(m_related_d[i]) - 1):
            for d3 in range(d2 + 1, len(m_related_d[i])):
                motifs.append([m_related_d[i][d3]+nm, m_related_d[i][d2]+nm, m_related_d[i][d1]+nm, i])



output_file = 'mddd_motif_SEED_{0}_FOLD_{1}'.format(SEED, FOLD)
            
output = open(output_file, 'w')
for line in motifs:
    line = [str(x) for x in line]
    line = '%s' % '\t'.join(line)
    output.write(line + "\n")
output.close()
