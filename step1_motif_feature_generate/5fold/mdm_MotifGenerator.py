# coding=utf-8
import copy
import numpy as np
import sys
import getopt
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


d_related_m = []
for j in range(nd):
    temp = []
    for i in range(nm):
        if AA[i][j] == 1:
            temp.append(i)
    d_related_m.append(temp)

motifs = []
for j in range(len(d_related_m)):
    for m1 in range(len(d_related_m[j]) - 1):
        for m2 in range(m1+1, len(d_related_m[j])):
            temp = [d_related_m[j][m1], j+nm, d_related_m[j][m2]]
            motifs.append(temp)


output_file = 'mdm_motif_SEED_{0}_FOLD_{1}'.format(SEED, FOLD)
            
output = open(output_file, 'w')
for line in motifs:
    line = [str(x) for x in line]
    line = '%s' % '\t'.join(line)
    output.write(line + "\n")
output.close()
