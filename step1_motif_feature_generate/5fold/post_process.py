# coding=utf-8
import sys
import getopt
import os
import copy
import numpy as np
import random
import pandas as pd
import time
import shutil


SEED = 0
FOLD = 0
MOTIF = ''
try:
    opts, args = getopt.getopt(sys.argv[1:], "s:f:m:", ['SEED=', 'FOLD=', 'MOTIF='])
except getopt.GetoptError:
    print('exception')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-s', '--SEED'):
        SEED = int(arg)
    if opt in ('-f', '--FOLD'):
        FOLD = int(arg)
    if opt in ('-m', '--MOTIF'):
        MOTIF = arg


dic_file = '{0}_step1_SEED_{1}_FOLD_{2}_dic'.format(MOTIF, SEED, FOLD)
tran_matrix_file = '{0}_step2_SEED_{1}_FOLD_{2}'.format(MOTIF, SEED, FOLD)



os.remove('{0}_motif_SEED_{1}_FOLD_{2}'.format(MOTIF, SEED, FOLD))
os.remove('{0}_step1_SEED_{1}_FOLD_{2}'.format(MOTIF, SEED, FOLD))

shutil.move(dic_file, '../../data/5_fold_motifs/'+dic_file)
shutil.move(tran_matrix_file, '../../data/5_fold_motifs/'+tran_matrix_file)
