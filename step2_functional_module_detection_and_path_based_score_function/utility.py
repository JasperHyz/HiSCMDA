import numpy as np
import pandas as pd
import scipy.sparse as sp
from copy import deepcopy
from sklearn.model_selection import KFold
import random

# threshold beta
T = 0.6



def tradition_5fold_index(seed, fold):
    random.seed(seed)
    np.random.seed(seed)
    all_associations = pd.read_csv('../data/all_miRNA_disease_pairs.csv')
    known_associations = all_associations.loc[all_associations['label'] == 1]
    sample_df = known_associations
    sample_df.reset_index(drop=True, inplace=True)

    kf = KFold(n_splits=5, shuffle=True, random_state=seed)
    train_index_all, test_index_all = [], []
    train_id_all, test_id_all = [], []
    for train_idx, test_idx in kf.split(sample_df):
        train_index_all.append(train_idx) 
        test_index_all.append(test_idx)

        train_id_all.append(np.array(sample_df.iloc[train_idx][['miRNA', 'disease', 'label']]))
        test_id_all.append(np.array(sample_df.iloc[test_idx][['miRNA', 'disease', 'label']]))
    # fold 1-5
    # 0:miRNA, 1:disease, 2:label; miRNA,disease begin with 0
    return train_id_all[fold-1], test_id_all[fold-1]


def balanced_5fold_index(seed, fold):
    random.seed(seed)
    np.random.seed(seed)

    all_associations = pd.read_csv('../data/all_miRNA_disease_pairs.csv')
    known_associations = all_associations.loc[all_associations['label'] == 1]
    unknown_associations = all_associations.loc[all_associations['label'] == 0]
    random_negative = unknown_associations.sample(n=known_associations.shape[0], random_state=seed, axis=0)
    sample_df = known_associations.append(random_negative)
    sample_df.reset_index(drop=True, inplace=True)

    kf = KFold(n_splits=5, shuffle=True, random_state=seed)
    train_index_all, test_index_all = [], []
    train_id_all, test_id_all = [], []
    for train_idx, test_idx in kf.split(sample_df):
        train_index_all.append(train_idx) 
        test_index_all.append(test_idx)

        train_id_all.append(np.array(sample_df.iloc[train_idx][['miRNA', 'disease', 'label']]))
        test_id_all.append(np.array(sample_df.iloc[test_idx][['miRNA', 'disease', 'label']]))
    # fold 1-5
    # 0:miRNA, 1:disease, 2:label; miRNA,disease begin with 0
    return train_id_all[fold-1], test_id_all[fold-1]
