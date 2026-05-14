import numpy as np
import sklearn as sk

def SPLIT_DATASET(DATASET, SEED, TEST_SIZE, TxID):

    X = DATASET.loc[:, (DATASET.columns != TxID)]
    y = DATASET.loc[:, TxID]

    X_train, X_test, \
    y_train, y_test = sk.model_selection.train_test_split(X, y, test_size = TEST_SIZE, random_state = SEED)

    return([X_train, X_test, y_train, y_test])
