import torch
import torch.distributions as dist
import numpy as np
import scipy.sparse as ssp
import networkx as nx
import random
import pandas as pd

def columns_equal(D, columns, value):
    # % Return boolean vector with 1s for all rows in table D for which all
    # % values in 'columns' match 'values'. 

    N = D.shape[0]
    ss = np.zeros(N)

    for c in columns:
        blocks = D[c].to_numpy()
        idx = np.where(np.isin(blocks, value))[0]
        ss[idx] = 1

    return ss