import torch
import torch.distributions as dist
import numpy as np
import scipy.sparse as ssp
import networkx as nx
import random
import pandas as pd
import utils
import os

from create_factorial_table import create_factorial_table
from PwcmpASAPScheduler import PwcmpASAPScheduler

# % Example of using PwcmpASAPScheduler

# % We have 3 scenes (images or contents), each processed by two methods and
# % each method is run using three different parameter values.
scene = ['A', 'B', 'C']
method = ['methodA', 'methodB']
param = [ 0, 1, 2 ]

sets = {
    'scene': scene, 
    'method': method, 
    'param': param
    }

# % Create a table with all combinations of the elements from the sets
condition_table = create_factorial_table( sets )

# % We want to compare relative quality within each scene: comapre all
# % methods at all parameter values with each other, but without any
# % cross-scene comparisons. 
sch = PwcmpASAPScheduler( 'test_results.csv', 'test_obs3', condition_table, ['scene'] )

N_batch = sch.get_pair_left()

for kk in range(N_batch): # % for each pairwise comparison in the batch
    
    stim_A, stim_B = sch.get_next_pair() # % Get the next pair to compare
    
    print( 'Compare {} with {}\n'.format(stim_A, stim_B) )
    # % display( condition_table([stim_A stim_B],:) );
    
    # % Simulate the response. In actual experiment, this is the place in which
    # % the observer is asked to choose one of the two conditions.
    i_a = np.where( condition_table['method'][stim_A] == method )
    i_b = np.where( condition_table['method'][stim_B] == method )    
    p_sel = np.random.normal()
    is_A_selected = np.random.normal()<=p_sel
    
    # % Record the result of the comparison. The result will be immediately
    # % written to the file. 
    sch.set_pair_result( is_A_selected )
