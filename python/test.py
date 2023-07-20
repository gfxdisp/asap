import numpy as np
from create_factorial_table import create_factorial_table

scene = ['A', 'B', 'C']
method = ['methodA', 'methodB']
param = [ 0, 1, 2 ]

sets = {
    'scene': scene, 
    'method': method, 
    'param': param
    }

df = create_factorial_table(sets)
print(df, df.dtypes)
