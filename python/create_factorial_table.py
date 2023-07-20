import numpy as np
import pandas as pd

def create_factorial_table(sets):
    # % Create a table with all combinations of passed sets
    # %
    # % T = create_factorial_table( set1, set2, ... )
    # % 
    # % Each set must be a named cell array or a numeric 1D array. The name of
    # % the passed variable will be used as the name of the column.
    # %
    # % Example: 
    # % set1 = { 'A', 'B' };
    # % set2 = [1 2];
    # % create_factorial_table( {'set1'=set1, 'set2'=set2} )

    X, Y, Z  = np.meshgrid(*sets.values())
    res = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    
    d = {}

    keys = [*sets]
    for i in range(len(sets)):
        d[str(keys[i])] = pd.Series(res[:, i], dtype=type(sets[keys[i]][0]))
    
    return pd.DataFrame(d)
