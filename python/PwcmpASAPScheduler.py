import torch
import torch.distributions as dist
import numpy as np
import scipy.sparse as ssp
import networkx as nx
import random
import pandas as pd
import utils
import os
import asap_gpu

class PwcmpASAPScheduler:
    # % Class for scheduling pairwise comparisons using ASAP algorithm
    # %
    # % This class handles loading and saving results between experiments runs
    # % and also splits the comparison into blocks. Each block of measurements is
    # % done idependently with no comparisons made between the blocks. 
    # %
    # % Details on ASAP: 
    # %  https://github.com/gfxdisp/asap
    # %  
    def __init__(
        self,
        result_file,
        observer,
        condition_table,
        block_columns,
        
        cmp_M = [], # % 
        N_blocks = 0,
        
        cond_index = [], # % [block, cond] matrix with the indices of each (block,cond) pair in condition_table    
        
        # % List of condition pairs to compare, 1st column - block index, 2nd and
        # % 3d column is the pair of indices of conditions WITHIN the bl0ck (not global)
        pair_list = [], 
        rand_ord = [],
        pair_index = 0, # % next pair to compare
    ):
        # % Creates new ASAP scheduler and initializes with the mesurements
        # % collected so far. 
        # %
        # % result file - a CSV file with the measurememts. It will create a
        # %               new file if the file does not exist
        # % observer - a string with the ID of the participant
        # % condition table - a populated table data structure with one or
        # %            more columns that describe the dimensions of measured conditions.
        # %            For example: a table with columns (scene, method, parameter_A, parameter_B)
        # % block_columns - a cell array with the names of the columns that define separate
        # %            blocks of pairwise comparison measurememts. Each block
        # %            is measured independently, i.e. no comparisons are
        # %            made between the conditions in each block. If no
        # %            blocks should be created (all conditions should be
        # %            compared with each other), pass an empty cell {} or skip this
        # %            parameter. 
        # %            Example: ['scene'] -
        # %            the comparisons will be made only between the
        # %            conditions for which the 'scene' attribute is the same.         

        self.observer = observer
        self.result_file = result_file
        self.condition_table = condition_table
        self.block_columns = block_columns
        self.cond_index = cond_index
        self.pair_list = pair_list
        self.pair_index = pair_index
        
        if len(block_columns) == 0:
            block_table = []
            self.N_blocks = 1
        else:
            # block_table = ["scene_1", "scene_2", ...]
            block_table = condition_table[block_columns].drop_duplicates().to_numpy()
            self.N_blocks = block_table.shape[0]
        
        # % Populate cond_index so that we can quickly find a condition
        # % based on its index value
        self.cmp_M = []
        for kk in range(self.N_blocks):
            if len(block_columns) == 0:
                ss = np.ones((condition_table.shape[0], 1))
                Dss = condition_table
            else:
                ss = utils.columns_equal(condition_table, block_columns, block_table[kk])
                Dss = condition_table[ss.astype(np.bool)]
            
            self.cmp_M.append(np.zeros((Dss.shape[0], Dss.shape[0])))
            if len( self.cond_index ) == 0:
                self.cond_index = np.zeros((self.N_blocks, Dss.shape[0]))
            
            ss_ind = np.where(ss > 0)[0]
            for ii in range(Dss.shape[0]):
                self.cond_index[kk, ii] = ss_ind[ii]      
        
        if not os.path.exists( self.result_file ):
            fh = open(self.result_file, 'w')

            header_str = 'observer'
            for kk in range(len(condition_table.columns)):
                vn = condition_table.columns[kk]
                if vn in block_columns:
                    header_str += ', {}'.format(str(vn))
                else:
                    header_str += ', {}_A, {}_B'.format(str(vn), str(vn))

            header_str += ', is_A_selected\n'
            fh.write(header_str)
            fh.close()
        else:
            # % Load existing results into comparison tables                
            R = pd.read_csv( self.result_file, sep=', ' )
            all_columns = condition_table.columns
            non_block_columns = all_columns[np.logical_not(np.isin(all_columns, block_columns))]
            non_block_columns_A = non_block_columns + '_A'
            non_block_columns_B = non_block_columns + '_B'

            if R.shape[0] == 0:
                return
            
            for kk in range(self.N_blocks): # % for each block

                
                if len( block_columns ) == 0:
                    Rss = R
                    Dss = condition_table
                else:
                    Rss = R[utils.columns_equal(R, block_columns, block_table[kk] ).astype(np.bool)]
                    Dss = condition_table[utils.columns_equal(condition_table, block_columns, block_table[kk] ).astype(np.bool)]
                
                N_missing = 0
                for rr in range(Rss.shape[0]):
                    ind_A = [idx for idx in range(Dss.shape[0]) if np.all(Dss[non_block_columns].to_numpy()[idx] == Rss[non_block_columns_A].iloc[[rr]].to_numpy()[0])]
                    ind_B = [idx for idx in range(Dss.shape[0]) if np.all(Dss[non_block_columns].to_numpy()[idx] == Rss[non_block_columns_B].iloc[[rr]].to_numpy()[0])]
                    
                    if len( ind_A ) == 0 or len( ind_B ) == 0:
                        N_missing = N_missing + 1                            
                        if N_missing < 2:
                            print( 'Results file contains a condition that is missing in the condition_table' )
                            print( Rss.iloc[rr] )
                        continue
                    if Rss['is_A_selected'].iloc[[rr]].to_numpy()[0]:
                        self.cmp_M[kk][ind_A,ind_B] += 1
                    else:
                        self.cmp_M[kk][ind_B,ind_A] += 1
                if N_missing > 0:
                    print( 'Results file contains {} answers that do not match any of the conditions in the condition_table'.format(N_missing) )

    # % Run ASAP to get the next batch of pairs
    def init_pair_list( self ):
        if len(self.pair_list) == 0:
            for kk in range(self.N_blocks):
                pairs = asap_gpu.ASAP( self.cmp_M[kk], mst_mode=True )
                pl_add = np.concatenate([np.ones((pairs.shape[0], 1))*kk, pairs], axis=1)
                if len(self.pair_list) == 0:            
                    self.pair_list = pl_add           
                else:
                    self.pair_list = np.concatenate((self.pair_list, pl_add), axis=0)                      

            self.rand_ord = np.random.permutation( self.pair_list.shape[0] )

    # % Get the number of pairwise comparison left in the current
    # % batch
    def get_pair_left( self ):
        self.init_pair_list()        
        N = self.pair_list.shape[0] - self.pair_index
        return N
    
    def get_next_pair( self ):

        self.pair_index = self.pair_index+1
        if self.pair_index >= self.pair_list.shape[0]:
            self.pair_list = [] # % We need to generate a new set of pairs
            self.pair_index = 1
        
        self.init_pair_list()

        ind = self.rand_ord[self.pair_index]
        block = self.pair_list[ind,0].astype(np.uint8)
        pair = self.pair_list[ind,1:].astype(np.uint8)
        stim_A = self.cond_index[block,pair[0]]
        stim_B = self.cond_index[block,pair[1]]
        
        return stim_A, stim_B
    
    def set_pair_result( self, is_A_selected ):

        ind = self.rand_ord[self.pair_index]
        block = self.pair_list[ind,0].astype(np.uint8)
        pair = self.pair_list[ind,1:].astype(np.uint8)
        stim_A = self.cond_index[block,pair[0]]
        stim_B = self.cond_index[block,pair[1]]
        
        if is_A_selected:
            self.cmp_M[block][pair[0],pair[1]] += 1
        else:
            self.cmp_M[block][pair[1],pair[0]] += 1
        
        fh = open( self.result_file, 'a' )
        write_str = '{}'.format( self.observer )
        for kk in range(len(self.condition_table.columns)):
            vn = self.condition_table.columns[kk]
            if vn in self.block_columns:
                write_str += self.write_value( self.condition_table[vn][stim_A] )
            else:
                write_str += self.write_value( self.condition_table[vn][stim_A] )
                write_str += self.write_value( self.condition_table[vn][stim_B] )

        write_str += ', {:d}\n'.format( is_A_selected )
        fh.write(write_str)
        fh.close()         
    
    def write_value( self, value ):
        write_str = ""

        if type(value) is np.ndarray and len(value)==1:
            value = value[0]
        write_str += ', {}'.format(value)

        return write_str
