classdef PwcmpASAPScheduler
% Class for scheduling pairwise comparisons using ASAP algorithm
%
% This class handles loading and saving results between experiments runs
% and also splits the comparison into blocks. Each block of measurements is
% done idependently with no comparisons made between the blocks. 
%
% Details on ASAP: 
%  https://github.com/gfxdisp/asap
%  

properties
    
    result_file;
    observer;
    condition_table;
    block_columns;
    block_same_conditions;
    
    cmp_M = []; % 
    cmp_inf = []; %
    N_blocks;
    
    cond_index = []; % [block, cond] matrix with the indices of each (block,cond) pair in condition_table    
    
    % List of condition pairs to compare, 1st column - block index, 2nd and
    % 3d column is the pair of indices of conditions WITHIN the black (not global)
    pair_list = []; 
    rand_ord = [];
    pair_index = 0; % next pair to compare
end

    methods
        
        % Creates new ASAP scheduler and initializes with the mesurements
        % collected so far. 
        %
        % result file - a CVS file with the measurememts. It will create a
        %               new file if the file does not exist
        % observer - a string with the ID of the participant
        % condition table - a populated table data structure with one or
        %            more columns that describe the dimensions of measured conditions.
        %            For example: a table with columns (scene, method, parameter_A, parameter_B)
        % block_columns - a cell array with the names of the columns that define separate
        %            blocks of pairwise comparison measurememts. Each block
        %            is measured independently, i.e. no comparisons are
        %            made between the conditions in each block. If no
        %            blocks should be created (all conditions should be
        %            compared with each other), pass an empty cell {} or skip this
        %            parameter. 
        %            Example: { 'scene' } -
        %            the comparisons will be made only between the
        %            conditions for which the 'scene' attribute is the same.         
        function sch = PwcmpASAPScheduler( result_file, observer, condition_table, block_columns, block_same_conditions )
        
            if ~exist( 'run_asap', 'var' )
                addpath( fullfile( pwd, '..', '..', 'asap', 'matlab', 'ASAP' ) );
            end
            
            if ~exist( 'block_columns', 'var' )
                block_columns = {};
            end

            if ~exist( 'block_same_conditions', 'var' )
                block_same_conditions = {};
            end
            
            sch.observer = observer;
            sch.result_file = result_file;
            sch.condition_table = condition_table;
            sch.block_columns = block_columns;
            sch.block_same_conditions = block_same_conditions; 
            
            if isempty( block_columns )
                block_table = [];
                sch.N_blocks = 1;
            else
                block_table = unique( condition_table(:,block_columns), 'rows' );
                sch.N_blocks = height(block_table);
            end
            
            % Populate cond_index so that we can quickly find a condition
            % based on its index value
            sch.cmp_M = cell( sch.N_blocks, 1 );
            sch.cmp_inf = cell( sch.N_blocks, 1); 

            for kk=1:sch.N_blocks
                if isempty( block_columns)
                    ss = true(height(condition_table),1);
                    Dss = condition_table;
                else
                    ss = columns_equal(condition_table, block_columns, block_table(kk,:) );
                    Dss = condition_table( ss, : );
                end
                sch.cmp_M{kk} = zeros(height(Dss));
                sch.cmp_inf{kk} = ones(height(Dss));

                % Define the edges that you do not want to compare between 

                if ~isempty(sch.block_same_conditions)
                    for bb = 1:length(sch.block_same_conditions)
                        for ii = 1:height(Dss)
                            for jj = 1:height(Dss)
                                if strcmp(Dss.(sch.block_same_conditions{bb}){ii}, Dss.(sch.block_same_conditions{bb}){jj})
                                    sch.cmp_inf{kk}(ii, jj) = 0 ; 
                                end 
                            end 
                        end
                    end 
                end 

                if isempty( sch.cond_index )
                    sch.cond_index = nan(sch.N_blocks,height(Dss));
                end
                ss_ind = find(ss);
                for ii=1:height(Dss)
                    sch.cond_index(kk,ii) = ss_ind(ii);
                end
            end            
            
            if ~isfile( sch.result_file )
                fh = fopen( sch.result_file, 'wb' );
                fprintf( fh, 'observer' );
                for kk=1:size(condition_table,2)
                    vn = condition_table.Properties.VariableNames{kk};
                    if ismember( vn, block_columns )
                        fprintf( fh, ', %s', vn );
                    else
                        fprintf( fh, ', %s_A, %s_B', vn, vn );
                    end
                end
                fprintf( fh, ', is_A_selected\n' );
                fclose( fh );
            else
                % Load existing results into comparison tables                
                R = readtable( sch.result_file );
                all_columns = condition_table.Properties.VariableNames;
                non_block_columns = all_columns(~ismember(all_columns,block_columns));
                non_block_columns_A = strcat( non_block_columns, '_A' );
                non_block_columns_B = strcat( non_block_columns, '_B' );
                for kk=1:sch.N_blocks % for each block
                    if isempty( block_columns )
                        Rss = R;
                        Dss = condition_table;
                    else
                        Rss = R(columns_equal(R, block_columns, block_table(kk,:) ), : );
                        Dss = condition_table( columns_equal(condition_table, block_columns, block_table(kk,:) ), : );
                    end
                    
                    N_missing = 0;
                    for rr=1:height(Rss)
                        ind_A = find( columns_equal( Dss, non_block_columns, Rss(rr,non_block_columns_A) ) );
                        ind_B = find( columns_equal( Dss, non_block_columns, Rss(rr,non_block_columns_B) ) );
                        if isempty( ind_A ) || isempty( ind_B )
                            N_missing = N_missing + 1;                            
                            if N_missing<2
                                warning( 'Results file contains a condition that is missing in the condition_table' );
                                display( Rss(rr,:) );
                            end
                            continue;
                        end
                        if Rss.is_A_selected(rr)
                            sch.cmp_M{kk}(ind_A,ind_B) = sch.cmp_M{kk}(ind_A,ind_B)+1;
                        else
                            sch.cmp_M{kk}(ind_B,ind_A) = sch.cmp_M{kk}(ind_B,ind_A)+1;
                        end
                    end
                    if N_missing>0
                        warning( 'Results file contains %d answers that do not match any of the conditions in the condition_table', N_missing );
                    end
                end                
                
            end
            
            
        end
        
        % Run ASAP to get the next batch of pairs
        function sch = init_pair_list( sch )
            if isempty( sch.pair_list )
                for kk=1:sch.N_blocks
                    pairs = run_asap( sch.cmp_M{kk}, sch.cmp_inf{kk}, 'mst');                    
                    pl_add = [ones(size(pairs,1),1)*kk pairs];                    
                    sch.pair_list = cat( 1, sch.pair_list, pl_add );                                        
                end            
                sch.rand_ord = randperm( size(sch.pair_list,1) );
            end
        end

        % Get the number of pairwise comparison left in the current
        % batch
        function [sch, N] = get_pair_left( sch )
            sch = sch.init_pair_list();            
            N = size(sch.pair_list,1) - sch.pair_index;
        end
        
        function [sch, stim_A, stim_B] = get_next_pair( sch )

            sch.pair_index = sch.pair_index+1;
            if sch.pair_index > size(sch.pair_list,1)
                sch.pair_list = []; % We need to generate a new set of pairs
                sch.pair_index = 1;
            end
            
            sch = sch.init_pair_list();

            ind = sch.rand_ord(sch.pair_index);
            block = sch.pair_list(ind,1);
            pair = sch.pair_list(ind,2:3);
            stim_A = sch.cond_index(block,pair(1));
            stim_B = sch.cond_index(block,pair(2));
            
        end       
        
        function [sch] = set_pair_result( sch, is_A_selected )

            ind = sch.rand_ord(sch.pair_index);
            block = sch.pair_list(ind,1);
            pair = sch.pair_list(ind,2:3);
            stim_A = sch.cond_index(block,pair(1));
            stim_B = sch.cond_index(block,pair(2));
            
            if is_A_selected
                sch.cmp_M{block}(pair(1),pair(2)) = sch.cmp_M{block}(pair(1),pair(2))+1;
            else
                sch.cmp_M{block}(pair(2),pair(1)) = sch.cmp_M{block}(pair(2),pair(1))+1;
            end
            
            fh = fopen( sch.result_file, 'a+' );
            fprintf( fh, '"%s"', sch.observer );
            for kk=1:size(sch.condition_table,2)
                vn = sch.condition_table.Properties.VariableNames{kk};
                if ismember( vn, sch.block_columns )
                    sch.write_value( fh, sch.condition_table.(vn)(stim_A) );
                else
                    sch.write_value( fh, sch.condition_table.(vn)(stim_A) );
                    sch.write_value( fh, sch.condition_table.(vn)(stim_B) );
                end
            end
            fprintf( fh, ', %d\n', is_A_selected );
            fclose( fh );            
            
        end
        
    end
    
    methods(Static)
        
        function write_value( fh, value )
            if iscell(value) && numel(value)==1
                value = value{1};
            end
            if isnumeric( value )
                if isinteger( value )
                    fprintf( fh, ', %d', value );
                else
                    fprintf( fh, ', %f', value );
                end
            elseif ischar( value )
                    fprintf( fh, ', "%s"', value );
            else
                error( 'Unsupported' );
            end
        end
    end

end