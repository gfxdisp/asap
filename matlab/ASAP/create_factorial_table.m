function T = create_factorial_table( varargin )
% Create a table with all combinations of passed sets
%
% T = create_factorial_table( set1, set2, ... )
% 
% Each set must be a named cell array or a numeric 1D array. The name of
% the passed variable will be used as the name of the column.
%
% Example: 
% set1 = { 'A', 'B' };
% set2 = [1 2];
% create_factorial_table( set1, set2 )


var_types = cell(nargin,1);
var_names = cell(nargin,1);
sz = zeros(nargin,1);
for kk=1:nargin
    var_types{kk} = class( varargin{kk} );
    if strcmp(var_types{kk}, 'cell' )
        var_types{kk} = 'cellstr';
    end      
    var_names{kk} = inputname( kk );
    sz(kk) = length(varargin{kk});
end
N = prod(sz);

T = table( 'Size', [N nargin], 'VariableTypes', var_types, 'VariableNames', var_names );

inds = cell(nargin,1);
for rr=1:N    
    [inds{1:nargin}] = ind2sub( sz, rr );
    for cc=1:nargin
        v = varargin{cc}(inds{cc});
        if isnumeric( v )
            T.(var_names{cc})(rr) = v;
        else
            T.(var_names{cc}){rr} = v{1};
        end
        %T(rr,cc) = varargin{cc}{inds{cc}};
    end
end

end