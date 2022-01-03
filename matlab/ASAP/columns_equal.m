function ss = columns_equal( D, columns, values )
% Return boolean vector with 1s for all rows in table D for which all
% values in 'columns' match ;values'. 

if istable( values )
    values = table2cell( values );
end

N = height(D);
if N==0 % In case we have an empty table
    ss = false( 0, 1 );
    return
end

ss = true( N, 1 );
for cc=1:length(columns)
    if isnumeric( D.(columns{cc})(1) )
        ss = ss & (D.(columns{cc}) == values{cc});
    else
        ss = ss & strcmp(D.(columns{cc}), values{cc});
    end
end

end