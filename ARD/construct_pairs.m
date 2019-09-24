function pairs=construct_pairs(M)
    added = [];
    pairs = [];
    for ii = 1:size(M,1)
        for jj =1:size(M,2)
            added = [added, M(ii,jj)];
            to_cmp =  [M(ii,:), M(:,jj)'];
            for kk=1:size(to_cmp,2)
                if ~ismember(to_cmp(kk),added)
                    pairs = [pairs;[M(ii,jj),to_cmp(kk)]];
                end
            end
        end
    end

end