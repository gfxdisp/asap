function [mapping] = id_to_idxes(N)
    mapping = zeros((N-1)*N/2,2);
    count = 1;
    for jj = 2:N
        for ii = 1:(jj-1)

            mapping(count,:) = [ii,jj];
            count = count+1;
        end
    end

end