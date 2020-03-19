function [G] = unroll_mat(M)
N = size(M,1);
G = [];
    for ii=1:N
        for jj=1:N
            G=[G; repmat([ii jj], M(ii,jj),1)];
            M(ii,jj)=0;
        end
    end
end