function [G] = unroll_mat(M)
% turn matrix m into an array with dimensionality (numb_comparisons, 2), where in row kk 
% tupple [3,4] means that in the kth comparison condition 3 was chosen over condition 4
N = size(M,1);
G = [];
    for ii=1:N
        for jj=1:N
            G=[G; repmat([ii jj], M(ii,jj),1)];
            M(ii,jj)=0;
        end
    end
end
