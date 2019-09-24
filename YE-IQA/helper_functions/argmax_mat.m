% This function returns the index [a,b] that maximizes a matrix M. It also
% optionally returns M(a,b). Note that the function assumes a unique
% maximum.
function [arg, maxval]=argmax_mat(M)
    [I,J]=size(M);
    maxval=-Inf;
    for i=1:I
        for j=1:J
            if M(i,j)>maxval
                arg=[i,j];
                maxval=M(i,j);
            end
        end
    end
end