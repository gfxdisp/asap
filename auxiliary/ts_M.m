function [Ms, Vs] = ts_M(M)
    G=[];
    N=size(M,1);
    for ii=1:N
        for jj=1:N
            G=[G; repmat([ii jj], M(ii,jj),1)];
            M(ii,jj)=0;
        end
    end
    [Ms,Vs]=ts_raw(N,G);
end