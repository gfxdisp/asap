function [r,c,numb] = calc_numb_cmp(N)

    K=1:N;
    divsrs = K(rem(N,K)==0);
    sdivsrs =size(divsrs,2);
    
    if (rem(sdivsrs,2)==0)
        r = divsrs(sdivsrs/2);
        c = divsrs(sdivsrs/2+1);
    else
        r = divsrs((sdivsrs-1)/2+1);
        c = divsrs((sdivsrs-1)/2+1);
    end

    numb = (r*r-2*r+r*c)*c/2;
end