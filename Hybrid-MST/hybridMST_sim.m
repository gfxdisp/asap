function [mstres] = hybridMST_sim(q,totalnum, M )
    warning('off','stats:glmfit:IllConditioned')
    N =  size(M,1);
    mstres = test_on_mst(N, totalnum,q,M);
    if sum(sum(M)) == 0
        mstres = mstres-0.5.*ones(size(mstres)) + diag(diag(0.5.*ones(size(mstres))));
    end
    warning('on','stats:glmfit:IllConditioned')
end