function [out_data,pairs_to_compare] = fill_trueskillmat(Nc,init_data,ii)
    betax = 1;
    
	init_data.n_iter=4;
    [out_data]=ts_offline(init_data);
    
    out_data.Pgs = [out_data.Pgs; zeros(1,2)];
    out_data.Mgs = [out_data.Mgs; zeros(1,2)];
    G = out_data.G;
    pairs_to_compare = [];
    
    maxval = -Inf;
    
    out_data.prob_cmps = zeros(Nc);
    for jj=1:(ii-1)
        den = (2*betax^2+out_data.Vs(ii)+out_data.Vs(jj));
        out_data.prob_cmps(ii,jj) = sqrt(2*betax^2/den)*exp(-(out_data.Ms(ii)-out_data.Ms(jj))^2/(2*den));

        if out_data.prob_cmps(ii,jj)>maxval
            maxval=out_data.prob_cmps(ii,jj);
            pairs_to_compare = [ii,jj];
        elseif out_data.prob_cmps(ii,jj)==maxval
            pairs_to_compare = [pairs_to_compare;[ii,jj]];
        else 
            ;
        end
    end
    pairs_to_compare = pairs_to_compare(randi(size(pairs_to_compare,1)),:);
end