function [kl_divs,out_data,arg] = compute_information_gain_mat(Nc,init_data)

    kl_divs=zeros(Nc);
	init_data.n_iter=4;
    [out_data]=ts_solve(init_data);
    
    out_data.n_iter = 2;
    out_data.Pgs = [out_data.Pgs; zeros(1,2)];
    out_data.Mgs = [out_data.Mgs; zeros(1,2)];
    G = out_data.G;
    arg = [];

    
    maxval = -Inf;
    out_data.prob_cmps = zeros(Nc);
    for ii=2:Nc
        for jj=1:(ii-1)
            pij=normcdf(out_data.Ms(ii)-out_data.Ms(jj), 0, sqrt(1+out_data.Vs(ii)+out_data.Vs(jj)) );
            if pij>rand()
                out_data.G = [G; ii,jj];
                [out_data_iijj]=ts_solve(out_data);
                kl1=kl_divergence_approx(out_data_iijj.Ms,out_data.Ms,out_data_iijj.Vs,out_data.Vs);

                out_data.G = [G; jj,ii];
                [out_data_jjii]=ts_solve(out_data);
                kl2=kl_divergence_approx(out_data_jjii.Ms,out_data.Ms,out_data_jjii.Vs,out_data.Vs);

                kl_gain=pij*kl1+(1-pij)*kl2;
                kl_divs(ii,jj)=kl_gain;
                if kl_gain>maxval
                    maxval=kl_gain;
                    arg = [ii,jj];
                elseif kl_gain==maxval
                    arg = [arg;[ii,jj]];
                else 
                    ;
                end
            else
                kl_divs(ii,jj) = -1;
            end
            
            out_data.prob_cmps(ii,jj) = min(pij,1-pij);
        end
        out_data.prob_cmps(ii,:) = out_data.prob_cmps(ii,:)./max(out_data.prob_cmps(ii,1:(ii-1)));
    end
    arg = arg(randi(size(arg,1)),:);
end