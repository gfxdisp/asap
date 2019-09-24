function [kl_divs,out_data,arg] = KLD_mat_approx(Nc,init_data)
    betax = 10;
    kl_divs=zeros(Nc);
	init_data.n_iter=4;
    [out_data]=ts_offline(init_data);
    
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
            if init_data.prob_cmps(ii,jj)>rand()
                [msiijj,vsiijj]=ts_online([ii,jj],out_data.Ms,out_data.Vs);
                kl1=norm_kl([msiijj(ii),msiijj(jj)],[out_data.Ms(ii),out_data.Ms(jj)],[vsiijj(ii),vsiijj(jj)],[out_data.Vs(ii),out_data.Vs(jj)]);

                [msjjii,vsjjii]=ts_online([jj,ii],out_data.Ms,out_data.Vs);
                kl2=norm_kl([msjjii(ii),msjjii(jj)],[out_data.Ms(ii),out_data.Ms(jj)],[vsjjii(ii),vsjjii(jj)],[out_data.Vs(ii),out_data.Vs(jj)]);

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
            den = (2*betax^2+out_data.Vs(ii)+out_data.Vs(jj));
            out_data.prob_cmps(ii,jj) = sqrt(2*betax^2/den)*exp(-(out_data.Ms(ii)-out_data.Ms(jj))^2/(2*den));
        end
        out_data.prob_cmps(ii,:) = out_data.prob_cmps(ii,:)./max(out_data.prob_cmps(ii,1:(ii-1)));
    end
    
    arg = arg(randi(size(arg,1)),:);
    
end
