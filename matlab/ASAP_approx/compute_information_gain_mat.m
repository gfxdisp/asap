function [kl_divs,out_data,arg] = compute_information_gain_mat(Nc,init_data)

    kl_divs=zeros(Nc);
    init_data.n_iter=4;
    
    % compute the scores from the current state of the pairwise comparison matrix
    [out_data]=ts_solve(init_data);
    
    % number of iterations for TrueSkill for the expected information gain calculation (can be set to lower number)
    % as we can use pre-computed messages
    out_data.n_iter = 2;
    
    % use pre-computed messages for message passing to speed up the computation
    out_data.Pgs = [out_data.Pgs; zeros(1,2)];
    out_data.Mgs = [out_data.Mgs; zeros(1,2)];
    G = out_data.G;
    arg = [];

    
    maxval = -Inf;
    out_data.prob_cmps = zeros(Nc);
    
    % Iterate over all possible pairs
    for ii=2:Nc
        for jj=1:(ii-1)
	
	    % pij for weighting the expected information gain
            pij=normcdf(out_data.Ms(ii)-out_data.Ms(jj), 0, sqrt(1+out_data.Vs(ii)+out_data.Vs(jj)) );
	    
	    % Selective EIG evaluation
            if init_data.prob_cmps(ii,jj)>rand()
	    
	        % Solve for the case that ii was selected over jj
                out_data.G = [G; ii,jj];
                [msiijj,vsiijj]=ts_online([ii,jj],out_data.Ms,out_data.Vs);
		
         		% Find kl between the current results from the experiment and those with ii over jj added
                kl1=kl_divergence_approx([msiijj(ii),msiijj(jj)],[out_data.Ms(ii),out_data.Ms(jj)],[vsiijj(ii),vsiijj(jj)],[out_data.Vs(ii),out_data.Vs(jj)]);
       	 	% Solve for the case that jj was selected over ii
                out_data.G = [G; jj,ii];
                [msjjii,vsjjii]=ts_online([jj,ii],out_data.Ms,out_data.Vs);                
                kl2=kl_divergence_approx([msjjii(ii),msjjii(jj)],[out_data.Ms(ii),out_data.Ms(jj)],[vsjjii(ii),vsjjii(jj)],[out_data.Vs(ii),out_data.Vs(jj)]);


                % Expected information gain
                kl_gain=pij*kl1+(1-pij)*kl2;
                kl_divs(ii,jj)=kl_gain;
		
		        % Save the pairs for the maximum attainable information gain
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
            
	    % update the probability of evaluating expected information gain
            out_data.prob_cmps(ii,jj) = min(pij,1-pij);
        end
        out_data.prob_cmps(ii,:) = out_data.prob_cmps(ii,:)./max(out_data.prob_cmps(ii,1:(ii-1)));
    end
    arg = arg(randi(size(arg,1)),:);
end
