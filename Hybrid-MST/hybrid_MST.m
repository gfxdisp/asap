function [pairs_to_compare] = hybrid_MST(M)
%% input
% M - matrix with comparisons collected so far

%% output
% pairs_to_compare - pairs to compare

    number_conditions = size(M,1);

    M = M+Initial_learning(M);


    [b_prior,stats] = paired_comparisons(M,'logit');

    covb = stats.covb; 
    mu_v = b_prior;
    n=30;
    eig = zeros(number_conditions, number_conditions);

    for i = 1:number_conditions-1
        for j = i+1:number_conditions
           mu(i,j) = mu_v(i)-mu_v(j);
           sigma(i,j) = sqrt(covb(i,i)+covb(j,j)-2.*covb(i,j));
           eig(i,j) = Gaussian_Hermite_BT(mu(i,j), sigma(i,j), n);
           eig(j,i) = eig(i,j);
        end
    end

    eig_inverse = 1./eig;
    eig_inverse(find(eig_inverse==Inf))= 0;
    [pairs_to_compare, ~] = prim(eig_inverse);
end
