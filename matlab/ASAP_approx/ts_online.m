function [mv, pv] = ts_online(G, mv, pv)

    N = size(G,1);
    
    for ii = 1:N
        beta = 1;
        c = sqrt(2*beta+pv(G(ii,1))+pv(G(ii,2)));
        term = (mv(G(ii,1))-mv(G(ii,2)))/c;
        normpdf_term = normpdf(term);
        normcdf_term = normcdf(term);
        fact_1 =v(normpdf_term,normcdf_term);
        fact_2 =w(normpdf_term,normcdf_term,term);
        mv(G(ii,1)) = mv(G(ii,1))+pv(G(ii,1))/c*fact_1;
        mv(G(ii,2)) = mv(G(ii,2))-pv(G(ii,2))/c*fact_1;
        pv(G(ii,1)) = pv(G(ii,1))*(1-pv(G(ii,1))*fact_2/c^2);
        pv(G(ii,2)) = pv(G(ii,2))*(1-pv(G(ii,2))*fact_2/c^2);
    end

    function [res] = v(pdft,cdft)
        res = pdft./cdft;
    end
    function [res] = w(pdft,cdft,x)
        res = v(pdft,cdft).*(v(pdft,cdft) + x);
    end
end
