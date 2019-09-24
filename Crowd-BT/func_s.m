function [obj,grad]=func_s(s, alpha, para, pair)

    s0=getopt(para, 's0', 0);
    reg_0=getopt(para, 'reg_0', 0);
    reg_s=getopt(para, 'reg_s', 0);
    reg_alpha=getopt(para, 'reg_alpha', 0);
    uni_weight=getopt(para, 'uni_weight', true);
    
    p=exp(s);
    p0=exp(s0);
    obj=-reg_0*(sum(log(p0./(p0+p)))+sum(log(p./(p0+p))));
    grad=2*reg_0*(p./(p0+p))-reg_0;

    for k=1:length(pair)
        
        if (uni_weight)
            s_k=1;
        else
            s_k=size(pair{k},1);
        end
        
        obj=obj-sum(log((alpha(k)*p(pair{k}(:,1))+(1-alpha(k))*p(pair{k}(:,2))))-...
            log(p(pair{k}(:,1))+p(pair{k}(:,2))))/s_k;
        for idx=1:size(pair{k},1)
            winner=pair{k}(idx,1);
            loser=pair{k}(idx,2);
            v=(p(winner)/(p(winner)+p(loser))...
                -alpha(k)*p(winner)/(alpha(k)*p(winner)+(1-alpha(k))*p(loser)))/s_k;
            grad(winner)=grad(winner)+v;
            grad(loser)=grad(loser)-v;
        end
    end

    obj=obj+reg_s/2*norm(s,2).^2+reg_alpha/2*norm(alpha,2).^2;
    grad=grad+reg_s.*s;

end

function [v] = getopt(options,opt,default)
if isfield(options,opt)
    if ~isempty(getfield(options,opt))
        v = getfield(options,opt);
    else
        v = default;
    end
else
    v = default;
end
end
