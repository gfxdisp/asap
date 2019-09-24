function [obj,grad, H]=func_alpha(alpha, s, para, pair)
  
    s0=getopt(para, 's0', 0);
    reg_0=getopt(para, 'reg_0', 0);
    reg_alpha=getopt(para, 'reg_alpha', 0);
    reg_s=getopt(para, 'reg_s', 0);
    uni_weight=getopt(para, 'uni_weight', true);
    
    n_anno=length(alpha);
    p=exp(s);
    p0=exp(s0);
    
    obj=-reg_0*(sum(log(p0./(p0+p)))+sum(log(p./(p0+p))));
    grad=zeros(n_anno,1);
    h=zeros(n_anno,1);
       
    for k=1:length(pair)
        if (uni_weight)
            s_k=1;
        else
            s_k=size(pair{k},1);
        end

        alpha_tmp=(alpha(k)*p(pair{k}(:,1))+(1-alpha(k))*p(pair{k}(:,2)));
        diff_tmp=p(pair{k}(:,1))-p(pair{k}(:,2));
        obj=obj-sum(log(alpha_tmp)-log(p(pair{k}(:,1))+p(pair{k}(:,2))))/s_k;        
        grad(k)=-sum(diff_tmp./alpha_tmp)/s_k+reg_alpha*alpha(k);
        h(k)=sum((diff_tmp./alpha_tmp).^2)/s_k+reg_alpha;
    end
    
    obj=obj+reg_s*norm(s,2).^2/2+reg_alpha*norm(alpha,2).^2/2;
    H=sparse(1:n_anno, 1:n_anno, h);
    
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