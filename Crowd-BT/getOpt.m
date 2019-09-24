function [v] = getOpt(options,opt,default)
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