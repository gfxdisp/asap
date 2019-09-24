% Augments a matrix A by adding a small constant \tau to all its
% zero-valued elements. See the discussion of unique minimizers/Ford's
% condition in the research paper for more information on this.
% Empirically tau=5e-3 works quite well. The paper suggests tau=1e-2.
function A=aug_mat(A,tau)
    A(A==0)=tau;
end