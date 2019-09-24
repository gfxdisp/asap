function [M] = remove_excess_and_simulate(design, M, q_mat, n_comp)
    N = size(M,1);
    numb_max = 0;
    if design<=2
        if design ==2 
            numb_max =1;
        end
        v_comp_perform = M(:);
        while sum(v_comp_perform) > n_comp
            % from which can we delete a comparison?
            feasible = v_comp_perform(v_comp_perform>numb_max);
            if numel(feasible)==0
                break
            end
            % remove the comparison
            ind_deleted = randi(numel(feasible));

            feasible(ind_deleted) = feasible(ind_deleted) - 1;
            v_comp_perform(v_comp_perform>numb_max) = feasible;
            
        end
        comp_perform = reshape(v_comp_perform,N,N);

        comp_perform = comp_perform + comp_perform';
        M = zeros(N,N);

        for ii=1:N
            for jj=(ii+1):N
                for count=1:comp_perform(ii,jj)
                    M = simulate_observer_choice(q_mat,ii,jj,M);
                end
            end
        end
    end
end