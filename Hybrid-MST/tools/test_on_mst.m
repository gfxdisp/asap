function pcm = test_on_mst(Number_stimuli, Number_observations, q, M)
%% test on mst model
    trial = sum(sum(M));
    trial_i = sum(sum(M));
    pcm = M;%zeros(Number_stimuli, Number_stimuli);
    full_pair_num = Number_stimuli*(Number_stimuli-1)/2;
    
    while(trial<(Number_observations))
        if sum(sum(pcm)) < full_pair_num % initialization will occupy 1 fpc number
            [pcm] = Active_Learning_MST(pcm, q,Number_observations);
            %[pcm] = Active_Learning_BT(pcm, q, 1);
        else
            [pcm] = Active_Learning_MST(pcm, q,Number_observations);
        end
        if trial_i==0
        trial = sum(sum(pcm))-full_pair_num; % due to initialization 
        else
            trial = sum(pcm(:));
        end
    end
    
end