function RES = Learn_pair(Number_stimuli, Number_observations, Score, std, error_rate)
for number_pair_required = [1 ceil(Number_stimuli/2)]
    
    trial = 0;
    pcm = zeros(Number_stimuli, Number_stimuli);
    round = 0; % one round represent one complete observation for one observer
    while(trial <=Number_observations)
        
        [pcm] = Active_Learning_BT(pcm, Score, std, number_pair_required, error_rate);
        trial = sum(sum(pcm));
        round = round + 1;
        
        RES{number_pair_required}.BT{round} = Learning_performance_evaluation(pcm, Score, std, 'logit');
        
        RES{number_pair_required}.trial_number{round} = trial;
    end
end