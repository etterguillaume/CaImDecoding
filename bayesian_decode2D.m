function [decoded_probabilities] = bayesian_decode2D(binarized_data, occupancy_vector, prob_being_active, tuning_map_data, cell_used)
%MSTRAIN_BAYESDECODER Summary of this function goes here
%   cell_used: array containing cells used to decode behavior. If 0, include all cells

%% Using only cells specified in cell_used variable
tuning_map_data = tuning_map_data(:,:,cell_used);
prob_being_active = prob_being_active(cell_used);
binarized_data = binarized_data(:,cell_used);

decoded_probabilities = zeros(size(tuning_map_data,1),size(tuning_map_data,2),size(binarized_data,1))*nan;

%% Estimate probabilities of being in a given state for each timestep
for step_i = 1:size(binarized_data,1)    
    bayesian_step_prob = [];
    for cell_i = 1:size(binarized_data,2)
        if binarized_data(step_i,cell_i) == 1
            active_tuning_map = tuning_map_data(:,:,cell_i);
            bayesian_step_prob(:,:,cell_i) = (active_tuning_map.*occupancy_vector)./prob_being_active(cell_i);
        elseif binarized_data(step_i,cell_i) == 0
            inactive_tuning_map = 1-tuning_map_data(:,:,cell_i);
            bayesian_step_prob(:,:,cell_i) = (inactive_tuning_map.*occupancy_vector)./(1-prob_being_active(cell_i));
        end
    end
    
    temp_decoded_probabilities = expm1(sum(log1p(bayesian_step_prob),3)); % This should be used instead of simple product to avoid numerical underflow
    decoded_probabilities(:,:,step_i) = temp_decoded_probabilities./sum(temp_decoded_probabilities(:),'omitnan');       
end

end

