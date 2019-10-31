function [decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, cell_used)
%MSTRAIN_BAYESDECODER Summary of this function goes here
%   cell_used: array containing cells used to decode behavior. If 0, include all cells

%% Using only cells specified in cell_used variable
tuning_curve_data = tuning_curve_data(:,cell_used);
prob_being_active = prob_being_active(cell_used);
binarized_data = binarized_data(:,cell_used);

decoded_probabilities = zeros(size(tuning_curve_data,1),size(binarized_data,1))*nan;

%% Estimate probabilities of being in a given state for each timestep
for step_i = 1:size(binarized_data,1)    
    bayesian_step_prob = [];
    for cell_i = 1:size(binarized_data,2)
        if binarized_data(step_i,cell_i) == 1
            active_tuning_curve = tuning_curve_data(:,cell_i);
            bayesian_step_prob(:,cell_i) = (active_tuning_curve.*occupancy_vector)./prob_being_active(cell_i);
        elseif binarized_data(step_i,cell_i) == 0
            inactive_tuning_curve = 1-tuning_curve_data(:,cell_i);
            bayesian_step_prob(:,cell_i) = (inactive_tuning_curve.*occupancy_vector)./(1-prob_being_active(cell_i));
        end
    end
    
    decoded_probabilities(:,step_i) = expm1(sum(log1p(bayesian_step_prob),2)); % This should be used instead of simple product to avoid numerical underflow
    decoded_probabilities(:,step_i) = decoded_probabilities(:,step_i)./sum(decoded_probabilities(:,step_i),'omitnan');    
    
end

end

