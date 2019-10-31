function [estimated_data] = encode_neuronal_data(actual_bin, occupancy_vector, prob_being_active, tuning_curve_data)
%REVERSE_DECODE Summary of this function goes here
%   Detailed explanation goes here

estimated_data = zeros(size(actual_bin,2), size(tuning_curve_data,2))*nan;

%% Estimate neuron probabilities of being in a given state for each timestep
for step_i = 1:length(actual_bin)  
    if ~isnan(actual_bin(step_i))
    estimated_data(step_i,:) = (tuning_curve_data(actual_bin(step_i),:).*prob_being_active)./occupancy_vector(actual_bin(step_i));
    end
end

end

