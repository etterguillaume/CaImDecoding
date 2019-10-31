function [noisy_tuning_curve_data] = create_noisy_tuning_curves(tuning_curve_data,portion_shuffle)
%CREATE_NOISY_TUNING_CURVES Summary of this function goes here
%   Detailed explanation goes here

numTuningDataPoints = numel(tuning_curve_data);
numTuningDataPointsToShuffle = round(portion_shuffle*numTuningDataPoints);

%numOccupancyPoints = numel(occupancy_vector);
%numOccupancyPointsToShuffle = round(portion_shuffle*numOccupancyPoints);

%numActiveProbPoints = numel(prob_being_active);
%numActiveProbPointsToShuffle = round(portion_shuffle*prob_being_active)

%% TO COMPLETE

max_joint_prob = max(tuning_curve_data(:));

indices_to_shuffle = [];

while length(indices_to_shuffle) < numTuningDataPointsToShuffle
    rand_index = floor(rand*numTuningDataPoints);
    if rand_index ~= 0
        indices_to_shuffle(end+1) = rand_index;
    end
    
end

noisy_tuning_curve_data = tuning_curve_data;

for shuffle_i = 1:numTuningDataPointsToShuffle
   noisy_tuning_curve_data(indices_to_shuffle(shuffle_i)) = rand*max_joint_prob;    
end

end

