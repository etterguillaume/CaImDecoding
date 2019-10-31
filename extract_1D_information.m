function [KL_divergence, PDF, occupancy_vector, prob_being_active, tuning_curve ] = extract_2D_information(binarized_trace, interp_behav_vec, ca_time, bin_vector, inclusion_vector)
%MSPLACE_CELL_OF_POSITIONS Analyse and plot spatial information
%   This function plots in space the position related to calcium 

% Copyright (C) 2017-2019 by Guillaume Etter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or any
% later version.
% Contact: etterguillaume@gmail.com

% binary_trace: logical vector representing neurons active/inactive periods
% 
% interp_behav_vec: m x 2  matrix representing behavioral state in two dimensions (eg position in a
% maze). Has to be the same size as binary_trace
%
% ca_time: time vector for both binary_trace and interp_behav_trace. If time correspondance
% is different between these two variables, you can use the fonction
% interp_behav to interpolate the behavioral trace
%
% OPTIONAL vectors:
% inclusion_vec: logical vector including (1) or excluding (0) corresponding timestamps. Has to be the same size as binary_trace.

%% Ignore excluded periods
ca_time = ca_time(inclusion_vector);
binarized_trace = binarized_trace(inclusion_vector);
interp_behav_vec = interp_behav_vec(inclusion_vector);

%% Create bin vectors
numBins = length(bin_vector);
prob_being_active = sum(binarized_trace)./length(ca_time); %Expressed in probability of firing (<1)

occupancyProbChance = 1/numBins;
jointProbChance = prob_being_active*occupancyProbChance;

%% Compute joint probabilities (of cell being active while being in a state bin)
tuning_curve = zeros(length(bin_vector)-1,1);
occupancy_vector = zeros(length(bin_vector)-1,1);

    for i = 1:length(bin_vector)-1
        binarized_spatial_vector = 0*ca_time;
        position_idx = find(interp_behav_vec>bin_vector(i) & interp_behav_vec < bin_vector(i+1));
        
        if ~isempty(position_idx)
            binarized_spatial_vector(position_idx)=1;
            occupancy_vector(i) = length(position_idx);
            activity_in_bin_idx = find(binarized_trace == 1 & binarized_spatial_vector == 1);
            joint_prob_active_in_bin = length(activity_in_bin_idx)/length(position_idx);
            tuning_curve(i) = joint_prob_active_in_bin;
        end
    end
    
    occupancy_vector = occupancy_vector./sum(occupancy_vector);
    
PDF = tuning_curve./sum(tuning_curve);
PDF(PDF==0) = eps; % This is to prevent log values to reach infinity during KL divergence computation

for i = 1:length(PDF)
KL_temp(i) = PDF(i)*log2(PDF(i)./occupancyProbChance);
end

KL_divergence = sum(KL_temp);
    
    
end



