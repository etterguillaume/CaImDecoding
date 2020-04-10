function [MI, posterior, occupancy_vector, prob_being_active, likelihood ] = extract_1D_information(binarized_trace, interp_behav_vec, bin_vector, inclusion_vector)
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
binarized_trace = binarized_trace(inclusion_vector);
interp_behav_vec = interp_behav_vec(inclusion_vector);

%% Create bin vectors
prob_being_active = sum(binarized_trace)./length(binarized_trace); % Expressed in probability of firing (<1)

%% Compute joint probabilities (of cell being active while being in a state bin)
likelihood = zeros(length(bin_vector)-1,1);
occupancy_vector = zeros(length(bin_vector)-1,1);
MI = 0;

for i = 1:length(bin_vector)-1
    binarized_spatial_vector = 0*binarized_trace;
    position_idx = find(interp_behav_vec >= bin_vector(i) & interp_behav_vec < bin_vector(i+1));
    
    if ~isempty(position_idx)
        binarized_spatial_vector(position_idx)=1;
        occupancy_vector(i) = length(position_idx)/length(binarized_trace);
        activity_in_bin_idx = find(binarized_trace == 1 & binarized_spatial_vector == 1);
        inactivity_in_bin_idx = find(binarized_trace == 0 & binarized_spatial_vector == 1);
        likelihood(i) = length(activity_in_bin_idx)/length(position_idx);
        
        joint_prob_active = length(activity_in_bin_idx)./length(binarized_trace);
        joint_prob_inactive = length(inactivity_in_bin_idx)./length(binarized_trace);
        prob_in_bin = length(position_idx)./length(binarized_trace);
        
        if joint_prob_active ~= 0
            MI = MI + joint_prob_active*log2(joint_prob_active./(prob_in_bin*prob_being_active));
        end
        if joint_prob_inactive ~= 0
            MI = MI + joint_prob_inactive*log2(joint_prob_inactive./(prob_in_bin*(1-prob_being_active)));
        end
    end
end

posterior = likelihood.*occupancy_vector/prob_being_active;

end



