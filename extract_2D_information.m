function [MI, posterior, occupancy_map, prob_being_active, prior ] = extract_2D_information(binarized_trace, interp_behav_vec, X_bin_vector, Y_bin_vector, inclusion_vector)
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
interp_behav_vec = interp_behav_vec(inclusion_vector,:);

%% Create bin vectors
prob_being_active = sum(binarized_trace)./length(binarized_trace);

%% Compute joint probabilities (of cell being active while being in a state bin)
prior = zeros(length(Y_bin_vector)-1,length(X_bin_vector)-1);
occupancy_map = zeros(length(Y_bin_vector)-1,length(X_bin_vector)-1);
MI = 0;

for y = 1:length(Y_bin_vector)-1
    for x = 1:length(X_bin_vector)-1
        binarized_spatial_vector = 0*binarized_trace;
        position_idx = find(interp_behav_vec(:,1) >= X_bin_vector(x) & interp_behav_vec(:,1) < X_bin_vector(x+1) & interp_behav_vec(:,2) >= Y_bin_vector(y) & interp_behav_vec(:,2) < Y_bin_vector(y+1));
        
        if ~isempty(position_idx)
            binarized_spatial_vector(position_idx)=1;
            occupancy_map(y,x) = length(position_idx)/length(binarized_trace);
            activity_in_bin_idx = find(binarized_trace == 1 & binarized_spatial_vector == 1);
            inactivity_in_bin_idx = find(binarized_trace == 0 & binarized_spatial_vector == 1);
            prior(y,x) = length(activity_in_bin_idx)/length(position_idx);
            
            joint_prob_active = length(activity_in_bin_idx)/length(binarized_trace);
            joint_prob_inactive = length(inactivity_in_bin_idx)/length(binarized_trace);
            prob_in_bin = length(position_idx)./length(binarized_trace);
            
            if joint_prob_active ~= 0
                MI = MI + joint_prob_active*log2(joint_prob_active./(prob_in_bin*prob_being_active));
            end
            if joint_prob_inactive ~= 0
                MI = MI + joint_prob_inactive*log2(joint_prob_inactive./(prob_in_bin*(1-prob_being_active)));
            end
        end
    end
end

posterior = prior.*occupancy_map/prob_being_active;

end



