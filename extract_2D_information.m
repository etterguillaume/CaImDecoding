function [KL_divergence, PDF, occupancy_vector, prob_being_active, tuning_map ] = extract_2D_information(binarized_trace, interp_behav_vec, ca_time, X_bin_vector, Y_bin_vector, inclusion_vector)
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
interp_behav_vec = interp_behav_vec(inclusion_vector,:);

%% Create bin vectors
numBins = length(X_bin_vector)*length(Y_bin_vector);
prob_being_active = sum(binarized_trace)./length(ca_time);

occupancyProbChance = 1/numBins;
jointProbChance = prob_being_active*occupancyProbChance;

%% Compute joint probabilities (of cell being active while being in a state bin)
tuning_map = zeros(length(Y_bin_vector)-1,length(X_bin_vector)-1);
occupancy_vector = zeros(length(Y_bin_vector)-1,length(X_bin_vector)-1);

    
    for y = 1:length(Y_bin_vector)-1
        for x = 1:length(X_bin_vector)-1
        binarized_spatial_vector = 0*ca_time;
        position_idx = find(interp_behav_vec(:,1)>X_bin_vector(x) & interp_behav_vec(:,1) < X_bin_vector(x+1) & interp_behav_vec(:,2)>Y_bin_vector(y) & interp_behav_vec(:,2) < Y_bin_vector(y+1));
        
        if ~isempty(position_idx)
            binarized_spatial_vector(position_idx)=1;
            occupancy_vector(y,x) = length(position_idx);
            activity_in_bin_idx = find(binarized_trace == 1 & binarized_spatial_vector == 1);
            joint_prob_active_in_bin = length(activity_in_bin_idx)/length(position_idx);
            tuning_map(y,x) = joint_prob_active_in_bin;
        end
        end
    end
    
    occupancy_vector = occupancy_vector./sum(occupancy_vector);

PDF = tuning_map./sum(tuning_map(:));
PDF(PDF==0) = eps; % This is to prevent log values to reach infinity during KL divergence computation

for i = 1:numel(PDF)
KL_temp(i) = PDF(i)*log2(PDF(i)./occupancyProbChance);
end

KL_divergence = sum(KL_temp);
    
    
end



