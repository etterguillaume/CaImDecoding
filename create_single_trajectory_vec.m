function [single_trajectory_vec] = create_single_trajectory_vec(interp_behav_vec)
%CREATE_SINGLE_TRAJECTORY_VEC Summary of this function goes here
%   Detailed explanation goes here

z_delta_x = diff(interp_behav_vec);
z_delta_x(end+1) = 0;
z_delta_x(isnan(z_delta_x)) = 0;
z_delta_x = zscore(z_delta_x);

right_trajectories = interp_behav_vec;
right_trajectories(z_delta_x<0.2) = NaN;

left_trajectories = interp_behav_vec;
left_trajectories(z_delta_x>-0.2) = NaN;

left_trajectories = 2*max(interp_behav_vec) - left_trajectories;

single_trajectory_vec = right_trajectories;
single_trajectory_vec(~isnan(left_trajectories)) = left_trajectories(~isnan(left_trajectories));               
                
end

