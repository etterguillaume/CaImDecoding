function [ velocity ] = extract_velocity(interp_behav_vec, ca_time)
%GE_VELOCITY Summary of this function goes here
%   Computes velocity from X,Y coordinates, and corresponding time vector.
% Output units depend on x,y units (eg mm) and time units (should be
% seconds). Results are usually expressed in mm/s. Always check if other
% functions convert to cm/s for ease of interpretation (most common in
% rodents studies)

velocity(1,:) = 0; % Velocity at t1 (coordinates x1,y1) to be zero. Delta_v1 is estimated for a displacement to coordinates x2,y2.

%% Extract velocity
if size(interp_behav_vec,2) == 1
    for i=2:length(ca_time)
        velocity(i,:) = abs(interp_behav_vec(i,1)-interp_behav_vec(i-1,1))/(ca_time(i)-ca_time(i-1));
    end
elseif size(interp_behav_vec,2) == 2
    for i=2:length(ca_time)
        velocity(i,:) = sqrt((interp_behav_vec(i,1)-interp_behav_vec(i-1,1)).^2 + (interp_behav_vec(i,2)-interp_behav_vec(i-1,2)).^2)/(ca_time(i)-ca_time(i-1));
    end
elseif isempty(interp_behav_vec)
    error('interp_behav_vec is empty.')
else
    error('interp_behav_vec contains too many dimensions. It should be either a m x 1 vector, or a m x 2 matrix')
end

%% Smoothing
% This helps reduce big jumps in velocities
velocity = smooth(velocity,round(1/mode(diff(ca_time))));

end