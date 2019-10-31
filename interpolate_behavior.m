function [interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time)
%INTERP_BEHAV Interpolates behavior data to match neural data
% This function takes behavior data and corresponding timestamps, as well
% as calcium time, and performs linear interpolation so that behavior data
% matches calcium imaging timestamps.

% Copyright (C) 2017-2019 by Guillaume Etter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or any
% later version.
% Contact: etterguillaume@gmail.com

% behav_vec: m x 1  vector or a m x 2 matrix representing behavioral state in one or two dimensions (eg position in a
% maze). m is corresponds to the total number of frames and has to be the
% same size as calcium_time
%
% behav_time: vector containing the timestamps (is s) corresponding to each
% behav frame.

% ca_time: vector containing the timestamps (is s) corresponding to each
% calcium imaging frame. ca_time should be longer or the same length as behav_time

%% Error checking
if size(behav_vec,1) ~= size(behav_time,1)
    error('The length of the behavioral vector is different from the length of its associated timestamps')
end

%% Perform linear interpolation
if size(behav_vec,2) == 1
    interp_behav_vec = interp1(behav_time,behav_vec(:,1),ca_time);
elseif size(behav_vec,2) == 2
    interp_behav_vec(:,1) = interp1(behav_time,behav_vec(:,1),ca_time);
    interp_behav_vec(:,2) = interp1(behav_time,behav_vec(:,2),ca_time);
elseif isempty(behav_vec)
    error('behav_vec is empty.')
else
    error('behav_vec contains too many dimensions. It should be either a m x 1 vector, or a m x 2 matrix')
end

end

