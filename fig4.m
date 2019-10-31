%% Load the open field data
load 'Open_field_data/ca_trace'
load 'Open_field_data/ca_time'
load 'Open_field_data/behav_vec'
load 'Open_field_data/behav_time'

%% Binarize calcium trace
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
[binarized_trace] = extract_binary(ca_trace, sampling_frequency, z_threshold);

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
interp_behav_vec(:,1) = interpolate_behavior(behav_vec(:,1), behav_time, ca_time); % in the X dimension
interp_behav_vec(:,2) = interpolate_behavior(behav_vec(:,2), behav_time, ca_time); % in the Y dimension

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running

%% Compute occupancy and joint probabilities
% You can use 'min(behav_vec)' and 'max(behav_vec)'  to estimate the
% boundaries of the behavior vector (in our case, location in space)
X_bin_vector = 3:3:42; % start : bin_size : end
X_bin_size = X_bin_vector(2) - X_bin_vector(1);
X_bin_centers_vector = X_bin_vector + X_bin_size/2;
X_bin_centers_vector(end) = [];
Y_bin_size = X_bin_size;
Y_bin_vector = X_bin_vector;
Y_bin_centers_vector = X_bin_centers_vector; % this will bin space in a square

[KL_divergence, PDF, occupancy_map, prob_being_active, tuning_map ] = extract_2D_information(binarized_trace, interp_behav_vec, ca_time, X_bin_vector, Y_bin_vector, inclusion_vector);
peak_joint_probability = max(tuning_map(:));

figure
subplot(3,1,1)
imagesc(X_bin_centers_vector,Y_bin_centers_vector,tuning_map)
daspect([1 1 1])
title 'Joint probability'
xlabel 'Position (cm)'
ylabel 'Position (cm)'
colorbar
subplot(3,1,2)
imagesc(X_bin_centers_vector,Y_bin_centers_vector,occupancy_map)
daspect([1 1 1])
title 'Relative occupancy'
xlabel 'Position (cm)'
ylabel 'Position (cm)'
colorbar
subplot(3,1,3)
surf(X_bin_centers_vector,Y_bin_centers_vector,PDF)
daspect([1 1 0.01])
title 'Probability density function'
xlabel 'Position (cm)'
ylabel 'Position (cm)'
colorbar


%% Shuffle data
numShuffles = 1000;
shuffled_data = zeros(length(Y_bin_centers_vector), length(X_bin_centers_vector), numShuffles);
shuffled_peak_joint_probability = [];

for k = 1:numShuffles
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(length(binarized_trace),1);
    
    % Permute the trace
    shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
    shuffled_binarized(end-random_ts+1:end) = binarized_trace(1:random_ts);
    
    [shuffled_KLD(k), shuffled_PDF(:,:,k), ~ , ~, shuffled_data(:,:,k)] = extract_2D_information(shuffled_binarized, interp_behav_vec, ca_time, X_bin_vector, Y_bin_vector, running_ts);
    shuffled_peak_joint_probability(k) = max(shuffled_data(:));
end

%% Compute z map
average_shuffled_map = mean(shuffled_data,3);
std_shuffled_map = std(shuffled_data,[],3);
z_map = (tuning_map-average_shuffled_map)./std_shuffled_map;

figure
imagesc(X_bin_centers_vector,Y_bin_centers_vector,z_map)
daspect([1 1 1])
colorbar
title 'Z-scored tuning map'
xlabel 'Position (cm)'
ylabel 'Position (cm)'
