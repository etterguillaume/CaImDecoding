%% Bayesian reconstruction of neuronal activity
%% Load the linear track data
load 'Linear_track_data/ca_data'
load 'Linear_track_data/ca_time'
load 'Linear_track_data/behav_vec'
load 'Linear_track_data/behav_time'

%% Binarization parameters
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise

%% Interpolate behavior
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1);
%% Create a behavior vector that discriminates between right and left trajectories
[right_indices] = isolate_direction(interp_behav_vec,'right');
[left_indices] = isolate_direction(interp_behav_vec,'left');

%% Compute velocities
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 5 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
both_inclusion_vector = running_ts; % Only include periods of running

right_inclusion_vector = right_indices; % Only include right runs
right_inclusion_vector(running_ts == 0) = 0; % Only include periods of running

left_inclusion_vector = left_indices; % Only include right runs
left_inclusion_vector(running_ts == 0) = 0; % Only include periods of running

%% Create a vector for location only
bin_size = 3;
% Make sure that your binning vector includes every data point of
% interp_behav_vec using the min/max function:
bin_vector = min(interp_behav_vec):bin_size:max(interp_behav_vec)+bin_size;
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

binarized_data = zeros(size(ca_data));
for cell_i = 1:size(ca_data,2)
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

%% Let us compute the tuning curves and corresponding MI values considering either left, right or both directions
numShuffles = 1000;
for k = 1:numShuffles
    random_ts = create_training_set(ca_time,'random',0.5);
    
    both_bs_ts = both_inclusion_vector;
    both_bs_ts = both_bs_ts == 1 & random_ts' == 1;
    
    right_bs_ts = right_inclusion_vector;
    right_bs_ts = right_bs_ts == 1 & random_ts' == 1;
    
    left_bs_ts = left_inclusion_vector;
    left_bs_ts = left_bs_ts == 1 & random_ts' == 1;
    
    %% Compute the actual tuning curve using a bootstrapped sample
    [both_MI(k), ~, ~, ~, both_tuning_curve(:,k)] = extract_1D_information(binarized_data(:,4), interp_behav_vec, bin_vector, both_bs_ts);
    [right_MI(k), ~, ~, ~, right_tuning_curve(:,k)] = extract_1D_information(binarized_data(:,4), interp_behav_vec, bin_vector, right_bs_ts);
    [left_MI(k), ~, ~, ~, left_tuning_curve(:,k)] = extract_1D_information(binarized_data(:,4), interp_behav_vec, bin_vector, left_bs_ts);
end

%% Find the 95% confidence interval
CI_idx_loc = 0.95*numShuffles/2;
median_idx = round(numShuffles/2);
upper_CI95_idx = median_idx+CI_idx_loc;
lower_CI95_idx = median_idx-CI_idx_loc;

sorted_both_tuning_curves = sort(both_tuning_curve,2);
sorted_right_tuning_curves = sort(right_tuning_curve,2);
sorted_left_tuning_curves = sort(left_tuning_curve,2);

% This will make sure that upper and lower bounds are withing the actual bootstrap data
upper_CI95_idx(upper_CI95_idx > numShuffles) = numShuffles;
upper_CI95_idx(upper_CI95_idx < 1) = 1;
lower_CI95_idx(lower_CI95_idx > numShuffles) = numShuffles;
lower_CI95_idx(lower_CI95_idx < 1) = 1;

both_upper_CI95 = sorted_both_tuning_curves(:,upper_CI95_idx);
both_lower_CI95 = sorted_both_tuning_curves(:,lower_CI95_idx);

right_upper_CI95 = sorted_right_tuning_curves(:,upper_CI95_idx);
right_lower_CI95 = sorted_right_tuning_curves(:,lower_CI95_idx);

left_upper_CI95 = sorted_left_tuning_curves(:,upper_CI95_idx);
left_lower_CI95 = sorted_left_tuning_curves(:,lower_CI95_idx);

figure
plot(bin_centers_vector, both_tuning_curve, 'color', [0.8 0 0.8], 'Linewidth', 0.3)
hold on
plot(bin_centers_vector, right_tuning_curve, 'color', [0.8 0 0], 'Linewidth', 0.3)
plot(bin_centers_vector, left_tuning_curve, 'color', [0 0 0.8], 'Linewidth', 0.3)

plot(bin_centers_vector, mean(both_tuning_curve,2), 'color', [0.4 0 0.4], 'Linewidth', 2)
plot(bin_centers_vector,both_upper_CI95, 'color', [0.4 0 0.4], 'Linewidth', 1)
plot(bin_centers_vector,both_lower_CI95, 'color', [0.4 0 0.4], 'Linewidth', 1)

plot(bin_centers_vector, mean(right_tuning_curve,2), 'color', [0.4 0 0], 'Linewidth', 2)
plot(bin_centers_vector,right_upper_CI95, 'color', [0.4 0 0], 'Linewidth', 1)
plot(bin_centers_vector,right_lower_CI95, 'color', [0.4 0 0], 'Linewidth', 1)

plot(bin_centers_vector, mean(left_tuning_curve,2), 'color', [0 0 0.4], 'Linewidth', 2)
plot(bin_centers_vector,left_upper_CI95, 'color', [0 0 0.4], 'Linewidth', 1)
plot(bin_centers_vector,left_lower_CI95, 'color', [0 0 0.4], 'Linewidth', 1)
