%% This example is dedicate to extracting tuning curves
% Used in figure 1

%% Load the data
load 'Linear_track_data/ca_trace'
load 'Linear_track_data/ca_time'
load 'Linear_track_data/behav_vec'
load 'Linear_track_data/behav_time'

%% Binarize calcium trace
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise
[binarized_trace] = extract_binary(ca_trace, sampling_frequency, z_threshold);

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1); % Make sure to interpolate or remove every NaN so that every timepoint has a corresponding behavioral state

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

%% Isolate one specific direction
% In this example, we will look only at when the mouse travels to the right
[direction_indices] = isolate_direction(interp_behav_vec,'right');

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 5 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = direction_indices; % Only include right runs
inclusion_vector(running_ts == 0) = 0; % Only include periods of running

%% Compute occupancy and joint probabilities
bin_size = 3;
% Make sure that your binning vector includes every data point of
% interp_behav_vec using the min/max function:
bin_vector = min(interp_behav_vec):bin_size:max(interp_behav_vec)+bin_size;
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

[MI, PDF, occupancy_vector, prob_being_active, tuning_curve ] = extract_1D_information(binarized_trace, interp_behav_vec, bin_vector, inclusion_vector);
tuning_curve_peak = max(tuning_curve);

figure
subplot(3,1,1)
plot(bin_centers_vector,tuning_curve,'Color', [0 0.1 0.8])
title 'Likelihood'
xlabel 'Location on the track (cm)'
ylabel 'Probability of activity given location'
subplot(3,1,2)
plot(bin_centers_vector,occupancy_vector,'Color', [0.1 0.8 0.1])
title 'Prior (occupancy)'
xlabel 'Location on the track (cm)'
ylabel 'Relative occupancy'
subplot(3,1,3)
plot(bin_centers_vector,PDF,'Color', [0.8 0.2 0])
title 'Posterior probability density function'
xlabel 'Location on the track (cm)'
ylabel 'Probability of being in location given activity'

%% Shuffle data
numShuffles = 1000;
shuffled_tuning_curve = zeros(length(bin_centers_vector), numShuffles);

for k = 1:numShuffles
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(length(binarized_trace),1);

    % Permute the trace
    shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
    shuffled_binarized(random_ts+1:end) = binarized_trace(1:end-random_ts);
    
    % Compute tuning curve
    [~, ~, ~, ~, shuffled_tuning_curve(:,k)] = extract_1D_information(shuffled_binarized, interp_behav_vec, bin_vector, running_ts);
end 

%% Compute significance - Method 1: pN
pvalue = sum(shuffled_tuning_curve > tuning_curve,2)/numShuffles; %  p-value, supra-threshold tests

significant_tuning_curve = tuning_curve;
significant_tuning_curve(pvalue > 0.05) = 0;

%% Compute significance - Method 2: bootstrapping
actual_bootstrap_tuning_curve = zeros(length(bin_centers_vector), numShuffles);
shuffled_bootstrap_tuning_curve = zeros(length(bin_centers_vector), numShuffles);
for k = 1:numShuffles
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(length(binarized_trace),1);
    % Permute the trace
    shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
    shuffled_binarized(random_ts+1:end) = binarized_trace(1:end-random_ts);
    
    bootstrap_ts = inclusion_vector;
    for i = 1:length(bootstrap_ts)
       if bootstrap_ts(i) == 1 & rand < 0.5
          bootstrap_ts(i) = 0;
       end
    end
    
    %% Compute the actual tuning curve using a bootstrapped sample
    [actual_bootstrap_MI(k), actual_bootstrap_PDF(:,k), ~, actual_bootstrap_prob_being_active(k), actual_bootstrap_tuning_curve(:,k) ] = extract_1D_information(binarized_trace, interp_behav_vec, bin_vector, bootstrap_ts);
    
    %% Compute the shuffled tuning curve using the same bootstrapped sample
    [shuffled_bootstrap_MI(k), shuffled_bootstrap_PDF(:,k), ~, shuffled_bootstrap_prob_being_active(k), shuffled_bootstrap_tuning_curve(:,k)] = extract_1D_information(shuffled_binarized, interp_behav_vec, bin_vector, bootstrap_ts);
end

figure
plot(bin_centers_vector,actual_bootstrap_tuning_curve, 'k')
hold on
plot(bin_centers_vector,shuffled_bootstrap_tuning_curve, 'r')

%% Find the 95% confidence interval
sorted_BS_tuning_curves = sort(actual_bootstrap_tuning_curve,2);
CI_idx_loc = 0.95*numShuffles/2;
median_idx = round(numShuffles/2);
upper_CI95_idx = median_idx+CI_idx_loc;
lower_CI95_idx = median_idx-CI_idx_loc;

% This will make sure that upper and lower bounds are withing the actual bootstrap data
upper_CI95_idx(upper_CI95_idx > numShuffles) = numShuffles;
upper_CI95_idx(upper_CI95_idx < 1) = 1;
lower_CI95_idx(lower_CI95_idx > numShuffles) = numShuffles;
lower_CI95_idx(lower_CI95_idx < 1) = 1;

upper_CI95 = sorted_BS_tuning_curves(:,upper_CI95_idx);
lower_CI95 = sorted_BS_tuning_curves(:,lower_CI95_idx);

figure
plot(bin_centers_vector, actual_bootstrap_tuning_curve, 'color', [0.8 0.8 0.8], 'Linewidth', 0.5)
hold on
plot(bin_centers_vector,tuning_curve, 'k', 'Linewidth', 2)
plot(bin_centers_vector,upper_CI95, 'r', 'Linewidth', 2)
plot(bin_centers_vector,lower_CI95, 'r', 'Linewidth', 2)
