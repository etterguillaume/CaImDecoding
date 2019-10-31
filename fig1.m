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

%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

%% Isolate one specific direction
% In this example, we will look only at when the mouse travels to the right
[direction_indices] = isolate_direction(interp_behav_vec,'right');

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = direction_indices; % Only include right runs
inclusion_vector(running_ts == 0) = 0; % Only include periods of running

%% Compute occupancy and joint probabilities
bin_vector = 0:3:100; % start : bin_size : end
bin_size = bin_vector(2) - bin_vector(1);
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

[KL_divergence, PDF, occupancy_vector, prob_being_active, tuning_curve ] = extract_1D_information(binarized_trace, interp_behav_vec, ca_time, bin_vector, inclusion_vector);
peak_joint_probability = max(tuning_curve);

figure
subplot(3,1,1)
plot(bin_centers_vector,tuning_curve,'Color', [0 0.1 0.8])
title 'Tuning curve'
xlabel 'Location on the track (cm)'
ylabel 'Joint probability'
subplot(3,1,2)
plot(bin_centers_vector,occupancy_vector,'Color', [0.1 0.8 0.1])
title 'Occupancy'
xlabel 'Location on the track (cm)'
ylabel 'Relative occupancy'
subplot(3,1,3)
plot(bin_centers_vector,PDF,'Color', [0.8 0.2 0])
title 'Probability density function'
xlabel 'Location on the track (cm)'
ylabel 'Probability'

%% Shuffle data
numShuffles = 1000;
shuffled_data = zeros(length(bin_centers_vector), numShuffles);

for k = 1:numShuffles
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(length(binarized_trace),1);
    
    % Permute the trace
    shuffled_binarized(1:random_ts) = binarized_trace(end-random_ts+1:end);
    shuffled_binarized(end-random_ts+1:end) = binarized_trace(1:random_ts);
    
    % Compute tuning curve
    [shuffled_KLD(k), shuffled_PDF(:,k), ~ ,  ~, shuffled_tuning_curve(:,k)] = extract_1D_information(shuffled_binarized, interp_behav_vec, ca_time, bin_vector, running_ts);
    shuffled_peak_joint_probability(k) = max(shuffled_tuning_curve(:,k));    
end

%% Compute z vector
average_shuffled_vector = mean(shuffled_tuning_curve,2);
std_shuffled_vector = std(shuffled_tuning_curve,[],2);

z_vector = (tuning_curve-average_shuffled_vector)./std_shuffled_vector;

figure
plot(bin_centers_vector, z_vector)
hold on
line([bin_centers_vector(1) bin_centers_vector(end)], [2 2], 'Color', 'r', 'LineWidth',5)
title 'Z-scored tuning curve'
xlabel 'Location on the track (cm)'
ylabel 'Z score'