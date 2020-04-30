%% This example is dedicate to decoding location on a linear track

%% Load the data
load 'Linear_track_data/ca_data'
load 'Linear_track_data/ca_time'
load 'Linear_track_data/behav_vec'
load 'Linear_track_data/behav_time'

%% Parameters used to binarize the calcium traces
sampling_frequency = 30; % This data set has been sampled at 30 images per second
z_threshold = 2; % A 2 standard-deviation threshold is usually optimal to differentiate calcium ativity from background noise

%% Interpolate behavior
% In most cases, behavior data has to be interpolated to match neural temporal
% activity assuming it has a lower sampling frequency than neural activity
[interp_behav_vec] = interpolate_behavior(behav_vec, behav_time, ca_time);
interp_behav_vec(end) = interp_behav_vec(end-1);
%% Compute velocities
% In this example, we will ignore moments when the mouse is immobile
[velocity] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 5 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Compute occupancy and joint probabilities
bin_size = 3;
% Make sure that your binning vector includes every data point of
% interp_behav_vec using the min/max function:
bin_vector = min(interp_behav_vec):bin_size:max(interp_behav_vec)+bin_size; % start : bin_size : end
bin_centers_vector = bin_vector + bin_size/2;
bin_centers_vector(end) = [];

%% Decoding
% First let's binarize traces from all cells
binarized_data = zeros(size(ca_data));
for cell_i = 1:size(ca_data,2)
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

%% Select the frames that are going to be used to train the decoder
training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.5; % Portion of the recording used to train the decoder for method 2 and 3

training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set

%% Create tuning curves for every cell
for cell_i = 1:size(binarized_data,2)
    [MI(cell_i), PDF(:,cell_i), occupancy_vector, prob_being_active(cell_i), tuning_curve_data(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, training_ts);
end

%% Plot the tunning curves
[~,max_index] = max(tuning_curve_data,[],1);
[~,sorted_index] = sort(max_index);
sorted_tuning_curve_data = tuning_curve_data(:,sorted_index);

figure
imagesc(bin_centers_vector,1:size(ca_data,2),sorted_tuning_curve_data')
daspect([1 1 1])
title 'Neuronal tuning curves'
xlabel 'Position on the track (cm)'
ylabel 'Cell ID'

%% Decode position
% First, let us establish the timestamps used for decoding.
decoding_ts = ~training_ts; % Training timestamps are excluded
decoding_ts(running_ts == 0) = 0; % Periods of immobility are excluded

% Minimal a priori (use to remove experimental a priori)
occupancy_vector = occupancy_vector./occupancy_vector*(1/length(occupancy_vector));

% Establish which cells are going to be used in the decoding process
cell_used = logical(ones(size(ca_data,2),1)); % Let us use every cell for now

[decoded_probabilities] = bayesian_decode1D(binarized_data, occupancy_vector, prob_being_active, tuning_curve_data, cell_used);

%% Let us now estimate the mouse location using the maximum a posteriori (MAP) value
[max_decoded_prob, decoded_bin] = max(decoded_probabilities,[],1);
decoded_position = bin_centers_vector(decoded_bin);

% Before looking at the error rate, we must first bin the actual data using the same bin vector used by
% the decoder
actual_bin = nan*interp_behav_vec;
actual_position = nan*interp_behav_vec;
for bin_i = 1:length(bin_vector)-1
    position_idx = find(interp_behav_vec>bin_vector(bin_i) & interp_behav_vec < bin_vector(bin_i+1));
    actual_bin(position_idx) = bin_i;
    actual_position(position_idx) = bin_centers_vector(bin_i);
end

figure
subplot(2,1,1)
imagesc(ca_time,bin_centers_vector,decoded_probabilities)
title 'Posterior probabilities'
xlabel 'Time (s)'
ylabel 'Position on the track (cm)'
ax1 = gca;
ax1.CLim = [0 0.1];
ax1.XLim = [447 452];
ax1.YDir = 'normal';
subplot(2,1,2)
plot(ca_time,actual_position)
hold on
plot(ca_time, decoded_position)
title 'Actual versus decoded position'
xlabel 'Time (s)'
ylabel 'Location on the track (cm)'
ax2 = gca;
ax2.XLim = [447 452]; % Let's plot a single trajectory
linkaxes([ax1 ax2], 'x')

%% Remove training timestamps to assess decoding error rate
decoded_bin(~decoding_ts) = nan;
decoded_position(~decoding_ts) = nan;
decoded_probabilities(:,~decoding_ts) = nan;

actual_bin(~decoding_ts) = nan;
actual_position(~decoding_ts) = nan;
actual_bin = actual_bin';
actual_position =  actual_position';

%% Compute decoding agreement
decoding_agreement_vector = double(decoded_bin == actual_bin);
decoding_agreement_vector(isnan(decoded_bin)) = nan;
decoding_agreement_vector(isnan(actual_bin)) = nan;
decoding_agreement_vector(isnan(decoding_agreement_vector)) = [];
decoding_agreement = sum(decoding_agreement_vector)./length(decoding_agreement_vector);

%% Compute decoding error
decoding_error = actual_position - decoded_position;
mean_decoding_error = mean(abs(decoding_error), 'omitnan');

%% Compute confusion matrix
confusion_matrix = zeros(length(bin_centers_vector),length(bin_centers_vector));

for actual_i = 1:length(bin_centers_vector)
   for decoded_i = 1:length(bin_centers_vector)
       confusion_matrix(actual_i,decoded_i) = sum(decoded_bin == decoded_i & actual_bin == actual_i)./sum(actual_bin == actual_i);
   end
end

% Plot the confusion matrix
figure
imagesc(bin_centers_vector, bin_centers_vector, confusion_matrix)
title 'Confusion matrix'
xlabel 'Actual position (cm)'
ylabel 'Decoded position (cm)'

