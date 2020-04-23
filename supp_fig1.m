load 'Linear_track_data/ca_data'
load 'Linear_track_data/ca_time'
load 'Linear_track_data/behav_vec'
load 'Linear_track_data/behav_time'

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
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;
[direction_indices] = isolate_direction(interp_behav_vec,'right');
inclusion_vector = direction_indices; % Only include right runs
inclusion_vector(running_ts == 0) = 0; % Only include periods of running
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

training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
training_set_portion = 0.5; % Portion of the recording used to train the decoder for method 2 and 3

for trial_i = 1:30
    bootstrap_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
    bootstrap_ts = bootstrap_ts == 1 & inclusion_vector' == 1; % Exclude periods of immobility from the traing set
    
    random_ts = ceil(rand*length(ca_time));
    shuffled_binarized = zeros(size(binarized_data));
    
    % Permute the trace
    shuffled_binarized(1:random_ts,:) = binarized_data(end-random_ts+1:end,:);
    shuffled_binarized(random_ts+1:end,:) = binarized_data(1:end-random_ts,:);
    
    for cell_i = 1:size(binarized_data,2)
        [MI(trial_i,cell_i), PDF(trial_i,:,cell_i), occupancy_vector, prob_being_active(trial_i,cell_i), tuning_curve_data(trial_i,:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, bin_vector, bootstrap_ts);
        [shuffled_MI(trial_i,cell_i), shuffled_PDF(trial_i,:,cell_i), ~, ~, shuffled_tuning_curve(trial_i,:,cell_i)] = extract_1D_information(shuffled_binarized(:,cell_i), interp_behav_vec, bin_vector,  bootstrap_ts);
    end
end

%% Let us now sort cells by their mean MI
[~, sorted_index] = sort(mean(MI),'descend');

figure
plot(MI(:,sorted_index)','color', [0.8 0 0]);
hold on
plot(shuffled_MI(:,sorted_index)', 'color', [0.4 0.4 0.4]);
xlabel 'Cell ID'
ylabel 'MI (bits)'




