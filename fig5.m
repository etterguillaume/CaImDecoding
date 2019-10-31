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

%% Create a behavior vector that discriminates between right and left trajectories
[single_trajectory_vec] = create_single_trajectory_vec(interp_behav_vec);

%% Compute velocities
[ velocity ] = extract_velocity(interp_behav_vec, ca_time);

%% Find periods of immobility
% This will be usefull later to exclude immobility periods from analysis
min_speed_threshold = 5; % 2 cm.s-1
running_ts = velocity > min_speed_threshold;

%% Create an inclusion vector to isolate only a specific behavior
inclusion_vector = running_ts; % Only include periods of running

%% Create a vector for location only
bin_vector_location = 0:3:100; % start : bin_size : end
bin_size = bin_vector_location(2) - bin_vector_location(1);
bin_centers_vector_location = bin_vector_location + bin_size/2;
bin_centers_vector_location(end) = [];

%% Create a vector for location and each direction
bin_vector_location_dir = 0:3:200; % start : bin_size : end
bin_size = bin_vector_location_dir(2) - bin_vector_location_dir(1);
bin_centers_vector_location_dir = bin_vector_location_dir + bin_size/2;
bin_centers_vector_location_dir(end) = [];

binarized_data = zeros(size(ca_data));
for cell_i = 1:size(ca_data,2)
    binarized_data(:,cell_i) = extract_binary(ca_data(:,cell_i), sampling_frequency, z_threshold);
end

%% Neuronal activity reconstruction
% First let's binarize traces from all cells
for trial_i = 1:30
   disp(['Trial ' num2str(trial_i) ' out of 30'])    
    %% Select the frames that are going to be used to train the decoder
    training_set_creation_method = 'random'; % 'odd', odd timestamps; 'first_portion', first portion of the recording; 3, 'random' random frames
    training_set_portion = 0.5; % Portion of the recording used to train the decoder for method 2 and 3
    
    training_ts = create_training_set(ca_time, training_set_creation_method, training_set_portion);
    training_ts(running_ts == 0) = 0; % Exclude periods of immobility from the traing set
    
    decoding_ts = ~training_ts; % Training timestamps are excluded
    decoding_ts(running_ts == 0) = 0; % Periods of immobility are excluded
    
    %% Create tuning curves for every cell
    for cell_i = 1:size(binarized_data,2)
        [KL_divergence_location(cell_i), PDF_location(:,cell_i), occupancy_vector_location, prob_being_active_location(cell_i), tuning_curve_data_location(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), interp_behav_vec, ca_time, bin_vector_location, training_ts);
        [KL_divergence_location_dir(cell_i), PDF_location_dir(:,cell_i), occupancy_vector_location_dir, prob_being_active_location_dir(cell_i), tuning_curve_data_location_dir(:,cell_i) ] = extract_1D_information(binarized_data(:,cell_i), single_trajectory_vec, ca_time, bin_vector_location_dir, training_ts);
    end
    
    %% Actual discretized location for location only
    actual_bin_location = nan*interp_behav_vec;
    actual_position_location = nan*interp_behav_vec;
    for bin_i = 1:length(bin_vector_location)-1
        position_idx = find(interp_behav_vec>bin_vector_location(bin_i) & interp_behav_vec < bin_vector_location(bin_i+1));
        actual_bin_location(position_idx) = bin_i;
        actual_position_location(position_idx) = bin_centers_vector_location(bin_i);
    end
    actual_bin_location(~decoding_ts) = nan;
    actual_position_location(~decoding_ts) = nan;
    actual_bin_location = actual_bin_location';
    actual_position_location =  actual_position_location';
    
    %% Actual discretized location for location and direction
    actual_bin_location_dir = nan*single_trajectory_vec;
    actual_position_location_dir = nan*single_trajectory_vec;
    for bin_i = 1:length(bin_vector_location_dir)-1
        position_idx = find(single_trajectory_vec>bin_vector_location_dir(bin_i) & single_trajectory_vec < bin_vector_location_dir(bin_i+1));
        actual_bin_location_dir(position_idx) = bin_i;
        actual_position_location_dir(position_idx) = bin_centers_vector_location_dir(bin_i);
    end
    actual_bin_location_dir(~decoding_ts) = nan;
    actual_position_location_dir(~decoding_ts) = nan;
    actual_bin_location_dir = actual_bin_location_dir';
    actual_position_location_dir =  actual_position_location_dir';
            
%% Encode information and predict neuronal activity
[estimated_data_location] = encode_neuronal_data(actual_bin_location, occupancy_vector_location, prob_being_active_location, tuning_curve_data_location);
estimated_data_location(estimated_data_location > 0.5) = 1;
estimated_data_location(estimated_data_location <= 0.5) = 0;
[agreement_vector_location] = bayesian_assess_model(binarized_data, estimated_data_location);
mean_agreement_location(trial_i) = mean(agreement_vector_location,'omitnan');

[estimated_data_location_dir] = encode_neuronal_data(actual_bin_location_dir, occupancy_vector_location_dir, prob_being_active_location_dir, tuning_curve_data_location_dir);
estimated_data_location_dir(estimated_data_location_dir > 0.5) = 1;
estimated_data_location_dir(estimated_data_location_dir <= 0.5) = 0;
[agreement_vector_location_dir] = bayesian_assess_model(binarized_data, estimated_data_location_dir);
mean_agreement_location_dir(trial_i) = mean(agreement_vector_location_dir,'omitnan');
end
