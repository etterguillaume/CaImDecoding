function [training_ts] = create_training_set(ca_time, training_set_creation_method, training_set_portion)
%CREATE_TRAINING_SET Summary of this function goes here
%   Detailed explanation goes here

training_ts = zeros(1,length(ca_time));

switch training_set_creation_method
    case 'odd'
    odd_vector = 1:length(ca_time);
    training_ts = rem(odd_vector, 2) == 1;
    case 'first_portion'
    first_portion = 1:ceil(length(ca_time*training_set_portion)/2);
    training_ts(first_portion) = 1;
    case 'random'
    training_ts_list = [];
    numTrainingTs = ceil(training_set_portion*length(ca_time));
    random_ts = round(rand*length(ca_time));
    while length(training_ts_list) < numTrainingTs
        while ismember(random_ts,training_ts_list) | random_ts == 0
            random_ts = round(rand*length(ca_time));
        end
        training_ts_list(end+1) = random_ts;
    end
    training_ts = 0*ca_time;
    training_ts(training_ts_list) = 1;
    otherwise
    error("Please define a valid method to create a training set: 'odd', 'first_portion', or 'random'");
end

training_ts = logical(training_ts);


end

