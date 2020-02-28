function [training_ts] = create_training_set(ca_time, training_set_creation_method, training_set_portion)
%CREATE_TRAINING_SET Summary of this function goes here
%   Detailed explanation goes here

training_ts = zeros(1,length(ca_time));
numFrames = length(ca_time);
half_ts = ceil(numFrames*training_set_portion);

switch training_set_creation_method
    case 'odd'
    odd_vector = 1:length(ca_time);
    training_ts = rem(odd_vector, 2) == 1;
    case 'first_portion'
    first_portion = 1:ceil(length(ca_time*training_set_portion)/2);
    training_ts(first_portion) = 1;
    case 'random'
    training_ts(1:half_ts) = 1;
    training_ts = logical(training_ts(randperm(numFrames)));
    otherwise
    error("Please define a valid method to create a training set: 'odd', 'first_portion', or 'random'");
end

training_ts = logical(training_ts);


end

