function [filt_decoded_probabilities] = bayesian_temporal_filter2D(decoded_probabilities, ca_time, window_size)
%BAYESIAN_TEMPORAL_FILTER Summary of this function goes here
%   Detailed explanation goes here

Fs = 1/mode(diff(ca_time));
window_size = round(window_size*Fs);
half_win = round(window_size/2);

%% Pad matrix with zeros
zero_pad = zeros(size(decoded_probabilities,1),size(decoded_probabilities,2),window_size);
padded_decoded_probabilities = cat(3,zero_pad,decoded_probabilities);
padded_decoded_probabilities = cat(3,padded_decoded_probabilities,zero_pad);

for step_i = 1:size(decoded_probabilities,3)
current_window = padded_decoded_probabilities(:,:,step_i+half_win: step_i-1+window_size+half_win);
temp_filt_decoded_probabilities = expm1(sum(log1p(current_window),3,'omitnan'));
filt_decoded_probabilities(:,:,step_i) = temp_filt_decoded_probabilities./sum(temp_filt_decoded_probabilities, 'omitnan');

end

end

