function [filt_decoded_probabilities] = bayesian_temporal_filter(decoded_probabilities, ca_time, window_size)
%BAYESIAN_TEMPORAL_FILTER Summary of this function goes here
%   Detailed explanation goes here

Fs = 1/mode(diff(ca_time));
window_size = round(window_size*Fs);
half_win = round(window_size/2);

%% Pad matrix with zeros
zero_pad = zeros(size(decoded_probabilities,1),window_size);
padded_decoded_probabilities = [zero_pad decoded_probabilities zero_pad];

for step_i = 1:size(decoded_probabilities,2)
current_window = padded_decoded_probabilities(:,step_i+half_win: step_i-1+window_size+half_win);
filt_decoded_probabilities(:,step_i) = expm1(sum(log1p(current_window),2,'omitnan'));
filt_decoded_probabilities(:,step_i) = filt_decoded_probabilities(:,step_i)./sum(filt_decoded_probabilities(:,step_i));
end

end

