function [active_agreement_vector] = bayesian_assess_model(binarized_data, estimated_data)
%BAYESIAN_ASSESS_MODEL Summary of this function goes here
%   Detailed explanation goes here

numCells = size(binarized_data,2);

for step_i = 1:size(binarized_data,1)
   if sum(isnan(estimated_data(step_i,:))) == 0
   active_agreement_vector(step_i) = sum(binarized_data(step_i,:) == 1 & estimated_data(step_i,:) == 1)./sum(binarized_data(step_i,:) == 1);
   else
   active_agreement_vector(step_i) = nan;
   end
end

end

