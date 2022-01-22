%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

function ExplainedVariance_value = PercExplainedVariance(X_data,X_reduced)

if max(abs(mean(X_data))) > 10^-6 
    [X_data, ~] = pca_normalise(X_data); 
end

if max(abs(mean(X_reduced))) > 10^-6 
    [X_reduced, ~] = pca_normalise(X_reduced); 
end


NumElementsOfX_data = size(X_data,1) * size(X_data,2);
Theta = pinv(X_reduced) * X_data;
Xhat = X_reduced * Theta;
MSEvalue = sum(sum((X_data - Xhat).^2))/NumElementsOfX_data;
ExplainedVariance_value = (1 - (MSEvalue / var(reshape(X_data,[],1),1))) * 100;