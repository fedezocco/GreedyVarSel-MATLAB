%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

function fp = FP(X)
%Compute the frame potential of the matrix X 
% X: input matrix

[n, p] = size(X);

% Normalize the matrix to have unit norm columns
for j = 1:p
    X(:, j) = X(:, j)./norm(X(:, j));
end

corr = X'*X; % Correlation of the columns of the matrix
fp = sum(sum(corr.^2));
end

