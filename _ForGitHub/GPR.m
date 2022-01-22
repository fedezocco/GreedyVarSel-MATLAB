%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

% This is an implementation of the algorithm proposed in: 
% Krause, A., Singh, A. and Guestrin, C., 2008. Near-optimal sensor placements in 
% Gaussian processes: Theory, efficient algorithms and empirical studies. Journal of 
% Machine Learning Research, 9(2).

function selected = GPR(X, Nc, start)
%   X: the data matrix
%   Nc: the number of components to be selected
%   start: optional - the initial set of features
[n, p] = size(X);
Xnormalized = X - mean(X);
N_origninal = Nc;

if nargin < 3
    start = [];
end

if length(start) >= Nc
    selected = start;
else
    if Nc == p
        Nc = Nc - 1;
    end
    S = start;
    Sbar = setdiff(1:p, start);
    Nc = Nc - length(start);
    for i = 1:Nc
        ratios = zeros(length(Sbar), 1);
        ind = 1;
        for j = Sbar
            Sbarnew = setdiff(Sbar, j); % Removes j from the set
            xbar = Xnormalized(:,j);
            sigmax = std(xbar);
            sigmaS = Xnormalized(:,S)'*Xnormalized(:,S)/n;
            sigmaSbar = Xnormalized(:, Sbarnew)'*Xnormalized(:,Sbarnew)/n;
            num = sigmax*sigmax - (Xnormalized(:, S)'*xbar/n)'*inv(sigmaS)*(Xnormalized(:, S)'*xbar/n);
            den = sigmax*sigmax - (Xnormalized(:, Sbarnew)'*xbar/n)'*inv(sigmaSbar)*(Xnormalized(:, Sbarnew)'*xbar/n);
            ratio = num/den;
            ratios(ind) = ratio;
            ind = ind + 1;
        end
        [v, maxIdx] = max(ratios);
        S = [S Sbar(maxIdx)];
        Sbar = setdiff(Sbar, Sbar(maxIdx));
    end
    if N_origninal == p
        S = [S Sbar(1)];
    end
    selected = S;
end
end

