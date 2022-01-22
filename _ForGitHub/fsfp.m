%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

% This is an implementation of the algorithm proposed in:
% Zocco, F. and McLoone, S., 2017. Mean squared error vs. frame potential 
% for unsupervised variable selection. In Intelligent Computing, Networked Control, and Their 
% Engineering Applications (pp. 353-362). Springer, Singapore.

function selected = fsfp(X,Nc,start)

[n, p] = size(X);

if nargin < 3
    start = [];
end

if length(start) >= Nc
    selected = start;
else
    S = start;
    Sbar = setdiff(1:p, start);
    Nc = Nc - length(start);
    for i = 1:Nc
        fp = zeros(length(Sbar),1);
        ind = 1;
        for j = Sbar
            Stemp = [S j];
            fp(ind) = FP(X) - FP(X(:,Stemp));
            ind = ind+1;
        end
        [val, maxInd] = max(fp);
        S = [S Sbar(maxInd)];
        Sbar = setdiff(Sbar, Sbar(maxInd));
    end
end
selected = S;
end

