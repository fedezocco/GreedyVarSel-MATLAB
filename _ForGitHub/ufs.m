%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

% This is an implementation of the algorithm proposed in:
% Whitley, D.C., Ford, M.G. and Livingstone, D.J., 2000. Unsupervised forward selection: a 
% method for eliminating redundant variables. Journal of chemical information and computer 
% sciences, 40(5), pp.1160-1168.

function selected = ufs(X, k)

[~,p] = size(X);
selected = [];

if k < 2
    error('k must be grater or equal to 2');
end

if k > p
    k = p;
end

% Selection of first 2 variables
X = zscore(X);
sqcorr = (X'*X).^2;
min_sqcorr = min(min(sqcorr));
[rmin, cmin] = ind2sub(size(sqcorr), find(sqcorr==min_sqcorr));
rmin = rmin(1);
cmin = cmin(1);
if rmin == cmin
    error('All correlations equal to one');
end
selected = [rmin cmin];

% Computation of an orthonormal basis c1 and c2 (i.e. Step 5 of the paper)
unselected = setdiff(1:p, selected); % Keeps track of the actual indices
c1 = X(:,rmin);
xb = X(:,cmin);
c2 = c1 - c1'*xb*c1; 
c2 = c2/vecnorm(c2);
c = [c1 c2];

% Selection of further variables (i.e. from Step 6)
for i = 3:k
    X_unselected = X(:, unselected);
    v_weights = zeros(1, size(X_unselected,2));
    M_projections = zeros(size(X_unselected,1), size(X_unselected,2));
    for q = 1:size(X_unselected,2)
        projection = zeros(size(c,1),1); 
        for j = 1:size(c,2)
          projection_j = c(:,j)'*X_unselected(:,q)*c(:,j);
          projection = projection + projection_j;   
        end
        v_weights(1,q) = vecnorm(projection);
        M_projections(:,q) = projection;        
    end
    
    [~, minidx] = min(v_weights);
    selected = [selected  unselected(minidx)]; % Step 7
    Xk = X(:, unselected(minidx));
    ck =  Xk - M_projections(:,minidx); % Step 8
    ck = ck/vecnorm(ck);
    c = [c ck]; 
    unselected = setdiff(1:p, selected);
end

