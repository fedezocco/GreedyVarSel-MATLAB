%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

% This is an implementation of the algorithm proposed in: 
% Wei, H.L. and Billings, S.A., 2006. Feature subset selection and ranking for data 
% dimensionality reduction. IEEE transactions on pattern analysis and machine intelligence, 
% 29(1), pp.162-166.


function [v_SelectedVariables] = fosmod(X_data, DesiredNumOfVariables) 

%% ALGORITHM:
if max(abs(mean(X_data))) > 10^-6 
   [X_data, ~] = pca_normalise(X_data); 
end

%% INITIALISATION:
[~,TotalNumOfVariables] = size(X_data);
v_weights = [];
X_data_unselected = X_data;
v_SelectedVariables = [];
v_trueIndex = 1:TotalNumOfVariables;

%%% First step (eq. (4), (5) and (6)):
C_ini = corr(X_data,X_data); % eq. (4)
meanCj = (sum(C_ini.^2,1))./TotalNumOfVariables; % eq. (5)
[~,SelectedVar_Index] = max(meanCj); % eq. (6)
v_SelectedVariables(1,1) = SelectedVar_Index;
X_data_unselected(:,SelectedVar_Index) = [];
Q_tilde = X_data(:,SelectedVar_Index);
v_trueIndex = setdiff(v_trueIndex,v_trueIndex(:,SelectedVar_Index));

%%% from the 2nd step:
if (DesiredNumOfVariables > 1)    
    
    for i = 2:DesiredNumOfVariables

        for j = 1:size(X_data_unselected,2) 
            alpha_j = X_data_unselected(:,j);
            
            for k = 1:size(Q_tilde,2) % eq. (7)
                w_j = ((alpha_j'*Q_tilde(:,k))*Q_tilde(:,k))/(Q_tilde(:,k)'*Q_tilde(:,k));
                v_wj(:,k) = w_j;
            end
            w = sum(v_wj,2);
            q_j = alpha_j - w; 
            v_Cj = corr(X_data,q_j);  % eq. (8)
            meanC_j = (sum(v_Cj.^2,1))/TotalNumOfVariables; % eq. (9) 
            v_weights(1,j) = meanC_j;
            w_j = [];
            q_j = [];
        end

        [~, SelectedVar_FalseIndex] = max(v_weights); % eq. (10) 
        v_SelectedVariables(1,i) = v_trueIndex(:,SelectedVar_FalseIndex);

        
        % Next step preparation:  
        v_weights = [];
        v_trueIndex = setdiff(v_trueIndex, v_trueIndex(:,SelectedVar_FalseIndex));
        
        for k = 1:size(Q_tilde,2) % eq. (7) - update of the matrix Q_tilde 
                w_m = ((X_data_unselected(:,SelectedVar_FalseIndex)'*Q_tilde(:,k))*Q_tilde(:,k))/(Q_tilde(:,k)'*Q_tilde(:,k));
                v_wm(:,k) = w_m;
        end
        
        w_m = sum(v_wm,2);
        q_m = (X_data_unselected(:,SelectedVar_FalseIndex)) - w_m;
        Q_tilde = [Q_tilde q_m];
        X_data_unselected(:,SelectedVar_FalseIndex) = [];
    end
end

