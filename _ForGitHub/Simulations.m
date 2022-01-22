%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

%% Source code of the following paper:
% F. Zocco, M. Maggipinto, G. A. Susto, S. McLoone; 2022; "Lazy FSCA for unsupervised variable 
% selection"; https://arxiv.org/abs/2103.02687.

%% The script performs Monte Carlo simulations on the dataset "datasetName" selected at the beginning.
% The datasets "Sim1" and "Sim2" are generated by "DataGenerator.m". "PitProps" is uploaded on this repository. 
% The datasets "Semiconductor" and "MofERdata" (called
% "PlasmaEtch" in the source paper) are not public for confidentiality
% reasons. The remaining datasets are not uploaded on this repository due to 
% memory constraints. The user can download them from the links provided in
% the source paper. 

% REFERENCES:
% [1] F. Zocco, M. Maggipinto, G. A. Susto, S. McLoone; 2022; "Lazy FSCA for unsupervised variable 
% selection"; https://arxiv.org/abs/2103.02687.    
% [2] L. Puggini and S. McLoone, "Forward selection component analysis:
% Algorithms and applications," IEEE Transactions on Pattern Analysis and
% Machine Intelligence, vol. 39, no. 12, pp. 2395-2408, 2017.
% [3] Zocco, F. and McLoone, S., 2017. Mean squared error vs. frame potential 
% for unsupervised variable selection. In Intelligent Computing, Networked Control, and Their 
% Engineering Applications (pp. 353-362). Springer, Singapore.
% [4] Cui, Y. and Dy, J.G., 2008. Orthogonal principal feature selection.
% [5] Krause, A., Singh, A. and Guestrin, C., 2008. Near-optimal sensor placements in 
% Gaussian processes: Theory, efficient algorithms and empirical studies. Journal of 
% Machine Learning Research, 9(2).
% [6] Wei, H.L. and Billings, S.A., 2006. Feature subset selection and ranking for data 
% dimensionality reduction. IEEE transactions on pattern analysis and machine intelligence, 
% 29(1), pp.162-166.
% [7] Whitley, D.C., Ford, M.G. and Livingstone, D.J., 2000. Unsupervised forward selection: a 
% method for eliminating redundant variables. Journal of chemical information and computer 
% sciences, 40(5), pp.1160-1168.
%=================================================================================================


clear all, warning off;

%% Experimental setup:
datasetNames = {'Sim1', 'Sim2', 'Medical', 'Semiconductor', 'PitProps', 'Music', 'GasBatch', 'Financial', 'MofERdata', 'wafer', 'YaleB', 'USPS'};
TestSelector = 1; % An integer > 1. See DataGenerator.m for details.
datasetName = datasetNames{TestSelector};
NumOfSimulationsPerCase = 2;  % Number of X_data to simulate. The outputs of this script are mean values obtained cosidering the results of the NumOfSimulationsPerCase X_data.
X_dataStore = DataGenerator(TestSelector,NumOfSimulationsPerCase);   
n_vars = 7;
        
  


%% =====================================================================================
NumOfSelectedVariables_consideredCases = 1:n_vars;

% Explained variance performance:
mean_ExpVar_fsca = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_ExpVar_LazyFsca = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_ExpVar_fscafsfp = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_ExpVar_OPFS = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_ExpVar_GPR = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_ExpVar_fosmod = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_ExpVar_pca = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_ExpVar_ufs = zeros(1,size(NumOfSelectedVariables_consideredCases,2));

mean_MI_fsca = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_MI_LazyFsca = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_MI_fscafsfp = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_MI_OPFS = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_MI_GPR = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_MI_fosmod = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_MI_ufs = zeros(1,size(NumOfSelectedVariables_consideredCases,2));

mean_FP_fsca = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_FP_LazyFsca = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_FP_fscafsfp = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_FP_OPFS = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_FP_GPR = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_FP_fosmod = zeros(1,size(NumOfSelectedVariables_consideredCases,2));
mean_FP_ufs = zeros(1,size(NumOfSelectedVariables_consideredCases,2));


% VARIABLE SELECTION 
if TestSelector > 2  % A real dataset is selected; NOTE: no train/test set split is done, hence train == test
    X_data = X_dataStore;
    compare;
else
    NumOfSimulations = size(X_dataStore,3);
    for i = 1:NumOfSimulations
        X_data = X_dataStore(:,:,i);
        compare;
    end
    % Variance explained:
    mean_ExpVar_fsca = mean_ExpVar_fsca/NumOfSimulations;
    mean_ExpVar_LazyFsca = mean_ExpVar_LazyFsca/NumOfSimulations;
    mean_ExpVar_fscafsfp = mean_ExpVar_fscafsfp/NumOfSimulations;
    mean_ExpVar_OPFS = mean_ExpVar_OPFS/NumOfSimulations;
    mean_ExpVar_GPR = mean_ExpVar_GPR/NumOfSimulations;
    mean_ExpVar_fosmod = mean_ExpVar_fosmod/NumOfSimulations;
    mean_ExpVar_pca = mean_ExpVar_pca/NumOfSimulations;
    mean_ExpVar_ufs = mean_ExpVar_ufs/NumOfSimulations;
    % Mutual information:
    mean_MI_fsca = mean_MI_fsca/NumOfSimulations;
    mean_MI_LazyFsca = mean_MI_LazyFsca/NumOfSimulations;
    mean_MI_fscafsfp = mean_MI_fscafsfp/NumOfSimulations;
    mean_MI_OPFS = mean_MI_OPFS/NumOfSimulations;
    mean_MI_GPR = mean_MI_GPR/NumOfSimulations;
    mean_MI_fosmod = mean_MI_fosmod/NumOfSimulations;
    mean_MI_ufs = mean_MI_ufs/NumOfSimulations;
    % Frame potential:
    mean_FP_fsca = mean_FP_fsca/NumOfSimulations;
    mean_FP_LazyFsca = mean_FP_LazyFsca/NumOfSimulations;
    mean_FP_fscafsfp = mean_FP_fscafsfp/NumOfSimulations;
    mean_FP_OPFS = mean_FP_OPFS/NumOfSimulations;
    mean_FP_GPR = mean_FP_GPR/NumOfSimulations;
    mean_FP_fosmod = mean_FP_fosmod/NumOfSimulations;
    mean_FP_ufs = mean_FP_ufs/NumOfSimulations;
end


%% Further metrics:
compute_k99AUC;


%% Plots
Plots;
