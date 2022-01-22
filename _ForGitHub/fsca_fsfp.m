%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

% In this algorithm, the first variable is selected by FSCA (see
% "fsca.m"), while the remaining variables are selected by FSFP (see "fsfp.m").

function [v_SelectedVariables] = fsca_fsfp(X_data, DesiredNumOfVariables)

[~,~,~,v_SelectedVariables_fsca] = fsca(X_data,1);

if (DesiredNumOfVariables > 1)
    v_SelectedVariables = fsfp(X_data, DesiredNumOfVariables, v_SelectedVariables_fsca); 
else
    v_SelectedVariables = v_SelectedVariables_fsca;
end
 
