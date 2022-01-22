%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

for j = 1:n_vars
    
    disp(['Number of variables ' num2str(j) '/' num2str(size(NumOfSelectedVariables_consideredCases,2))]);
    %% Selection step:
    % (1) FSCA algorithm:
    tic;
    [~,~,~,list_SelectedVariables_fsca] = fsca(X_data,j);
    compTime_fsca(j) = toc;
    
    % (2) Lazy-greedy implementation of FSCA:
    tic;
    [~,~,~,list_SelectedVariables_LazyFsca] = lazy_fsca(X_data,j); 
    compTime_LazyFsca(j) = toc;
    
    % (3) FSFP-FSCA algorithm:
    tic;
    [list_SelectedVariables_fscafsfp] = fsca_fsfp(X_data,j);
    compTime_fscafsfp(j) = toc;
    
    % (4) OPFS algorithm:
    tic;
    [~,~,~,list_SelectedVariables_OPFS] = OPFS(X_data,j);
    compTime_OPFS(j) = toc;
    
    % (5) GPR (Krause-based) algorithm:
    tic;
    [list_SelectedVariables_GPR] = GPR(X_data,j);
    compTime_GPR(j) = toc;
    
    % (6) FOS-MOD algorithm:
    tic;
    [list_SelectedVariables_fosmod] = fosmod(X_data,j);
    compTime_fosmod(j) = toc;
    
    % (7) UFS algorithm:
    tic;
    if j > 1
        [list_SelectedVariables_ufs] = ufs(X_data,j);
    else
        [list_SelectedVariables_ufs] = ufs(X_data,n_vars);
    end
    compTime_ufs(j) = toc;

    
    %% Selection performance:
    
    % "mean_ExpVar_pca" can be evaluated with pca_svd.m. Currently, variable
    % selection algorithms only are compared.  

    % (1) FSCA:
    X_fsca = X_data(:, list_SelectedVariables_fsca(1:j));
    X_unsel = X_data; X_unsel(:, list_SelectedVariables_fsca(1:j)) = [];
    mean_ExpVar_fsca(j) = mean_ExpVar_fsca(j) + PercExplainedVariance(X_data,X_fsca);
    mean_MI_fsca(j) = mean_MI_fsca(j) + mutual_info(X_fsca, X_unsel);
    mean_FP_fsca(j) = mean_FP_fsca(j) + FP(X_fsca);
    % (2) Lazy FSCA:
    X_LazyFsca = X_data(:, list_SelectedVariables_LazyFsca(1:j));
    X_unsel = X_data; X_unsel(:, list_SelectedVariables_LazyFsca(1:j)) = [];
    mean_ExpVar_LazyFsca(j) = mean_ExpVar_LazyFsca(j) + PercExplainedVariance(X_data,X_LazyFsca);
    mean_MI_LazyFsca(j) = mean_MI_LazyFsca(j) + mutual_info(X_LazyFsca, X_unsel);
    mean_FP_LazyFsca(j) = mean_FP_LazyFsca(j) + FP(X_LazyFsca);
    % (3) FSFP-FSCA:
    X_fscafsfp = X_data(:, list_SelectedVariables_fscafsfp(1:j));
    X_unsel = X_data; X_unsel(:, list_SelectedVariables_fscafsfp(1:j)) = [];
    mean_ExpVar_fscafsfp(j) = mean_ExpVar_fscafsfp(j) + PercExplainedVariance(X_data,X_fscafsfp);
    mean_MI_fscafsfp(j) = mean_MI_fscafsfp(j) + mutual_info(X_fscafsfp, X_unsel);
    mean_FP_fscafsfp(j) = mean_FP_fscafsfp(j) + FP(X_fscafsfp);
    % (4) OPFS:
    X_OPFS = X_data(:, list_SelectedVariables_OPFS(1:j));
    X_unsel = X_data; X_unsel(:, list_SelectedVariables_OPFS(1:j)) = [];
    mean_ExpVar_OPFS(j) = mean_ExpVar_OPFS(j) + PercExplainedVariance(X_data,X_OPFS);
    mean_FP_OPFS(j) = mean_FP_OPFS(j) + FP(X_OPFS);
    mean_MI_OPFS(j) = mean_MI_OPFS(j) + mutual_info(X_OPFS, X_unsel);
    % (5) GPR (i.e. ITFS):
    X_GPR = X_data(:, list_SelectedVariables_GPR(1:j));
    X_unsel = X_data; X_unsel(:, list_SelectedVariables_GPR(1:j)) = [];
    mean_ExpVar_GPR(j) = mean_ExpVar_GPR(j) + PercExplainedVariance(X_data,X_GPR);
    mean_FP_GPR(j) = mean_FP_GPR(j) + FP(X_GPR);
    mean_MI_GPR(j) = mean_MI_GPR(j) + mutual_info(X_GPR, X_unsel);
    % (6) FOS-MOD:
    X_fosmod = X_data(:, list_SelectedVariables_fosmod(1:j));
    X_unsel = X_data; X_unsel(:, list_SelectedVariables_fosmod(1:j)) = [];
    mean_ExpVar_fosmod(j) = mean_ExpVar_fosmod(j) + PercExplainedVariance(X_data,X_fosmod);
    mean_FP_fosmod(j) = mean_FP_fosmod(j) + FP(X_fosmod);
    mean_MI_fosmod(j) = mean_MI_fosmod(j) + mutual_info(X_fosmod, X_unsel);
    % (7) UFS:
    X_ufs = X_data(:, list_SelectedVariables_ufs(1:j));
    X_unsel = X_data; X_unsel(:, list_SelectedVariables_ufs(1:j)) = [];
    mean_ExpVar_ufs(j) = mean_ExpVar_ufs(j) + PercExplainedVariance(X_data,X_ufs);
    mean_FP_ufs(j) = mean_FP_ufs(j) + FP(X_ufs);
    mean_MI_ufs(j) = mean_MI_ufs(j) + mutual_info(X_ufs, X_unsel);
end

