%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

k80fosmod = find(mean_ExpVar_fosmod >= 80, 1);
k90fosmod = find(mean_ExpVar_fosmod >= 90, 1);
k99fosmod = find(mean_ExpVar_fosmod >= 99, 1);
k95fosmod = find(mean_ExpVar_fosmod >= 95, 1);

k80fsca = find(mean_ExpVar_fsca >= 80, 1);
k90fsca = find(mean_ExpVar_fsca >= 90, 1);
k99fsca = find(mean_ExpVar_fsca >= 99, 1);
k95fsca = find(mean_ExpVar_fsca >= 95, 1);

k80LazyFsca = find(mean_ExpVar_LazyFsca >= 80, 1);
k90LazyFsca = find(mean_ExpVar_LazyFsca >= 90, 1);
k99LazyFsca = find(mean_ExpVar_LazyFsca >= 99, 1);
k95LazyFsca = find(mean_ExpVar_LazyFsca >= 95, 1);

k80fscafsfp = find(mean_ExpVar_fscafsfp >= 80, 1);
k90fscafsfp = find(mean_ExpVar_fscafsfp >= 90, 1);
k99fscafsfp = find(mean_ExpVar_fscafsfp >= 99, 1);
k95fscafsfp = find(mean_ExpVar_fscafsfp >= 95, 1);

k80GPR = find(mean_ExpVar_GPR >=80, 1);
k90GPR = find(mean_ExpVar_GPR >= 90, 1);
k99GPR = find(mean_ExpVar_GPR >= 99, 1);
k95GPR = find(mean_ExpVar_GPR >= 95, 1);

K80OPFS = find(mean_ExpVar_OPFS >= 80, 1);
K90OPFS = find(mean_ExpVar_OPFS >= 90, 1);
K99OPFS = find(mean_ExpVar_OPFS >= 99, 1);
K95OPFS = find(mean_ExpVar_OPFS >= 95, 1);

K80ufs = find(mean_ExpVar_ufs >= 80, 1);
K90ufs = find(mean_ExpVar_ufs >= 90, 1);
K99ufs = find(mean_ExpVar_ufs >= 99, 1);
K95ufs = find(mean_ExpVar_ufs >= 95, 1);

K80 = [k80fsca k80LazyFsca k80fscafsfp K80OPFS k80GPR k80fosmod K80ufs];
K90 = [k90fsca k90LazyFsca k90fscafsfp K90OPFS k90GPR k90fosmod K90ufs];
K99 = [k99fsca k99LazyFsca k99fscafsfp K99OPFS k99GPR k99fosmod K99ufs];
K95 = [k95fsca k95LazyFsca k95fscafsfp K95OPFS k95GPR k95fosmod K95ufs];
AUC = [sum(mean_ExpVar_fsca) sum(mean_ExpVar_LazyFsca) sum(mean_ExpVar_fscafsfp) ... 
    sum(mean_ExpVar_OPFS) sum(mean_ExpVar_GPR) sum(mean_ExpVar_fosmod) ...
    sum(mean_ExpVar_ufs)]/(100*length(mean_ExpVar_ufs));
