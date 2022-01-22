%% Implemented by Marco Maggipinto and Federico Zocco; last update: 21/01/2022

plot(1:n_vars,mean_ExpVar_fsca(1:n_vars),'-ob','linewidth',5,'markersize',7,'MarkerFaceColor',[0 0 1]);
hold on;
plot(1:n_vars,mean_ExpVar_LazyFsca(1:n_vars),'-.ok','linewidth',5,'markersize',7,'MarkerFaceColor',[0 0 1]);
hold on;
plot(1:n_vars,mean_ExpVar_fscafsfp(1:n_vars),'-sg','linewidth',5,'markersize',7,'MarkerFaceColor',[0 1 0]);
hold on;
plot(1:n_vars,mean_ExpVar_OPFS(1:n_vars),':dc','linewidth',5,'markersize',7,'MarkerFaceColor',[0.75 0 0.75]);
hold on;
plot(1:n_vars,mean_ExpVar_GPR(1:n_vars),'-.*y','linewidth',5,'markersize',7,'MarkerFaceColor',[1 1 0]);
hold on;
plot(1:n_vars,mean_ExpVar_fosmod(1:n_vars),'--om','linewidth',5,'markersize',7);
hold on;
plot(1:n_vars,mean_ExpVar_ufs(1:n_vars),'-or','linewidth',5,'markersize',7);
hold on;
grid on
% axis([1 15 0 110]);
set(gca,'xtick',1:3:n_vars,'fontsize',30)
xlabel(['Number of variables out of ' num2str(size(X_dataStore,2))],'FontSize',30)
ylabel('Explained variance','FontSize',30)
leg = legend({'FSCA','L-FSCA','FP-FSCA','OPFS','ITFS','FOS-MOD', 'UFS'},'FontSize',30);
xlim([1 n_vars]);
