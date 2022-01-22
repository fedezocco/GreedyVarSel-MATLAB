%% (c) Sean McLoone, s.mcloone@ieee.org, Nov 2017

% This is an implementation of the algorithm proposed in: 
% Cui, Y. and Dy, J.G., 2008. Orthogonal principal feature selection.

function [S, M, VarEx, compId]=OPFS_SVD(X, Nc, ~);

% [M, F, VarEx, CompId]=OPFS_SVD(X, Nc);
% X = n x m matrix (measurements x variables)
% NB: X is assumed to be normalised -- i.e. mean of each column is zero.
% Nc = number of components to select


% Algorithm requires X to have zero mean columns
mX=mean(X);
if max(abs(mX)) > 10^-6
    %if nargin <3 fprintf('\nWarning: Data not zero mean ... detrending\n'); end
    X=X-ones(size(X,1),1)*mX;
end

L=size(X,2);

Y=X;
VT=var(Y(:));
compId=[];
VarEx=[];
YhatP=0;
S=[];
M=[];

if nargin <2
    Nc=1;
end

PVE=zeros(Nc,1);

for j=1:Nc
    [p1,t1,ic]=pca_svd(Y,1);
    
    EFS=zeros(L,1);
    for i=1:L;
        x=Y(:,i);
        EFS(i)=((x'*t1)^2)/((x'*x)+eps);
    end
    [v Id]=max(EFS);
    x=Y(:,Id);
 
    th=pinv(x)*Y;
    Yhat= x*th;
    Y=Y-Yhat;
    
    S=[S x];
    M=[M  th'];
    YhatP=YhatP+Yhat; 
    compId=[compId Id];
    VarEx=[VarEx var(YhatP(:))/VT*100];

end


