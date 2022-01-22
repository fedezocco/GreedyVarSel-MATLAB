%% (c) Sean McLoone, s.mcloone@ieee.org, Nov 2017

% REFERENCES:
% [1] L. Puggini and S. McLoone, "Forward selection component analysis:
% Algorithms and applications," IEEE Transactions on Pattern Analysis and
% Machine Intelligence, vol. 39, no. 12, pp. 2395-2408, 2017.
% [2] I. T. Jolliffe, Principal component analysis; Second edition.
% Springer, 2002.

function [P, T, VE]= pca_svd(X,Nc)

% [P, T, VE]=pca_svd(X,Nc)
%
% Principal Component Decomposition of X, i.e.  X=T*P' computed using SVD
% where X is X = m measurements x v variables
% Nc is the max number of principal components to include in the
% decomposition - if not specified defaults to r =rank(X).
% T = scores matrix  (m x min(Nc,r)
% P = loading matrix (v x min(Nc,r)
% VE = Variance explained
%
% NB: This assumes X = m x v, 
% where m= no of samples and v = no of variables


[U,S,V]=svd(X,0);
P=V(:,1:Nc);
T=U*S;
T=T(:,1:Nc);
s=diag(S);
VE=cumsum(s.*s)/sum(s.*s)*100;
VE=VE(1:Nc);


