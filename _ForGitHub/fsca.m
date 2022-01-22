%% (c) Seán McLoone, s.mcloone@ieee.org, updated December 2021 (V3.0)

function [SelComp, CompWeights, VarEx, SelVars]=fsca(X, Nc, ~);

% [S, M, VarEx, CompId]=FSCA(X, Nc);
% X = m x v matrix (measurements x variables)
% NB: X is assumed to be normalised -- i.e. mean of each column is zero.
% Nc = number FSC components to compute
%
% X ~= S*M'; where S are FS components and M is the linear weightings
% VarEx gives the accumulative variance explained and CompId lists the
% coordinates of the corresponding variables.
%
% Computes FSC decompostion by selecting successive components as those
% which explain the most variance across all the data after the contributions
% of the previously selected components has been removed.
%
% Note: This is equivalent to selecting successive components as those
% which in combination with the previously selected components explain the
% most variance across all the data. (i.e the implementation of FSCA given in FSCV)
% 
% Reference: 
% L. Puggini, S. McLoone, Forward Selection Component Analysis: Algorithms
% and Applications, IEEE Transactions on Pattern Analysis and Machine Intelligence, 
% Vol. 39(12), pp. 2395-2408, December 2017, DOI: 10.1109/TPAMI.2017.2648792.



if nargin <2
    Nc=1;
end

if Nc > min(size(X))
    Nc = min(size(X));
end

% Algorithm requires X to have zero mean columns
mX=mean(X);
if max(abs(mX)) > 10^-6
    if nargin <3 fprintf('\nWarning: Data not zero mean ... detrending\n'); end
    X=X-ones(size(X,1),1)*mX;
end

Nvars=size(X,2);

Y=X;
TR=trace(Y'*Y); %Used to compute variance explained

%initalise storage variables
SelVars=[];
VarEx=[];
YhatP=0;
SelComp=[];
CompWeights=[];
VEX=0;

EFS=zeros(Nvars,1); %intialise storage for Rayleigh quotient values.
AvailVar=1:1:Nvars;  % Variables to select from

for j=1:Nc 
    for i=1:1:length(AvailVar)
        x=Y(:,AvailVar(i));
        r=Y'*x;
        EFS(i)=r'*r/(x'*x);  %Rayleigh quotient for x_i
    end
    [v Id]=max(EFS); %select variable with max RQ.

    VEX=VEX+100*v/TR; % accumulated variance explained 
                      % 100*v/TR = variance expalined by selected component
                      
    %Deflate matrix
    x=Y(:,AvailVar(Id));
    th=pinv(x)*Y;
    Yhat= x*th;
    Y=Y-Yhat;
    
    %store results and prepare for next iteration
    SelComp=[SelComp x];
    CompWeights=[CompWeights  th'];
    VarEx=[VarEx VEX];
    SelVars=[SelVars AvailVar(Id)];
    AvailVar(Id)=[];
    EFS(Id)=[]; 
    
end

end


