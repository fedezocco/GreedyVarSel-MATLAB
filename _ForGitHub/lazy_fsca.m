%% (c) Seán McLoone, s.mcloone@ieee.org, updated December 2021 (V3.0)

function [SelComp, CompWeights, VarEx, SelVars]=lazy_fsca(X, Nc,~);

% [S, M, VarEx, CompId]=lazy_fsca(X, Nc);
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
% of the previously selected components has been removed, as approximated
% by the lazy greedy search algorithm.
%
% Note: This is equivalent to selecting successive components as those
% which in combination with the previously selected components explain the
% most variance across all the data. (i.e the implementation of FSCA given in FSCV)
% 
% References: 
% [1] F. Zocco, M. Maggipinto, G. A. Susto, S. McLoone; 2022; "Lazy FSCA for unsupervised variable 
% selection"; https://arxiv.org/abs/2103.02687.
% [2] L. Puggini, S. McLoone, Forward Selection Component Analysis: Algorithms
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
    if nargin<3 fprintf('\nWarning: Data not zero mean ... detrending\n');  end
    X=X-ones(size(X,1),1)*mX;
end

if nargin <3 
    fprintf('\nWarning: Data not zero mean ... detrending\n'); 
    mX=mean(X);
    X=X-ones(size(X,1),1)*mX;
end
    


Nvars=size(X,2);

Y=X;
TR=trace(Y'*Y);

%First component
 for i=1:Nvars;
        x=Y(:,i);
        r=Y'*x;
        EFS(i)=r'*r/(x'*x);  %Rayleigh quotient for x_i  
 end
 VEX=100*(EFS/TR);  %Convert to variance explained.
 
 [dGs Is]=sort(VEX,'descend');  %order from max to min 
 SelVars=Is(1);  %first entry is the optimum one
 
 %deflate matrix
 x=Y(:,SelVars);
 th=pinv(x)*Y;
 YhatP= x*th;
 Y=Y-YhatP; 

 %initalise storage variables 
 SelComp=x; 
 CompWeights=th';
 VarEx=dGs(1);  % equivalent to VarEx=var(YhatP(:))/VT*100;

for j=2:Nc
    top=j; IsB=Is(top); % used to track how far through the list we are
    dGsB=0; loop=1; dGsL=VarEx(1);
    while loop==1
        x=Y(:,Is(top));  
        r=Y'*x;
        dGs(top)=(100/TR)*(r'*r)/(x'*x);  %increment in variance contribution with inclusion of this variable

        if dGs(top) > dGsB   %remember best exact result and its location
            dGsB=dGs(top);
            IsB=Is(top);
        end
       
        if dGs(top) < dGsL   %remember worst exact result and its location
            dGsL=dGs(top);
        end
         
        %If the best exact increment found to date is not greater than the upper bound on the
        %next variable contribution move to the next variable and evaluate
        %it exactly
        
        if top == Nvars    %reached the end of the line -- best index in the already checked entries.
            loop=0; 
        else
            if dGsB>dGs(top+1)
                loop=0; %found the best value - stop searching
            else
                top=top+1; % check the next variable in the list
            end
        end
    end
    
    %Deflate matrix for selected x
    x=Y(:,IsB);
    th=pinv(x)*Y;
    Yhat= x*th;
    Y=Y-Yhat;
   
    
    %store results
    SelComp=[SelComp x];
    CompWeights=[CompWeights  th'];
    SelVars=[SelVars IsB];
    VarEx=[VarEx dGsB+VarEx(end)];
    
   %resort the top of the list
   if top < Nvars
      while dGsL < dGs(top+1)  %first find range that needs to be resorted
         top=top+1;
         if top == Nvars
             break;
         end
     end
   end
   [dGs(1:top) Inew]=sort(dGs(1:top),'descend');    
   Is(1:top)=Is(Inew);
   
end


end


