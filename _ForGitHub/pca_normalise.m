%% (c) Sean McLoone, March 2006

function [d, m, s]= pca_normalise(x,flag)

% [d, m, s]= pca_normalise(x,flag)
%
% Removes the mean from each of the columns in x (and scales each column by
% the standard deviation if flag >0)
% Returns the normalised data in d, the means in m and the standard
% deviations in s.
%  x is N data values by M variables
%  m is a 1 x M vector
%  s is a 1 x M vector

if nargin <2
    flag=0;
end

m=mean(x);
s=std(x)+2*eps; %to avoid problems with s=0.

    
d=x-ones(size(x,1),1)*m;

if flag >0
    d=d./(ones(size(x,1),1)*s);
end
