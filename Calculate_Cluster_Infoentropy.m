function H = calcualate_infoentropy(L)
%INFOENTROPY Entropy of information.
% H = INFOENTROPY(L) returns the entropy of information for an N-by-1
% integer array of classification data, L. 
%
% Copyright (2009) Sandia Corporation. Under the terms of Contract 
% DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
% certain rights in this software.

% Demo
% L=[2 2 1 3 4 4 4 4]
% L=L';

N = length(L);
k = unique(L);
H = 0;
  
% loop over the unique elements of L
for i = 1:size(k,1)
    % the probability of a given classification index occurring in L
    pk = sum(k(i) == L)/N;
    H = H - pk*log(pk);
end