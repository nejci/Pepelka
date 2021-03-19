function dvi = voi(L1,L2)
%VOI Variation of information.
% DVI = VOI(L1,L2) returns the variation of information shared by two
% N-by-1 integer arrays of classification data, L1 and L2.
%
% DVI is 0 for identical clusterings and 1 for completely different
% clusterings.
%
% Copyright (2009) Sandia Corporation. Under the terms of Contract 
% DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
% certain rights in this software.

% entropies for individual classification arrays
h1 = infoentropy(L1);
h2 = infoentropy(L2);

% mutual information shared by the classification arrays
i12 = mutualinfo(L1,L2);

% variation of inofrmation
dvi = h1 + h2 - 2*i12;

function H = infoentropy(L)
%INFOENTROPY Entropy of information.
% H = INFOENTROPY(L) returns the entropy of information for an N-by-1
% integer array of classification data, L. 
%
% Copyright (2009) Sandia Corporation. Under the terms of Contract 
% DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
% certain rights in this software.

N = length(L);
k = unique(L);
H = 0;
  
% loop over the unique elements of L
for i = 1:size(k,1)
    % the probability of a given classification index occurring in L
    pk = sum(k(i) == L)/N;
    H = H - pk*log(pk);
end

function I = mutualinfo(L1,L2)
%MUTUALINFO Mutual information.
% I = MUTUALINFO(L1,L2) returns the mutual information shared by two N-by-1
% integer arrays of classification data, L1 and L2. 
%
% Copyright (2009) Sandia Corporation. Under the terms of Contract 
% DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
% certain rights in this software.

N = length(L1);
k1 = unique(L1);
k2 = unique(L2);
I = 0;

% loop over the unique elements of L1
for i = 1:size(k1,1)
    % loop over the unique elements of L2
    for j = 1:size(k2,1)
        % the mutual probability of two classification indices occurring in
        % L1 and L2
        pij = sum((L1 == k1(i)).*(L2 == k2(j)))/N;
        % the probability of a given classification index occurring in L1
        pi = sum(L1 == k1(i))/N;
        % the probability of a given classification index occurring in L2
        pj = sum(L2 == k2(j))/N;
        if (pij > 0)
            I = I + pij*log(pij/(pi*pj));
        end
    end
end