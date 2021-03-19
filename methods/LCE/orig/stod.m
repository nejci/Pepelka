function d = stod(S)
%==========================================================================
% FUNCTION: d = stod(S)
% DESCRIPTION: This function converts similarity values to distance values
%              and change matrix's format from square to vector (input
%              format for linkage function)
%
% INPUTS:   S = N-by-N similarity matrix
%
% OUTPUT:   d = a distance vector
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

s = [];
for a = 1:length(S)-1 %change matrix's format to be input of linkage fn
    s = [s S(a,[a+1:end])];
end
d = 1 - s; %compute distance (d = 1-sim)