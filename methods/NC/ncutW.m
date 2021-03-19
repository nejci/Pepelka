function [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(W,nbcluster,options)
% [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(W,nbcluster);
% 
% Calls ncut to compute NcutEigenvectors and NcutEigenvalues of W with nbcluster clusters
% Then calls discretisation to discretize the NcutEigenvectors into NcutDiscrete
% Timothee Cour, Stella Yu, Jianbo Shi, 2004

% Added by Nejc:
% options.offset = 5e-1; %offset in the diagonal of W
% options.verbose = 0; %0 for verbose off mode, 1,2,3 for verbose on modes
% options.maxiterations = 100; %max number of iterations in eigensolver
% options.eigsErrorTolerance = 1e-6; %error tolerance in eigensolver
% options.valeurMin=1e-6; %truncates any values in W less than valeurMin

% compute continuous Ncut eigenvectors
[NcutEigenvectors,NcutEigenvalues] = ncut(W,nbcluster,options);

% compute discretize Ncut vectors
[NcutDiscrete,NcutEigenvectors] = discretisation(NcutEigenvectors);


NcutDiscrete = full(NcutDiscrete);