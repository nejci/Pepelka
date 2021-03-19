function [phiAll,gamaAll,resultAlpha,resultBeta]=learnBCE(X,oldAlpha,oldBeta,lap,Q)
%
% BCE learning
% 
%   k = number of ensemble clusters
%   N = number of base clusterings
%   M = number of data points
%
% Input:
%   X:          M*N, base clustering results
%   oldAlpha:   k*1, model parameter for Dirichlet distribution
%   oldBeta:    cell with N elements for N base clusterings, each element is a
%               k*q matrix, i.e., k parameters for a q-dimensional discrete distribution
%   lap:        smoothing parameter
%   Q:          cell with N elements, each is the number of clusters in base
%               clustering results         
% Output:
%   phiAll:     k*N*M, variational parameters for discrete distributions
%   gamaAll:    k*M, variational parameters for Dirichlet distributions
%--------------------------------------------------------------------

[M,N] = size(X);
k=length(oldAlpha);

% initial value and variables for iteration
alpha_t=oldAlpha;
beta_t=oldBeta;
epsilon=0.01;
time=500;
e=100;
t=1;

% start learning iterations
%disp(['learning BCE'])
while e>epsilon && t<time
    
    % E-step
    for s=1:M
        sample=X(s,:);
        [estimatedPhi,estimatedGama]=bceEstep(alpha_t,beta_t,sample);
        phiAll(:,:,s)=estimatedPhi;
        gamaAll(:,s)=estimatedGama;    
    end

    % M-step
    [alpha_tt,beta_tt]=bceMstep(alpha_t,phiAll,gamaAll,X,Q,lap);
    
    % error
    upvalue=0;downvalue=0;
    for index=1:length(Q)
        upvalue=upvalue+sum(sum(abs(beta_t{index}-beta_tt{index})));
        downvalue=downvalue+sum(sum(beta_t{index}));
    end
    e=upvalue/downvalue;
    %disp(['t=',int2str(t),', error=',num2str(e)]);
  
    % update
    alpha_t=alpha_tt;
    beta_t=beta_tt;
        
    t=t+1;

end

resultAlpha=alpha_t;
resultBeta=beta_t;

