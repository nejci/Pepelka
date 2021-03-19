function [alpha,beta]=bceMstep(alpha,phi,gama,X,Q,lap)
%
% M-step of BCE
%
%   k = number of ensemble clusters
%   N = number of base clusterings
%   M = number of data points
%
% Input:
%   phi:    k*N*M, variational parameters for discrete distributions
%   gama:   k*M, variational parameters for Dirichlet distributions
%   alpha:  k*1, model parameters for the Dirichlet distribution
%   X:      M*N, base clustering results
%   Q:      cell with N elements, each is the number of clusters in base
%           clustering results
%   lap:    laplacian smoothing parameter
%
% Output:
%   alpha:  k*1, model paramter for Dirichlet distribution
%   beta:   cell with N elements for N base clusterings, each element is a
%           k*q matrix, i.e., k parameters for a q-dimensional discrete distribution   
% -------------------------------------------------

[k,N,M]=size(phi);

%-------update beta----------
for ind=1:N
    beta{ind}=zeros(k,Q(ind));
end 


for ind=1:N
  for q=1:Q(ind)
    temp=zeros(k,N);
    for s=1:M
        x=X(s,:);
        filter=(ones(k,1)*(x==q));
        temp=temp+phi(:,:,s).*filter;
    end
    beta{ind}(:,q)=temp(:,ind); 
  end
end

% smoothing
for ind=1:N
    beta{ind}=beta{ind}+lap;
    beta{ind}=beta{ind}./(sum(beta{ind},2)*(ones(1,Q(ind)))) ;
end

% -------update alpha-----------
alpha_t=alpha;
epsilon=0.001;
time=500;

t=0;
e=100;
psiGama=psi(gama);
psiSumGama=psi(sum(gama,1));
while e>epsilon&&t<time
    g=sum((psiGama-ones(k,1)*psiSumGama),2)+M*(psi(sum(alpha_t))-psi(alpha_t));
    h=-M*psi(1,alpha_t);
    z=M*psi(1,sum(alpha_t));
    c=sum(g./h)/(1/z+sum(1./h));
    delta=(g-c)./h;

    % line search
    eta=1;
    alpha_tt=alpha_t-delta;
    while (length(find(alpha_tt<=0))>0)
        eta=eta/2;
        alpha_tt=alpha_t-eta*delta;
    end
    e=sum(abs(alpha_tt-alpha_t))/sum(alpha_t);
    
    alpha_t=alpha_tt;

    t=t+1;
end
alpha=alpha_t;