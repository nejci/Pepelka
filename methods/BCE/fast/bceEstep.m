function [phi_t,gama_t]=bceEstep(alpha,beta,x)
% 
% E step of BCE
%
%   k = number of ensemble clusters
%   N = number of base clusterings
%   M = number of data points
%
% Input:
%   alpha:  k*1, parameter for Dirichlet distribution
%   beta:   cell with N elements for N base clusterings, each element is a
%           k*q matrix, i.e., k parameters for a q-dimensional discrete distribution
%   x:      1*N, base clustering results for one data point. 0 indicates
%           missing base clustering results
%
% output:
%   phi_t:  k*N, variational parameter for discrete distribution
%   gama_t: k*1, variational parameter for Dirichlet distribution
%-------------------------------------------------


k=length(alpha);
N=size(x,2);
mask = x~=0;
V=sum(mask); 
filter=ones(k,1)*mask; 

% initial value for variational parameters
phi_t=ones(k,N)/k.*filter;
gama_t=alpha+V/k;

% variables for iteration
epsilon=0.01;
time=500;
e=100;
t=1;


% for i=1:k
%     for n=1:N
%         if x(n)~=0
%            tempBeta(i,n)=beta{n}(i,x(n));
%         else 
%            tempBeta(i,n)=-1;
%         end
%     end
% end

tempBeta = ones(k,N)*(-1);
for i=1:k
    for n=1:N
        tempBeta(i,n)=beta{n}(i,x(n));
    end
end


% Continue iteration, if the error is larger than the threshold, or
% iteration time is smaller than the predefined steps.
while e>epsilon && t<time   
    % new phi
    phi_tt=exp((psi(gama_t)-psi(sum(gama_t)))*ones(1,N)).*tempBeta;
    phi_tt=phi_tt./(ones(k,1)*sum(phi_tt+realmin,1));
    phi_tt=phi_tt.*filter;
    
    % new gamma
    gama_tt=alpha+sum(phi_tt,2);
    
    % error of the iteration
    e1=sum(sum(abs(phi_tt-phi_t)))/sum(sum(phi_t));
    e2=sum(abs(gama_tt-gama_t))/sum(gama_t);
    e=max(e1,e2);
    
    % update the variational parameters
    phi_t=phi_tt;
    gama_t=gama_tt;
    % disp(['t=',int2str(t),', e1,e2,e:',num2str(e1),',',num2str(e2),',',num2str(e)]);
    t=t+1;
end


