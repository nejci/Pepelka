function [labels, moreInfo]=pplk_clustererCS(data,K,params)
%
% testCSC je funkcija, ki vhodne podatke v data razvrsti v K gru� po algoritmu Cauchy-Schwartz divergence IC clustering 
%
% VHODI:
%   data - matrika vhodnih podatkov (velikosti Nxd, kjer je N �tevilo vzorcev in d njihova razse�nost)
%   K - �tevilo gru�, v katere razvr��amo vzorce
%   sigma - �irina Gaussove jedrne funkcije, �e ni podana, se oceni s
%   Silvermanovim pravilom
%   Kin - za�etno �tevilo skupin
%   Nin - za�etno �tevilo sosedov
%
% IZHODI:
%   labels - vektor oznak razvr��enih vzorcev; uporabljena so naravna �tevila.
%

sigma = [];
Kin = 20;
Nin = 10;

if exist('params','var') && isstruct(params)
    if isfield(params,'CS_sigma')
        sigma = params.CS_sigma;
    end
    if isfield(params,'CS_Kin')
        Kin = params.CS_Kin;    
    end
    if isfield(params,'CS_Nin')
        Nin = params.CS_Nin;    
    end
end

N = size(data,1);

%check Kin and Nin in relation to N and K
%Final number of clusters should not be greater than Initial number of
%clusters.
Kin = max(Kin,K);

% If Kin*Nin exceeds N, adjust Nin
%Rule of thumb: Kin*Nin = 2/3 N
if Kin*Nin > N
	% disp('!!! CS: Violated rule: Kin*Nin < N !!!')
	% Kin=floor(sqrt(N*0.6667));	    
    Nin = ceil(0.6667*N/Kin);    
end

oldPath=chdir('CS');

%Silverman's rule of thumb for sigma 
if(isempty(sigma))
   sigma=estimate_sigma(data);
   fprintf('Sigma estimation: %f\n', sigma);
end



ticID = tic();
labels = CSC(data,sigma,K,Kin,Nin);
time=toc(ticID);

labels=labels';
moreInfo.time=time;
moreInfo.options.sigma=sigma;
moreInfo.options.Kin=Kin;
moreInfo.options.Nin=Nin;

chdir(oldPath);
