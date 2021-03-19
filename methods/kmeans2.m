function labels=kmeans2(data,K,params)
%
% testKMEANS je funkcija, ki vhodne podatke v data razvrsti v K gru� po algoritmu K-means 
%
% VHODI:
%   data - matrika vhodnih podatkov (velikosti Nxd, kjer je N �tevilo vzorcev in d njihova razse�nost)
%   K - �tevilo gru�, v katere razvr��amo vzorce
%   maxIter - parameter, ki dolo�a, koliko iteracij bo algoritem naredil.{100}
%   fLog - fid dnevni�ke datoteke
%
% IZHODI:
%   labels - vektor oznak razvr��enih vzorcev; uporabljena so naravna �tevila.
%
%
% Avtor: Nejc Ilc, Fakulteta za ra�unalni�tvo in informatiko, 2008
% Zadnja sprememba: 27. november 2008

maxIter=params.KM_maxIter;
replicates=params.KM_nRuns;

    
tic;
[labels,c,sumd]=kmeans(data,K, 'distance','sqEuclidean','start','sample','replicates', replicates, 'maxiter',maxIter, 'emptyaction','drop','display','notify');
time=toc;

%moreInfo.time=time;