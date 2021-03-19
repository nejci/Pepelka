function labels=kmeans2(data,K,params)
%
% testKMEANS je funkcija, ki vhodne podatke v data razvrsti v K gruè po algoritmu K-means 
%
% VHODI:
%   data - matrika vhodnih podatkov (velikosti Nxd, kjer je N število vzorcev in d njihova razsežnost)
%   K - število gruè, v katere razvršèamo vzorce
%   maxIter - parameter, ki doloèa, koliko iteracij bo algoritem naredil.{100}
%   fLog - fid dnevniške datoteke
%
% IZHODI:
%   labels - vektor oznak razvršèenih vzorcev; uporabljena so naravna števila.
%
%
% Avtor: Nejc Ilc, Fakulteta za raèunalništvo in informatiko, 2008
% Zadnja sprememba: 27. november 2008

maxIter=params.KM_maxIter;
replicates=params.KM_nRuns;

    
tic;
[labels,c,sumd]=kmeans(data,K, 'distance','sqEuclidean','start','sample','replicates', replicates, 'maxiter',maxIter, 'emptyaction','drop','display','notify');
time=toc;

%moreInfo.time=time;