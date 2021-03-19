function [sigma]=estimate_sigma(data)
%data must be matrix of shape: patterns x dimensions
%function returns estimated sigma for Parzen windowing using Silverman's
%rule of a thumb, as described at:
%http://www.phys.uit.no/~robertj/ROBERT_NOBS/PDF/04244695.pdf
%Author: Nejc Ilc

[N,d]=size(data);

%calculation of sample covariance matrix of data
C=cov(data);

%sum of diagonal elements of covariance matrix
S=trace(C);

%deviation of samples; bolje je, da ni korena - empirièno
std_data=sqrt((1/d)*S);

sigma=std_data*(4/((2*d+1)*N))^(1/(d+4));