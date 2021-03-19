% This script obtains the pdf of the random variable resulting from the sum of
% indepedent weigthed random varaibles. The resulting ppdf is obtained
% using the characteristic function.
%
% Henry Carrillo, University of Zaragoza, Spain
% http://www.hcarrillo.com/
% -------------------------------------------------------------------------
function y = betaChfSum(res, AB)
% Set support
supportPdf = size(AB,1);
x = 0:res:supportPdf;

% Individual Beta pdfs
fBeta = zeros(supportPdf,length(x));
for i=1:supportPdf
    fBeta(i,:) = betapdf(x, AB(i,1), AB(i,2));
end

% Individual Characteristics functions and sum of independent RV
fsumrv = 1/supportPdf;
for i=1:supportPdf
    fChf = fft(fBeta(i,:));
    fsumrv = fsumrv.*fChf;
end

% Compute pdf
y = ifft(fsumrv);
% Reduce to [0..supportPdf] support
y = y(1:length(x));
% Normalize (so that all values sum to 1/res)
y = y / (sum(y) * res);
end
