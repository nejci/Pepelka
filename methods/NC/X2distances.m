function distances = X2distances(X,Sigma)
%Timothee Cour, 2004
n = size(X,1);
if exist('Sigma','var') && ~isempty(Sigma)
    X = X\sqrtm(Sigma);
end

temp = sum(X.*X,2);
temp = repmat(temp,1,n);
distances = -2*(X*X') + temp + temp';