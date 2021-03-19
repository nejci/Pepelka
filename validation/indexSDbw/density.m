function dens = density(X,u,stdev)
% X - data points
% u - center of hyper-sphere
% stdev - radius of hyper-sphere

D = sqrt(sum(bsxfun(@minus,X,u).^2,2));
dens = sum(D <= stdev);

end