function y = loggausspdf(X, mu, sigma)
% log pdf of Gaussian
% Written by Mo Chen (mochen@ie.cuhk.edu.hk). March 2009.
[d,k] = size(mu);

if size(sigma,1)==d && size(sigma,2)==d && k==1
    X = bsxfun(@minus,X,mu);
    [R,p]= chol(sigma);
    if p ~= 0
        error('ERROR: sigma is not SPD (not positive definite).');
    end
    q = sum((R'\X).^2,1);  % quadratic term (M distance)
    c = d*log(2*pi)+2*sum(log(diag(R)));   % normalization constant
    y = -(c+q)/2;
elseif size(sigma,1)==d && size(sigma,2)==k
    lambda = 1./sigma;
    q = bsxfun(@plus,X'.^2*lambda-2*X'*(mu.*lambda),sum((mu.^2).*lambda,1)); % M distance
    c = (d*log(2*pi)+sum(log(sigma),1))/(-2); % normalization constant
    y = bsxfun(@plus,q/(-2),c);
elseif size(sigma,1)==1 && size(sigma,2)==k
    X2 = repmat(sum(X.^2,1)',1,k);
    D = bsxfun(@plus,X2-2*X'*mu,sum(mu.^2,1));
    q = bsxfun(@times,D,1./sigma);  % M distance
    c = d*log(2*pi*sigma)/(-2);          % normalization constant
    y = bsxfun(@plus,q/(-2),c);
else
    error('Parameters mismatched.');
end
