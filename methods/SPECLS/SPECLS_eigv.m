function V = SPECLS_eigv(data,Knn)
% Precompute eigenvalues of data to speed-up computation of multiple
% repetitions.

if ~exist('Knn','var') || isempty(Knn)
    Knn = 7;
end

% centralize and scale the data
data = bsxfun(@minus,data,mean(data,1));
data = data/max(abs(data(:)));

% Build affinity matrix A
D = dist2(data,data);    % squared Euclidean distance
[~,A] = scale_dist_safe(D,Knn); % Locally scaled affinity matrix

% Zero out diagonal
N = size(data,1);
A(1:N+1:N^2) = 0;
% ZERO_DIAG = ~eye(size(data,1));
% A = A.*ZERO_DIAG;

% Compute the Laplacian
dd = 1./(sum(A,1)+eps);
dd = sqrt(dd);
DD = diag(dd);
L = DD*A*DD;

if sum(isinf(L(:))) > 0 || sum(isnan(L(:))) > 0
    fprintf(1,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
end

%%%%%%% Compute eigenvectors
%[V,~] = svd(L);
[~,~,V] = svdecon(L); % faster
