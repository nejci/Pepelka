function [p_corr, p_dist] = pearson_corr(x,y,mode,absolute)

if ~exist('mode','var') || isempty(mode)
	mode='';
end
if ~exist('absolute','var') || isempty(absolute)
	absolute=0;
end

% MATLAB kmeans
% pripravi podatke
% X = [x;y];
% [n,p]=size(X);
% Xn = X - repmat(mean(X,2),1,p);
% Xnorm = sqrt(sum(Xn.^2, 2));
% Xn = Xn ./ Xnorm(:,ones(1,p));
% 
% % ko zelis izracunati koeficient, samo zmnozi dva vzorca
% pc2 = X(1,:) * X(2,:)';

switch mode
	case 'uncentered'
		%uncentered Pearson correlation (without mean subtraction)
		p_corr = sum(x .* y) / ( sqrt(sum(x.^2)) * sqrt(sum(y.^2)) );

	case 'sqCentered'
		x_norm = x - mean(x);
		y_norm = y - mean(y);
		p_corr = (sum(x_norm .* y_norm) / ( sqrt(sum(x_norm.^2)) * sqrt(sum(y_norm.^2)) ))^2;
		
	case 'sqUncentered'
		p_corr = (sum(x .* y) / ( sqrt(sum(x.^2)) * sqrt(sum(y.^2)) ))^2;
		
	otherwise
		x_norm = x - mean(x);
		y_norm = y - mean(y);
		p_corr = sum(x_norm .* y_norm) / ( sqrt(sum(x_norm.^2)) * sqrt(sum(y_norm.^2)) );
end

if absolute
	p_corr = abs(p_corr);
end

% Distance from correlation
p_dist = max(1 - p_corr, 0);
%dist = (1-pc)/2;
%dist = 1 - pc^2;

end