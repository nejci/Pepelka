function c = pplk_nchoosek2(n,k)
% Using Matlab code of function nchoosek but without overhead.
% If n is a vector, c is also a vector, where c(i) = nchoosek(n(i),k)
n = n(:);
N = length(n);
tolerance = 1e15;

c = zeros(N,1);
maskNan = k > n;
c(maskNan) = NaN; % K out of range

k = ones(N,1).*k;
k(maskNan) = NaN;

mask = k > (n./2);
k(mask) = n(mask)-k(mask);

mask = k <= 1;
c(mask) = n(mask).^k(mask);

ind = find(~mask & ~maskNan);

for i = ind'
    nums = (n(i)-k(i)+1):n(i);
    dens = 1:k(i);
    nums = nums./dens;
    c(i) = round(prod(nums));
end

if any(c > tolerance)
    warning(message('MATLAB:nchoosek:LargeCoefficient', sprintf( '%e', tolerance ), log10( tolerance )));
end
