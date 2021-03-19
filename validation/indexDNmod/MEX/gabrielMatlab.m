function d = gabrielMatlab(distMat,tol)

n = size(distMat,1);
d = triu(distMat,1);
% Cycle thru all possible pairs of points
for i = 1:(n-1)
    for j = (i+1):n
        for k = 1:n
            if (k~=i && k~=j)
                if (distMat(i,k)+distMat(j,k)-tol <= distMat(i,j))
                    d(i,j) = 0;
                    break;
                end
            end
        end
    end
end