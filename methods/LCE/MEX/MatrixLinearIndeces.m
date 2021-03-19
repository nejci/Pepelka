M = 4;
numPairs = M*(M-1)/2;

% With diagonal
for i = 0:numPairs-1+M;
    
    r = (-2*M - 1 + sqrt( (4*M*(M+1) - 8*i - 7) )) / -2;
    if r == fix(r)
        r = r-1;
    end
    r = fix(r);
    
    c = i - M * r + r*(r+1) / 2;
    c = fix(c);
    fprintf(1,'i: %d, r: %d, c: %d\n',i,r,c);
end

% Without diagonal - upper, one based indexing
for i = 1:numPairs;
    c = round(floor(-.5 + .5 * sqrt(1 + 8 * (i-1))) + 2);
    r = round(c * (3 - c) / 2 + (i-1));
    fprintf(1,'i: %d, r: %d, c: %d\n',i,r,c);
end

% Without diagonal - upper, zero based indexing
for i = 0:numPairs-1;
    c = round(floor(-.5 + .5 * sqrt(1 + 8 * i)) + 2);
    r = round(c * (3 - c) / 2 + i);
    c = c-1;
    r = r-1;
    fprintf(1,'i: %d, r: %d, c: %d\n',i,r,c);
end

