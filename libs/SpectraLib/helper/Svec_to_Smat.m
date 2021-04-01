function [S]=Svec_to_Smat(Svec,dim)

S = zeros(dim,dim);

loc = 1;
for j=1:dim-1
  S(j,j+1:dim) = Svec(loc:loc+dim-j-1)';
  loc = loc + dim-j;
end
S = S + S';
    
end