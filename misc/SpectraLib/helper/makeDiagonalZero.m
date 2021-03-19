function newA=makeDiagonalZero(A)
% newA=makeDiagonalZero(A)
%
% Makes the diagonal of the matrix zero. 
% Is useful when you are creating the similarity matrix someHow and want to make the diagonal entry zero. 

  newA=A-diag(diag(A));

