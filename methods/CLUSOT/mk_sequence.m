function [non_seq]=mk_sequence(non_seq)
% [seq]=mk_sequence(non_seq)
%
% mk_sequence forms a sequence of numbers from a non-sequence,
% e.g. [1 1 4 1 3 1] is transformed to [1 1 3 1 2 1].
%
% Input:
% non_seq - (matrix) non-sequence of numbers
%
% Output:
% seq - (vector)
% 
% $Id: mk_sequence.m 541 2005-04-06 12:10:17Z dome $
% D. Brugger, 24 March 2005
% util/mk_sequence.m

if(nargin == 0)
  test_mk_sequence();
  return;
end

k=1;
while(k < max(max(non_seq)))
  while(isempty(find(non_seq == k)))
    ind = find(non_seq > k);
    non_seq(ind) = non_seq(ind) - 1;
  end
  k = k + 1;
end

function test_mk_sequence()
% Test case #1 - vector
non_seq=[1 1 1 4 1 1 1 1 1];
eseq=   [1 1 1 2 1 1 1 1 1];
seq=mk_sequence(non_seq)
check_equal(eseq,seq,'eseq','seq');

% Test case #2 - vector
non_seq=[1 1 1 4 1 3 3 1 1];
eseq=   [1 1 1 3 1 2 2 1 1];
seq=mk_sequence(non_seq)
check_equal(eseq,seq,'eseq','seq');

% Test case #3 - matrix
non_seq=[0 0 1 0; 0 3 3 0];
eseq=   [0 0 1 0; 0 2 2 0];
seq=mk_sequence(non_seq)
check_equal(eseq,seq,'eseq','seq');

% Test case #4 - matrix
non_seq=[0 0 0 0; 0 0 1 0; 0 5 0 4];
eseq=   [0 0 0 0; 0 0 1 0; 0 3 0 2];
seq=mk_sequence(non_seq)
check_equal(eseq,seq,'eseq','seq');

fprintf('test_mk_sequence succeded\n');
