% function dout = mapdense(din,column)
%
% copyright (c) 1998-2002 by Alexander Strehl

function dout = mapdense(din,column)

if ~exist('column'),
  column = 1;
end;

dout = mapasc(sortrows(din,column),column);
