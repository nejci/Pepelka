function [new_s]=remove_empty(s)
% [new_s]=remove_empty(s)
%
% remove_empty removes all emtpy matrices in cell array s
%
% Input:
% s - (cell array of matrices)
%
% Output:
% new_s - (cell array of matrices) with all emtpy entries removed,
%           e.g. all entries with index k where s{k}=[]
% 
% $Id: remove_empty.m 584 2005-04-18 14:09:59Z dome $
% D. Brugger, 18 April 2005
% util/remove_empty.m

if(nargin == 0)
  test_remove_empty();
  return;
end

new_s={};
pos1=1; pos2=1;
for k=1:size(s,1)
  for l=1:size(s,2)
    if(~isempty(s{k,l}))
      new_s{pos1,pos2}=s{k,l};
      pos2 = pos2 + 1;
    end
  end
  pos2 = 1;
  pos1 = pos1 + 1;
end

function test_remove_empty()
% test case #1
s{1}=[1 2 3];
s{2}=[];
s{3}=[2];
enew_s{1}=[1 2 3];
enew_s{2}=[2];
new_s=remove_empty(s)
check_matching(enew_s,new_s);
% test case #2
clear s;
s{1}=[];
s{2}=[];
enew_s={};
new_s=remove_empty(s)
check_matching(enew_s,new_s);
% test case #3
clear s;
s{1}=[1 2 3; 4 5 6];
s{2}=[0 0 1; 1 2 3; 0 0 0];
s{3}=[];
clear enew_s;
enew_s{1}=[1 2 3; 4 5 6];
enew_s{2}=[0 0 1; 1 2 3; 0 0 0];
new_s=remove_empty(s)
check_matching(enew_s,new_s);
fprintf('test_remove_empty succeded\n');
