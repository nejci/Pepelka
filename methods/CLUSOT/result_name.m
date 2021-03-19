function [r]=result_name(map, filename)
% [r]=result_name(map, filename)
%
% nn_result_name returns a filename for result of a som 
% training. Naming conventions described in
% $DIPLOM_ROOT/results/proben1/readme.txt are used.
%
% Inputs:
% map      - trained som
% filename - name of dataset, e.g. 'diabetes1'
% 
% $Id: result_name.m 421 2005-03-09 16:57:24Z dome $
% D. Brugger, 01 February 2005
% util/result_name.m

if(nargin == 0)
   test_result_name();
   return
end

dims = map.topol.msize;
[i,rl,tl]= map.trainhist.trainlen;
[init,alg0,alg1] = map.trainhist.algorithm;
[dump,alpha_type,dump] = map.trainhist.alpha_type;
shape = map.topol.shape;
alg='';
if(~strcmp(alg0,alg1))
  alg = [alg0 '-' alg1];
else
  alg = alg0;
end

if(strcmp(shape,'sheet'))
  r = [filename '_' init '_' alg '_' num2str(dims(1)) 'x' num2str(dims(2)) '-' map.topol.lattice ...
       '_' map.neigh '_' num2str(rl) '_' num2str(tl) '_' alpha_type ...
       '.mat'];
else
  r = [filename '_' init '_' alg '_' num2str(dims(1)) 'x' num2str(dims(2)) '-' map.topol.lattice ...
       '_' map.neigh '_' num2str(rl) '_' num2str(tl) '_' alpha_type ...
       '_' shape '.mat'];
end

function test_result_name()
D=rand(100,3);
map = som_make(D,'training',[50 50]);
% expected
er = 'diabetes_lininit_batch_8x6-hexa_gaussian_50_50_inv.mat';
r = result_name(map, 'diabetes');

if(~ strcmp(er, r))
  error('!!! test_result_name failed. Got %s, expected %s !!!', r, ...
	er)
end
fprintf('test_result_name succeded\n');
