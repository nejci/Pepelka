function dr=clusot_root()
% dr=clusot_root()
%
% clusot_root determines the root directory of
% the toolbox
%
% Output:
% dr (string)
% 
% $Id: clusot_root.m 619 2005-04-27 07:35:11Z dome $
% D. Brugger, 19 March 2005
% util/clusot_root.m
dr = strrep(which('clusot_root'),[filesep 'clusot_root.m'],'');

