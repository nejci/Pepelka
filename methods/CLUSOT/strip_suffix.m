function [r]=strip_suffix(filename)
% [r]=strip_suffix(filename)
%
% strip_sufffix removes suffix, e.g. '.mat'
% from filename.
%
% Input:
% filename
%
% Output:
% r - filename w/o suffix
% 
% $Id: strip_suffix.m 374 2005-02-04 11:44:06Z dome $
% D. Brugger, 04 February 2005
% util/strip_suffix.m

tmp=regexp(filename,'\.','start');
r = filename(1:tmp(end)-1);
