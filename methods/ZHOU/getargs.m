function [errmsg,varargout] = getargs(pnames,defaults,varargin)
%
% resolve the varargin arguments for km.m and clusterer_ensemble.m
%
% ATTN: This package is free for academic usage. The code was developed by Mr. W. Tang (wtang314@yahoo.com). You can run
% it at your own risk. For other purposes, please contact Prof. Zhi-Hua Zhou (zhouzh@nju.edu.cn)
%
% ATTN2: This package was developed by Mr. W. Tang (wtang314@yahoo.com). For any problem concerning the code, please feel
% free to contact Mr. Tang.
%

errmsg = '';
nparams = length(pnames);
varargout = defaults;
unrecog = {};
nargs = length(varargin);

if mod(nargs,2)~=0
    errmsg = sprintf('wrong number of arguments');
else
    for j=1:2:nargs
        pname = varargin{j};
        if ~ischar(pname)
            emsg = sprintf('parameter name must be in text');
            break;
        end
        i = strmatch(lower(pname),pnames);
        if isempty(i)
            if nargout > nparams+1
                unrecog((end+1):(end+2)) = {varargin{j} varargin{j+1}};
            else
                errmsg = sprintf('invalid parameter name:  %s',pname);
                break;
            end
        elseif length(i)>1
            errmsg = sprintf('ambiguous parameter name:  %s',pname);
            break;
        else
            varargout{i} = varargin{j+1};
        end
    end
end
varargout{nparams+1} = unrecog;
