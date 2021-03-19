function value = fieldvalue(x, field, default)
% Written by Mo Chen (mochen@ie.cuhk.edu.hk). March 2009.
if isfield(x,field)
    value = x.(field);
else
    value = default;
end