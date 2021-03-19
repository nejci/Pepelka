% function s = checks(s)
% 
% DESCRIPTION
%   Checks a similarity matrix for validity and fixes
%   detected problems
%
% Copyright (c) 1998-2002 by Alexander Strehl


function s = checks(s)
if ~isempty(s),
   if sum(sum(~isreal(s))),
      warning('checks: complex similarities found');
      s = real(s);
      warning('checks: using real component');
   end;
   if size(s,1)~=size(s,2)
      warning('checks: s is not square');
      s = s(1:min(size(s)),1:min(size(s)));
      warning('checks: using indrawn square');
   end;
   mas = max(max(s));
   mis = min(min(s));
   if (mas>1)||(mis<0),
      warning(['checks: similarity more than 1 or less than 0 detected: values ' num2str(mis) ' to ' num2str(mas) ]); 
      s(s<0)=0;
      s(s>1)=1;
      warning('checks: bounded');
   end;
   if sum(sum(isinf(s)|isnan(s))),
      warning('checks: non-finite similarity detected !!! (serious)'); 
      if 0,
         s(find(isinf(s)|isnan(s))) = 0; % hangs the computer - no idea why...
      else
        [a,b] = find(isfinite(s));
         c = isfinite(s);
         s = sparse(a,b,s(c)); 
      end;
      warning('checks: made zero !!! (serious)');
   end;
   if (sum(sum(s~=s'))>0),
      warning('checks: s is not symmetric');
      s = (s+s')./2;     
      warning('checks: symmetrised');
   end;
   if (size(s,1)==size(s,2)),
      if sum(diag(s)~=1)
         warning('checks: self-similarity s(i,i) not always 1');
         for i=1:size(s,1),
            s(i,i) = 1;
         end;
         warning('checks: diagonal made 1');
      end;
   end;
else
   warning('checks: empty similarity matrix');	
end;
