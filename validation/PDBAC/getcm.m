function [confus,numcorrect,precision,recall,F] = getcm (actual,pred,classes)

% GETCM : gets confusion matrices, precision, recall, and F scores
% [confus,numcorrect,precision,recall,F] = getcm (actual,pred,[classes])
%
% actual is a N-element vector representing the actual classes
% pred is a N-element vector representing the predicted classes
% classes is a vector with the numbers of the classes (by default, it is 1:k, where k is the
%    largest integer to appear in actual or pred.
%
% dinoj@cs.uchicago.edu , Apr 2005, modified July 2005
% [confus,numcorrect,precision,recall,F] = getcm ([1 2 3 2 2 3 1 2 3],[1 2 3 3 2 1 1 2 3],1:3)

if size(actual,2) == 1
    actual = actual';
end
if size(pred,2) == 1
    pred = pred';
end
if nargin < 3
    classes = 1:max(max(actual),max(pred));
end
nC = classes(end);
numcorrect = sum(actual==pred);
confus = zeros(nC);
for i=1:nC
    d = actual==i; % d ones where points are with class a
    confus(i,:) = sum(bsxfun(@eq,pred(d)',classes),1);
end

% Slower for large nC and N
% for i=1:nC
%     a = classes(i);
%     d = target==a;     % d ones where points are with class a
%     for j=1:nC
%         confus(i,j) = length(find(pred(d)==classes(j)));
%     end
% end



precision=zeros(nC,1);
recall=zeros(nC,1);
F=zeros(nC,1);
for i=1:nC
    S = sum(confus(i,:));
    if nargout>=4
        if S
            recall(i) = confus(i,i) / sum(confus(i,:));
        else
            recall(i) = 0;
        end
    end
    S =  sum(confus(:,i));
    if nargout>=3
        if S
            precision(i) = confus(i,i) / S;
        else
            precision(i) = 0;
        end
    end
    if nargout>=5
        if (precision(i)+recall(i))
            F(i) = 2 * (precision(i)*recall(i)) / (precision(i)+recall(i));
        else
            F(i) = 0;
        end
    end
end

