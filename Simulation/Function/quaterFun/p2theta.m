function theta =  p2theta(usefulC,value,fixpoint)

len = size(usefulC,2);
if fixpoint == 1
    usefulC = usefulC(:,2:end);
elseif fixpoint == len
    usefulC = usefulC(:,1:end-1);
else
    usefulC = [usefulC(:,1:fixpoint-1),usefulC(:,fixpoint+1:end)];
end
len = len -1;
bb = usefulC;
multiveValue = ones(size(usefulC,1),1);
% a = zeros(size(bb));
% tic
% V = zeros(size(usefulC,1),len);
for i = 1:len
    multiveValue = multiveValue.*value(bb(:,i),i);
% value(bb,(1:5));
%       a(:,i) =value(bb(:,i),i);
%       V(:,i) =  value(bb(:,i),i);
%       multiveValue = prod(a,2);
end
% a = value(bb);
% toc
theta = (sum(multiveValue));
end