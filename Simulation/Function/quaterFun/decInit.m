function [nonZerosElements,s] = decInit(H,maxlen)
[n,k] = size(H');
% len = 6;
% totalCases = 4^maxlen;
% UsefulcodeLen = 4^(maxlen-2);
% quarterArrayString = 0:1:totalCases-1;
% quarterArrayString = dec2base(quarterArrayString,4) - '0';
% valueSize = [4,maxlen];
% index = (1:maxlen);
s = cell(k,maxlen);
nonZerosElements = zeros(k,maxlen);
for i  = 1:k
    hi = H(i,:);
    nonZerosElement = find(hi~=0);
    nonZerosElements(i,1:length(nonZerosElement)) = nonZerosElement;
    h = (hi(nonZerosElement));
    %     sydromValue =  gfMulti(quarterArrayString,h);
    %     if i == 76
    %         1
    %     end
    %     if length(h) > 7
    %         1
    %     end
    len = length(h);
    totalCases = 4^len;
    UsefulcodeLen = 4^(len-2);
    quarterArrayString = 0:1:totalCases-1;
    quarterArrayString = dec2base(quarterArrayString,4) - '0';
    sydromValue = gf(quarterArrayString,2)*gf(h',2);
    usefulC  = quarterArrayString(sydromValue == 0,:);
    %     usefulCodeWord0 = usefulC(usefulC(:,j)==0,:);
    %     usefulCs(i).c = usefulC;
    valueSize = [4,len];
    index = (1:len);
    for j = 1:len
        usefulCodeWord0 = usefulC(usefulC(:,j)==0,:)+1;
        %         usefulCs(i,j).usefulCodeWord0 = usefulCodeWord0;
        usefulCodeWord1 = usefulC(usefulC(:,j)==1,:)+1;
        %         usefulCs(i,j).usefulCodeWord1 = usefulCodeWord1;
        usefulCodeWord2 = usefulC(usefulC(:,j)==2,:)+1;
        %         usefulCs(i,j).usefulCodeWord2 = usefulCodeWord2;
        usefulCodeWord3 = usefulC(usefulC(:,j)==3,:)+1;
        %         usefulCs(i,j).usefulCodeWord3 = usefulCodeWord3;
        repMat = repmat(index,UsefulcodeLen*4,1);
        s{i,j} = uint8(sub2ind(valueSize,[usefulCodeWord0;usefulCodeWord1;usefulCodeWord2;usefulCodeWord3],repMat));
        ss = s{i,j};
        if max(ss(:))>32
            1;
        end
    end
end
end