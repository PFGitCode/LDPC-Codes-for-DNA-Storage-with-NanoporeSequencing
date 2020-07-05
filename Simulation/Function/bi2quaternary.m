function encodedData = bi2quaternary(encodedData1,encodedData2)
len = length(encodedData1);
encodedData = zeros(1,len);
for i = 1:len
    if encodedData1(i) ==0 && encodedData2(i) == 0
        encodedData(i) = 0;
    elseif encodedData1(i) ==0 && encodedData2(i) == 1
        encodedData(i) = 1;
    elseif encodedData1(i) ==1 && encodedData2(i) == 0
        encodedData(i) = 2;
    elseif encodedData1(i) ==1 && encodedData2(i) == 1
        encodedData(i) = 3;
    end
end
end