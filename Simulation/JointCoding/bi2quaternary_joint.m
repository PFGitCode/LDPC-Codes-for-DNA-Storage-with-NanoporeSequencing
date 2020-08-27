function encodedData = bi2quaternary_joint(biencodedData)
len = length(biencodedData)/2;
encodedData = zeros(1,len);
for i = 1:len
    if biencodedData(i) ==0 && biencodedData(i+len) == 0
        encodedData(i) = 0;
    elseif biencodedData(i) ==0 && biencodedData(i+len) == 1
        encodedData(i) = 1;
    elseif biencodedData(i) ==1 && biencodedData(i+len) == 0
        encodedData(i) = 2;
    elseif biencodedData(i) ==1 && biencodedData(i+len) == 1
        encodedData(i) = 3;
    end
end
end