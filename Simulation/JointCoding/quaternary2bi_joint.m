function biData = quaternary2bi_joint(quaterData)
biData1 = zeros(size(quaterData));
biData2 = zeros(size(quaterData));
for i = 1:length(quaterData)
    switch quaterData(i)
        case 0
            biData1(i) = 0;
            biData2(i) = 0;
        case 1
            biData1(i) = 0;
            biData2(i) = 1;
        case 2
            biData1(i) = 1;
            biData2(i) = 0;
        case 3
            biData1(i) = 1;
            biData2(i) = 1;
    end
end
biData = [biData1,biData2];
end