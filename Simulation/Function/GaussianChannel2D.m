function [dataTime, dataCurrent] = GaussianChannel2D(data, mu,matrix)
muA = mu(1:2);
muC = mu(3:4);
muT = mu(5:6);
muG = mu(7:8);

matrixA = matrix(:,1:2);
matrixC = matrix(:,3:4);
matrixT = matrix(:,5:6);
matrixG = matrix(:,7:8);

dataTime = zeros(size(data));
dataCurrent = zeros(size(data));

for i = 1:length(data)
    if data(i) == 0
        output = mvnrnd(muA,matrixA);
        dataTime(i) = output(1);
        dataCurrent(i) = output(2);
    elseif data(i) == 1
        output = mvnrnd(muC,matrixC);
        dataTime(i) = output(1);
        dataCurrent(i) = output(2);
    elseif data(i) == 2
        output = mvnrnd(muT,matrixT);
        dataTime(i) = output(1);
        dataCurrent(i) = output(2);
    elseif data(i) == 3
        output = mvnrnd(muG,matrixG);
        dataTime(i) = output(1);
        dataCurrent(i) = output(2);
%     elseif data(i) == -1
%         dataTime(i) = -100;
%         dataCurrent(i) = -100;
    end
end
end