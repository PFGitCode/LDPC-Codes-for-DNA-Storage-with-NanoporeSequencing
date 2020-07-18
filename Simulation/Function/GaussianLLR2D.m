function llr = GaussianLLR2D(inputdataTime,inputdataCurrent,mu,matrix)
    n = length(inputdataTime);
    llr = zeros(2*n,1);
%     llr2 = zeros(n,1);
    
    muA = mu(1:2);
    muC = mu(3:4);
    muT = mu(5:6);
    muG = mu(7:8);

    matrixA = matrix(:,1:2);
    matrixC = matrix(:,3:4);
    matrixT = matrix(:,5:6);
    matrixG = matrix(:,7:8);
    
    for i = 1:n
        yTime = inputdataTime(i);
        yCurrent = inputdataCurrent(i);
        a = mvnpdf([yTime,yCurrent],muA,matrixA);
        c = mvnpdf([yTime,yCurrent],muC,matrixC);
        g = mvnpdf([yTime,yCurrent],muG,matrixG);
        t = mvnpdf([yTime,yCurrent],muT,matrixT);
        llr(i) = log((c+a)/(g+t));
        llr(i + n)= log((t+a)/(g+c));
        if isnan(llr(i))
            1;
        end
    end
end