function decodeWord = quaternaryDecoderSoft(receivedData,H,mu,matrix,nonZerosElements,ss) % P = [Pac,Pat,Pag,Pct,Pcg,Ptg]
muA = mu(1:2);
muC = mu(3:4);
muT = mu(5:6);
muG = mu(7:8);

matrixA = matrix(:,1:2);
matrixC = matrix(:,3:4);
matrixT = matrix(:,5:6);
matrixG = matrix(:,7:8);

[n,k] = size(H');
decodeWord = zeros(1,length(receivedData));
P = zeros(4,n,1);
q = zeros(4,k,n);
theta0 = zeros(k,n);
theta1 = zeros(k,n);
theta2 = zeros(k,n);
theta3 = zeros(k,n);
for i = 1 : n
    yTime = receivedData(1,i);
    yCurrent = receivedData(2,i);
    P(1,i) = mvnpdf([yTime,yCurrent],muA,matrixA);
    P(2,i) = mvnpdf([yTime,yCurrent],muC,matrixC);
    P(3,i) = mvnpdf([yTime,yCurrent],muT,matrixT);
    P(4,i) = mvnpdf([yTime,yCurrent],muG,matrixG);

    q(1,:,i) = P(1,i);
    q(2,:,i) = P(2,i);
    q(3,:,i) = P(3,i);
    q(4,:,i) = P(4,i);
end
for iter = 1:50
    for i  = 1:k
        nonZerosElement = nonZerosElements(i,:);
        nonZerosElement = nonZerosElement(nonZerosElement~=0);
        len = length(nonZerosElement);
        UsefulcodeLen = 4^(len-2);
        valueq = q(:,:,nonZerosElement);
        valulei = reshape(valueq(:,i,:),[4,len]);
        for j = 1:len
            mul = prod(valulei(double(ss{i,j})),2);
            theta0(i,nonZerosElement(j)) = sum(mul(1:UsefulcodeLen)/valulei(1,j));
            theta1(i,nonZerosElement(j)) = sum(mul(UsefulcodeLen+1:UsefulcodeLen*2)/valulei(2,j));
            theta2(i,nonZerosElement(j)) = sum(mul(UsefulcodeLen*2+1:UsefulcodeLen*3)/valulei(3,j));
            theta3(i,nonZerosElement(j)) = sum(mul(UsefulcodeLen*3+1:UsefulcodeLen*4)/valulei(4,j));
        end
        clearvars valulei valueq mul;
    end
    finalValue = zeros(1,4);
    for j = 1:n
        finalValue(1) = decisionMake(theta0(:,j),P(1,j));
        finalValue(2) = decisionMake(theta1(:,j),P(2,j));
        finalValue(3) = decisionMake(theta2(:,j),P(3,j));
        finalValue(4) = decisionMake(theta3(:,j),P(4,j));
        [value,index] = max(finalValue);
        decodeWord(j) = index - 1;
    end
    clearvars finalValue;
    
    hf = gf(H,2);
    decodef = gf(decodeWord,2);
    re = decodef*hf';
    
    re = re.x;
    ifend = 1;
    for i = 1:k
        if(re(i) ~= 0)
            ifend = 0;
        end
    end
    if ifend == 1
        break;
    end
    
    for j = 1:n
        hj =  H(:,j);
        nonZerosElement = find(hj~=0);
        for i = 1:length(nonZerosElement)
            q(1,nonZerosElement(i),j) = theta2p(nonZerosElement(i),theta0(:,j),P(1,j));
            q(2,nonZerosElement(i),j) = theta2p(nonZerosElement(i),theta1(:,j),P(2,j));
            q(3,nonZerosElement(i),j) = theta2p(nonZerosElement(i),theta2(:,j),P(3,j));
            q(4,nonZerosElement(i),j) = theta2p(nonZerosElement(i),theta3(:,j),P(4,j));
            qSum = q(1,nonZerosElement(i),j)+q(2,nonZerosElement(i),j)+q(3,nonZerosElement(i),j)+q(4,nonZerosElement(i),j);
            q(1,nonZerosElement(i),j) = q(1,nonZerosElement(i),j)/qSum;
            q(2,nonZerosElement(i),j) = q(2,nonZerosElement(i),j)/qSum;
            q(3,nonZerosElement(i),j) = q(3,nonZerosElement(i),j)/qSum;
            q(4,nonZerosElement(i),j) = q(4,nonZerosElement(i),j)/qSum;
        end
    end
end
% P0 = signma;
end


