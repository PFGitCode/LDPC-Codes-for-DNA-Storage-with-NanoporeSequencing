function [data1,data2] = method3Soft(llr1,llr2, H1,H2,receivedSignal,mu,matrix)
muA = mu(1:2);
muC = mu(3:4);
muT = mu(5:6);
muG = mu(7:8);

mA = matrix(:,1:2);
mC = matrix(:,3:4);
mT = matrix(:,5:6);
mG = matrix(:,7:8);

[k, n] = size(H1);
data1 = zeros(n,1);
data2 = zeros(n,1);
p1 = zeros(n,1);
p2 = zeros(n,1);
Lq1 = zeros(k,n);
Lq2 = zeros(k,n);
Lr1 = zeros(k,n);
Lr2 = zeros(k,n);
LQ1 = zeros(n,1);
LQ2 = zeros(n,1);

pdfA = zeros(1,length(llr1));
pdfT = zeros(1,length(llr1));
pdfC = zeros(1,length(llr1));
pdfG = zeros(1,length(llr1));
Lu = zeros(1,length(llr1));
Lu1 = zeros(1,length(llr1));

for i = 1:length(llr1)
    p1(i) = llr1(i);%llr of p0 and p1
    p2(i) = llr2(i);%llr of p0 and p1
    Lq1(H1(:,i)==1,i) = p1(i);%initial variabl to check
    Lq2(H2(:,i)==1,i) = p2(i);%initial variabl to check
end

for i = 1:length(llr1)
    pdfA(i) = mvnpdf(receivedSignal(:,i)',muA,mA);
    pdfC(i) = mvnpdf(receivedSignal(:,i)',muC,mC);
    pdfT(i) = mvnpdf(receivedSignal(:,i)',muT,mT);
    pdfG(i) = mvnpdf(receivedSignal(:,i)',muG,mG);
    Lu(i) = pdfA(i) + pdfG(i);
    Lu1(i) = pdfC(i) + pdfT(i);
end

for iter = 1:50
    %---------------------message from checl to variable----------
    for j = 1:k
        nonZerosElement1 = find(H1(j,:)~=0);
        nonZerosElement2 = find(H2(j,:)~=0);
        sumj1 = prod(tanh(0.5*Lq1(j,nonZerosElement1)));
        sumj2 = prod(tanh(0.5*Lq2(j,nonZerosElement2)));
        temp1 = sumj1./(tanh(0.5.*Lq1(j,nonZerosElement1)));
        temp2 = sumj2./(tanh(0.5.*Lq2(j,nonZerosElement2)));
        for t = 1:length(temp1)
            if abs(temp1(t)) == 1
                Lr1(j,nonZerosElement1(t)) = 2*19.07*temp1(t);
            else
                Lr1(j,nonZerosElement1(t)) = 2*atanh(temp1(t));
            end
        end
        for t = 1:length(temp2)
            if abs(temp2(t)) == 1
                Lr2(j,nonZerosElement2(t)) = 2*19.07*temp2(t);
            else
                Lr2(j,nonZerosElement2(t)) = 2*atanh(temp2(t));
            end
        end
    end
    %---------get the codeword and check if it is right(H*data=0)-----------
    for i = 1:n
        sumi1 = sum(Lr1(:,i));
        sumi2 = sum(Lr2(:,i));
        LQ1(i) = p1(i) + sumi1;
        LQ2(i) = p2(i) + sumi2;
        if (pdfT(i) > pdfA(i) && pdfT(i) > pdfG(i) && receivedSignal(2,i) < muG(2) &&  receivedSignal(1,i) < muA(1))
            temp21 = tanh(0.5*(sumi2))*tanh(0.5*log((Lu(i))/(Lu1(i))));
            if abs(temp21) == 1
                sumi21 = 2*19.07*temp21;
            else
                sumi21 = 2*atanh(temp21);
            end
            temp12 = tanh(0.5*(sumi1))*tanh(0.5*log((Lu(i))/(Lu1(i))));
            if abs(temp21) == 1
                sumi12 = 2*19.07*temp12;
            else
                sumi12 = 2*atanh(temp12);
            end
            LQ1(i) = LQ1(i)+sumi21;
            LQ2(i) = LQ2(i)+sumi12;
        end
        data1(i) = LQ1(i) < 0;
        data2(i) = LQ2(i) < 0;
    end
    re1 = H1*data1;
    re2 = H2*data2;
    ifend1 = 1;
    ifend2 = 1;
    for i = 1:length(re1)
        if(mod(re1(i),2) == 1)
            ifend1 = 0;
        end
        if(mod(re2(i),2) == 1)
            ifend2 = 0;
        end
    end
    if(ifend1&&ifend2)
        break;
    end
    %-------------variable to check message------------------------------------
    for i = 1:n
        nonZerosElementi1 = find(H1(:,i)~=0);
        nonZerosElementi2 = find(H2(:,i)~=0);
        sumi1 = sum(Lr1(:,i));
        sumi2 = sum(Lr2(:,i));
        Lq1(nonZerosElementi1,i) = p1(i) + sumi1 - Lr1(nonZerosElementi1,i);
        Lq2(nonZerosElementi2,i) = p2(i) + sumi2 - Lr2(nonZerosElementi2,i);
        if (pdfT(i) > pdfA(i) && pdfT(i) > pdfG(i) && receivedSignal(2,i) < muG(2) &&  receivedSignal(1,i) < muA(1))
            temp21 = tanh(0.5*(sumi2))*tanh(0.5*log((Lu(i))/(Lu1(i))));
            if abs(temp21) == 1
                sumi21 = 2*19.07*temp21;
            else
                sumi21 = 2*atanh(temp21);
            end
            %             sumi21 = 2*atanh(temp21);
            temp12 = tanh(0.5*(sumi1))*tanh(0.5*log((Lu(i))/(Lu1(i))));
            if abs(temp21) == 1
                sumi12 = 2*19.07*temp12;
            else
                sumi12 = 2*atanh(temp12);
            end
            Lq1(nonZerosElementi1,i) = Lq1(nonZerosElementi1,i)+sumi21;
            Lq2(nonZerosElementi2,i) = Lq2(nonZerosElementi2,i)+sumi12;
        end
    end
end
end
