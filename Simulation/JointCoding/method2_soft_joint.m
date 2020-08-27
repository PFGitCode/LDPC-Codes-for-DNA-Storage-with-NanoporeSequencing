%I use LLR and SPA to do the LDPC decode
%reference
%1.Channel code classical and Modern by William E.Ryan and Shu Lin
%Chapter5.4
%2. Maltab Document https://www.mathworks.com/help/comm/ref/ldpcdecoder.html
function data = method2_soft_joint(llr, H,receivedSignal,mu,matrix)
muA = mu(1:2);
muC = mu(3:4);
muT = mu(5:6);
muG = mu(7:8);

mA = matrix(:,1:2);
mC = matrix(:,3:4);
mT = matrix(:,5:6);
mG = matrix(:,7:8);

[k, n] = size(H);
data = zeros(n,1);
p = zeros(n,1);
Lq = zeros(k,n);
Lr = zeros(k,n);
LQ = zeros(n,1);

pdfA = zeros(1,length(llr)/2);
pdfT = zeros(1,length(llr)/2);
pdfC = zeros(1,length(llr)/2);
pdfG = zeros(1,length(llr)/2);

for i = 1:n
    p(i) = llr(i);%llr of p0 and p1
    Lq(H(:,i)==1,i) = p(i);%initial variabl to check
end

for i = 1:length(llr)/2
    pdfA(i) = mvnpdf(receivedSignal(:,i)',muA,mA);
    pdfC(i) = mvnpdf(receivedSignal(:,i)',muC,mC);
    pdfT(i) = mvnpdf(receivedSignal(:,i)',muT,mT);
    pdfG(i) = mvnpdf(receivedSignal(:,i)',muG,mG);
end

for iter = 1:50
    %--------------------message from checl to variable---------
    for j = 1:k
        nonZerosElement = find(H(j,:)~=0);
        sumj = prod(tanh(0.5*Lq(j,nonZerosElement)));
        temp = sumj./(tanh(0.5.*Lq(j,nonZerosElement)));
        for t = 1:length(temp)
            if abs(temp(t)) == 1
                Lr(j,nonZerosElement(t)) = 2*19.07*temp(t);
            else
                Lr(j,nonZerosElement(t)) = 2*atanh(temp(t));
            end
        end
    end
    %---------get the codeword and check if it is right(H*data=0)-----------
    for i = 1:n/2
        sumi1 = sum(Lr(:,i));
        sumi2 = sum(Lr(:,i+n/2));
        LQ(i) = p(i) + sum(Lr(:,i));
        LQ(i + n/2) = p(i + n/2) + sum(Lr(:,i+n/2));
        if (pdfT(i) > pdfA(i) && pdfT(i) > pdfG(i) && receivedSignal(2,i) < muG(2) &&  receivedSignal(1,i) < muA(1))
            Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
            Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
            Pz1L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp(-(sumi1)/2);
            Pz0L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp((sumi1)/2);
            sumi12 = log(((pdfA(i))/(pdfA(i)+pdfC(i))*Pz0L1+(pdfT(i)/(pdfG(i)+pdfT(i)))*Pz1L1)/(pdfC(i)/(pdfC(i)+pdfA(i))*Pz0L1+(pdfG(i)/(pdfG(i)+pdfT(i)))*Pz1L1));
            sumi21 = log(((pdfA(i))/(pdfA(i)+pdfT(i))*Pz0L2+(pdfC(i)/(pdfC(i)+pdfG(i)))*Pz1L2)/(pdfT(i)/(pdfA(i)+pdfT(i))*Pz0L2+(pdfG(i)/(pdfC(i)+pdfG(i)))*Pz1L2));
            LQ(i) = LQ(i)+sumi21;
            LQ(i+n/2) = LQ(i+n/2)+sumi12;
        end
        data(i) = double(LQ(i) < 0);
        data(i+n/2) = double(LQ(i+n/2) < 0);
    end
    re = H*data;
    ifend = 1;
    for i = 1:length(re)
        if(mod(re(i),2) == 1)
            ifend = 0;
        end
    end
    if(ifend)
        break;
    end
    %-------------variable to check message------------------------------------
    for i = 1:n/2
        sumi1 = sum(Lr(:,i));
        sumi2 = sum(Lr(:,i+n/2));
        hi1 = H(:,i);
        hi2 = H(:,i+n/2);
        nonZerosElementi1 = find(hi1~=0);
        nonZerosElementi2 = find(hi2~=0);
        Lq(nonZerosElementi1,i) = (p(i)+sumi1) - Lr(nonZerosElementi1,i)+1.0e-10;
        Lq(nonZerosElementi2,i+n/2) = p(i+n/2)+sumi2- Lr(nonZerosElementi2,i+n/2) +1.0e-10;
        if (pdfT(i) > pdfA(i) && pdfT(i) > pdfG(i) && receivedSignal(2,i) < muG(2) &&  receivedSignal(1,i) < muA(1))
            Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
            Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
            Pz1L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp(-(sumi1)/2);
            Pz0L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp((sumi1)/2);
            sumi12 = log(((pdfA(i))/(pdfA(i)+pdfC(i))*Pz0L1+(pdfT(i)/(pdfG(i)+pdfT(i)))*Pz1L1)/(pdfC(i)/(pdfC(i)+pdfA(i))*Pz0L1+(pdfG(i)/(pdfG(i)+pdfT(i)))*Pz1L1));
            sumi21 = log(((pdfA(i))/(pdfA(i)+pdfT(i))*Pz0L2+(pdfC(i)/(pdfC(i)+pdfG(i)))*Pz1L2)/(pdfT(i)/(pdfA(i)+pdfT(i))*Pz0L2+(pdfG(i)/(pdfC(i)+pdfG(i)))*Pz1L2));
            Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i) + sumi21;
            Lq(nonZerosElementi2,i+n/2) = Lq(nonZerosElementi2,i+n/2) + sumi12;
        end
    end
    %     toc
end
end


