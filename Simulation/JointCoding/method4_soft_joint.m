%I use LLR and SPA to do the LDPC decode
%reference
%1.Channel code classical and Modern by William E.Ryan and Shu Lin
%Chapter5.4
%2. Maltab Document https://www.mathworks.com/help/comm/ref/ldpcdecoder.html
function data = method4_soft_joint(llr, H,receivedSignal,mu,matrix)
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

for i = 1:length(llr)
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
    %---------------------message from checl to variable----------
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
        LQ(i) = sumi1;
        LQ(i + n/2) = sumi2;
        sumi21 = maxStar(log(pdfA(i))+sumi2, log(pdfC(i))) - maxStar(log(pdfT(i))+sumi2, log(pdfG(i)));
        sumi12 = maxStar(log(pdfA(i))+sumi1, log(pdfT(i))) - maxStar(log(pdfC(i))+sumi1, log(pdfG(i)));
        LQ(i) = LQ(i)+sumi21;
        LQ(i + n/2) = LQ(i + n/2)+sumi12;
        data(i) = LQ(i) < 0;
        data(i+n/2) = LQ(i+n/2) < 0;
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
        nonZerosElementi1 = find(H(:,i)~=0);
        nonZerosElementi2 = find(H(:,i+n/2)~=0);
        sumi1 = sum(Lr(:,i));
        sumi2 = sum(Lr(:,i+n/2));
        Lq(nonZerosElementi1,i) = sumi1 - Lr(nonZerosElementi1,i);
        Lq(nonZerosElementi2,i+n/2) = sumi2 - Lr(nonZerosElementi2,i+n/2);
        sumi21 = maxStar(log(pdfA(i))+sumi2, log(pdfC(i))) - maxStar(log(pdfT(i))+sumi2, log(pdfG(i)));
        sumi12 = maxStar(log(pdfA(i))+sumi1, log(pdfT(i))) - maxStar(log(pdfC(i))+sumi1, log(pdfG(i)));
        Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i)+sumi21;
        Lq(nonZerosElementi2,i+n/2) = Lq(nonZerosElementi2,i+n/2)+sumi12;
    end
end
end

function s = maxStar(x1,x2)
x = [x1,x2];
y = max(x1,x2);
x = bsxfun(@minus,x,y);
s = y + log(sum(exp(x)));
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end
end



