%I use LLR and SPA to do the LDPC decode
%reference
%1.Channel code classical and Modern by William E.Ryan and Shu Lin
%Chapter5.4
%2. Maltab Document https://www.mathworks.com/help/comm/ref/ldpcdecoder.html
function data = method1DecoderHard(llr, H,receivedSignal)
[k, n] = size(H);
data = zeros(n,1);
p = zeros(n,1);
Lq = zeros(k,n);
Lr = zeros(k,n);
LQ = zeros(n,1);

alpha = 0.3;
for i = 1:length(llr)
    p(i) = llr(i);%llr of p0 and p1
    Lq(H(:,i)==1,i) = p(i);%initial variabl to check
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
        LQ(i) = p(i) + sum(Lr(:,i));
        LQ(i + n/2) = p(i + n/2) + sum(Lr(:,i + n/2));
        if (receivedSignal(i) == 1||receivedSignal(i) == 2)
            LQ(i) = LQ(i) - alpha*(sumi2);
            LQ(i + n/2) =  LQ(i + n/2) - alpha*(sumi1);
        end
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
        Lq(nonZerosElementi1,i) = p(i) +sumi1 - Lr(nonZerosElementi1,i)+1.0e-10;
        Lq(nonZerosElementi2,i+n/2) = p(i+n/2) +sumi2 - Lr(nonZerosElementi2,i+n/2)+1.0e-10;
        if (receivedSignal(i) == 1||receivedSignal(i) == 2)
            Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i) - alpha*(sumi2);
            Lq(nonZerosElementi2,i+n/2) = Lq(nonZerosElementi2,i+n/2) - alpha*(sumi1);
        end
    end
end
end

