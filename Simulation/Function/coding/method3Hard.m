%I use LLR and SPA to do the LDPC decode
%reference
%1.Channel code classical and Modern by William E.Ryan and Shu Lin
%Chapter5.4
%2. Maltab Document https://www.mathworks.com/help/comm/ref/ldpcdecoder.html
function [data1,data2] = method3Hard(llr1,llr2, H1,H2,receivedSignal,P)
Pac = P(1);
Pat = P(2);
Pag = P(3);
Pct = P(4);
Pcg = P(5);
Ptg = P(6);

Paa = 1 - Pag - Pac - Pat;
Pgg = 1 - Pag - Pcg - Ptg;

Lu1 = (Pac + Pcg);
Lu2 = (Pat + Ptg);
% Lu1 = 0.35;
% Lu2 = 0.55;
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

for i = 1:length(llr1)
    p1(i) = llr1(i);
    p2(i) = llr2(i);
    Lq1(H1(:,i)==1,i) = p1(i);%initial variabl to check
    Lq2(H2(:,i)==1,i) = p2(i);
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
        LQ1(i) = p1(i) + sum(Lr1(:,i));
        LQ2(i) = p2(i) + sum(Lr2(:,i));
        if (receivedSignal(i) == 1)
            sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu1)/(1-Lu1))));
            sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu1)/(1-Lu1))));
            LQ1(i) = LQ1(i)+sumi21;
            LQ2(i) = LQ2(i)+sumi12;
        elseif (receivedSignal(i) == 2)
            sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu2)/(1-Lu2))));
            sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu2)/(1-Lu2))));
            LQ1(i) = LQ1(i)+sumi21;
            LQ2(i) = LQ2(i)+sumi12;
%         elseif (receivedSignal(i) == 3)
%             sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu3)/(1-Lu3))));
%             sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu3)/(1-Lu3))));
%             LQ(i) = LQ(i)+sumi21;
%             LQ(i+n/2) = LQ(i+n/2)+sumi12;
%         elseif (receivedSignal(i) == 0)
%             sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu0)/(1-Lu0))));
%             sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu0)/(1-Lu0))));
%             LQ(i) = LQ(i)+sumi21;
%             LQ(i+n/2) = LQ(i+n/2)+sumi12;
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
    for i = 1:n/2
        nonZerosElementi1 = find(H1(:,i)~=0);
        nonZerosElementi2 = find(H2(:,i)~=0);
        sumi1 = sum(Lr1(:,i));
        sumi2 = sum(Lr2(:,i));
        Lq1(nonZerosElementi1,i) = p1(i)+sumi1- Lr1(nonZerosElementi1,i);
        Lq2(nonZerosElementi2,i) = p2(i)+sumi2- Lr2(nonZerosElementi2,i);
        if (receivedSignal(i) == 1)
            sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu1)/(1-Lu1))));
            sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu1)/(1-Lu1))));
            
            Lq1(nonZerosElementi1,i) = Lq1(nonZerosElementi1,i)+sumi21;
            Lq2(nonZerosElementi2,i) =Lq2(nonZerosElementi2,i)+sumi12;
        elseif (receivedSignal(i) == 2)
            sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu2)/(1-Lu2))));
            sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu2)/(1-Lu2))));
            
            Lq1(nonZerosElementi1,i) = Lq1(nonZerosElementi1,i)+sumi21;
            Lq2(nonZerosElementi2,i) =Lq2(nonZerosElementi2,i)+sumi12;
%         elseif (receivedSignal(i) == 3)
%             sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu3)/(1-Lu3))));
%             sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu3)/(1-Lu3))));
%             
%             Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i)+sumi21;
%             Lq(nonZerosElementi2,i+n/2) =Lq(nonZerosElementi2,i+n/2)+sumi12;
%         elseif (receivedSignal(i) == 0)
%             sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu0)/(1-Lu0))));
%             sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu0)/(1-Lu0))));
%             
%             Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i)+sumi21;
%             Lq(nonZerosElementi2,i+n/2) =Lq(nonZerosElementi2,i+n/2)+sumi12;
        end
    end
end
end




