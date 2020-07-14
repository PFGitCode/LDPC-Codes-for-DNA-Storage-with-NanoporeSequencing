%I use LLR and SPA to do the LDPC decode
%reference
%1.Channel code classical and Modern by William E.Ryan and Shu Lin
%Chapter5.4
%2. Maltab Document https://www.mathworks.com/help/comm/ref/ldpcdecoder.html
function data = method3Hard(llr1,llr2, H1,H2,receivedSignal,P)
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
[k, n] = size(H);
data = zeros(n,1);
p = zeros(n,1);
Lq = zeros(k,n);
Lr = zeros(k,n);
LQ = zeros(n,1);

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
        if (receivedSignal(i) == 1)
            sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu1)/(1-Lu1))));
            sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu1)/(1-Lu1))));
            LQ(i) = LQ(i)+sumi21;
            LQ(i+n/2) = LQ(i+n/2)+sumi12;
        elseif (receivedSignal(i) == 2)
            sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu2)/(1-Lu2))));
            sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu2)/(1-Lu2))));
            LQ(i) = LQ(i)+sumi21;
            LQ(i+n/2) = LQ(i+n/2)+sumi12;
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
        Lq(nonZerosElementi1,i) = p(i)+sumi1- Lr(nonZerosElementi1,i);
        Lq(nonZerosElementi2,i+n/2) = p(i+n/2)+sumi2- Lr(nonZerosElementi2,i+n/2);
        if (receivedSignal(i) == 1)
            sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu1)/(1-Lu1))));
            sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu1)/(1-Lu1))));
            
            Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i)+sumi21;
            Lq(nonZerosElementi2,i+n/2) =Lq(nonZerosElementi2,i+n/2)+sumi12;
        elseif (receivedSignal(i) == 2)
            sumi21 = 2*atanh(tanh(0.5*(sumi2))*tanh(0.5*log((Lu2)/(1-Lu2))));
            sumi12 = 2*atanh(tanh(0.5*(sumi1))*tanh(0.5*log((Lu2)/(1-Lu2))));
            
            Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i)+sumi21;
            Lq(nonZerosElementi2,i+n/2) =Lq(nonZerosElementi2,i+n/2)+sumi12;
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




