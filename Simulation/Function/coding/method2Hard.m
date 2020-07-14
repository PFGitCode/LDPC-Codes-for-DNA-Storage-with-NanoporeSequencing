%I use LLR and SPA to do the LDPC decode
%reference
%1.Channel code classical and Modern by William E.Ryan and Shu Lin
%Chapter5.4
%2. Maltab Document https://www.mathworks.com/help/comm/ref/ldpcdecoder.html
function data = method2Hard(llr,H,receivedSignal,Pac,Pat,Pag,Pct,Pcg,Ptg)
[k, n] = size(H);
data = zeros(n,1);
p = zeros(n,1);
Lq = zeros(k,n);
Lr = zeros(k,n);
LQ = zeros(n,1);
% P1 = Lu;
% P2 = Pct;
Pcc = 1-Pac - Pcg-Pct;
Ptt = 1-Pat - Pct-Ptg;
for i = 1:length(llr)
    p(i) = llr(i);%llr of p0 and p1
    Lq(H(:,i)==1,i) = p(i);%initial variabl to check
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
        if (receivedSignal(i) == 1)
            Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
            Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
            Pz1L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp(-(sumi1)/2);
            Pz0L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp((sumi1)/2);
            sumi21 = log(((Pac)/(Pac+Pcc)*Pz0L2+(Pct/(Pcg+Pct))*Pz1L2)/((Pcc)/((Pcc+Pac))*Pz0L2+(Pcg/(Pcg+Pct))*Pz1L2));
            sumi12 = log(((Pac)/((Pac+Pct))*Pz0L1+((Pcc)/(Pcc+Pcg))*Pz1L1)/((Pct)/((Pac+Pct))*Pz0L1+(Pcg/(Pcc+Pcg))*Pz1L1));
            %             sumi21 = log(((P1/(P1+P2))*Pz0L2+(1-2*P1-P2)/((1-P1-P2))*Pz1L2)/(((P2/(P1+P2))/((1-P1-P2)))*Pz0L2+(P1/(1-P1-P2))*Pz1L2));
            %             sumi12 = log(((P1)/((1-P1-P2))*Pz0L1+(P2/(P1+P2))*Pz1L1)/((1-2*P1-P2)/((1-P1-P2))*Pz0L1+(P1/(P1+P2))*Pz1L1));
            LQ(i) = LQ(i)+sumi21;
            LQ(i+n/2) = LQ(i+n/2)+sumi12;
        elseif  receivedSignal(i) == 2
            Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
            Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
            Pz1L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp(-(sumi1)/2);
            Pz0L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp((sumi1)/2);
            sumi21 = log(((Pat)/(Pat+Pct)*Pz0L2+(Ptt/(Ptg+Ptt))*Pz1L2)/((Pct)/((Pct+Pat))*Pz0L2+(Ptg/(Ptg+Ptt))*Pz1L2));
            sumi12 = log(((Pat)/((Pat+Ptt))*Pz0L1+((Pct)/(Pct+Ptg))*Pz1L1)/((Ptt)/((Pat+Ptt))*Pz0L1+(Ptg/(Pct+Ptg))*Pz1L1));
            LQ(i) = LQ(i)+sumi21;
            LQ(i+n/2) = LQ(i+n/2)+sumi12;
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
        sumi1 = sum(Lr(:,i));
        sumi2 = sum(Lr(:,i+n/2));
        hi1 = H(:,i);
        hi2 = H(:,i+n/2);
        nonZerosElementi1 = find(hi1~=0);
        nonZerosElementi2 = find(hi2~=0);
        Lq(nonZerosElementi1,i) = (p(i)+sumi1) - Lr(nonZerosElementi1,i)+1.0e-10;
        Lq(nonZerosElementi2,i+n/2) = p(i+n/2)+sumi2- Lr(nonZerosElementi2,i+n/2) +1.0e-10;
        if (receivedSignal(i) == 1)
            Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
            Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
            Pz1L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp(-(sumi1)/2);
            Pz0L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp((sumi1)/2);
            sumi21 = log(((Pac)/(Pac+Pcc)*Pz0L2+(Pct/(Pcg+Pct))*Pz1L2)/((Pcc)/((Pcc+Pac))*Pz0L2+(Pcg/(Pcg+Pct))*Pz1L2));
            sumi12 = log(((Pac)/((Pac+Pct))*Pz0L1+((Pcc)/(Pcc+Pcg))*Pz1L1)/((Pct)/((Pac+Pct))*Pz0L1+(Pcg/(Pcc+Pcg))*Pz1L1));
            %             sumi21 = log(((P1/(P1+P2))*Pz0L2+(1-2*P1-P2)/((1-P1-P2))*Pz1L2)/(((P2/(P1+P2))/((1-P1-P2)))*Pz0L2+(P1/(1-P1-P2))*Pz1L2));
            %             sumi12 = log(((P1)/((1-P1-P2))*Pz0L1+(P2/(P1+P2))*Pz1L1)/((1-2*P1-P2)/((1-P1-P2))*Pz0L1+(P1/(P1+P2))*Pz1L1));
            Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i) + sumi21;
            Lq(nonZerosElementi2,i+n/2) = Lq(nonZerosElementi2,i+n/2) + sumi12;
        elseif receivedSignal(i) == 2
            Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
            Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
            Pz1L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp(-(sumi1)/2);
            Pz0L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp((sumi1)/2);
            sumi21 = log(((Pat)/(Pat+Pct)*Pz0L2+(Ptt/(Ptg+Ptt))*Pz1L2)/((Pct)/((Pct+Pat))*Pz0L2+(Ptg/(Ptg+Ptt))*Pz1L2));
            sumi12 = log(((Pat)/((Pat+Ptt))*Pz0L1+((Pct)/(Pct+Ptg))*Pz1L1)/((Ptt)/((Pat+Ptt))*Pz0L1+(Ptg/(Pct+Ptg))*Pz1L1));
            Lq(nonZerosElementi1,i) = Lq(nonZerosElementi1,i) + sumi21;
            Lq(nonZerosElementi2,i+n/2) = Lq(nonZerosElementi2,i+n/2) + sumi12;
        end
    end
    %     toc
end
end


