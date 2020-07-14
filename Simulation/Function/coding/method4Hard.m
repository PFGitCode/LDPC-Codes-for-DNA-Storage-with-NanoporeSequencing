%I use LLR and SPA to do the LDPC decode
%reference
%1.Channel code classical and Modern by William E.Ryan and Shu Lin
%Chapter5.4
%2. Maltab Document https://www.mathworks.com/help/comm/ref/ldpcdecoder.html
% P = [Pac,Pat,Pag,Pct,Pcg,Ptg]
function  [data1,data2]  = method4Hard(llr1,llr2, H1,H2,receivedSignal,P) 
Pac = P(1);
Pat = P(2);
Pag = P(3);
Pct = P(4);
Pcg = P(5);
Ptg = P(6);

P00a = 1-Pac-Pat-Pag;
P01a = Pac;
P10a = Pat;
P11a = Pag;

P00g = Pag;
P01g = Pcg;
P10g = Ptg;
P11g = 1-Pag-Pcg-Ptg;

P00c = Pac;
P01c = 1-Pac-Pcg-Pct;
P10c = Pct;
P11c = Pcg;

P00t = Pat;
P01t = Pct;
P10t = 1-Pat-Ptg-Pct;
P11t = Ptg;

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
    p1(i) = llr1(i);%llr of p0 and p1
    p2(i) = llr2(i);%llr of p0 and p1
    Lq1(H1(:,i)==1,i) = p1(i);%initial variabl to check
    Lq2(H2(:,i)==1,i) = p2(i);%initial variabl to check
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
        LQ1(i) = sumi1 ;
        LQ2(i) =  sumi2;
        if(receivedSignal(i) == 1)
            sumi21 = maxStar(log(P00c)+sumi2, log(P01c)) - maxStar(log(P10c)+sumi2, log(P11c));
            sumi12 = maxStar(log(P00c)+sumi1, log(P10c)) - maxStar(log(P01c)+sumi1, log(P11c));
        elseif(receivedSignal(i) == 2)
            sumi21 = maxStar(log(P00t)+sumi2, log(P01t)) - maxStar(log(P10t)+sumi2, log(P11t));
            sumi12 = maxStar(log(P00t)+sumi1, log(P10t)) - maxStar(log(P01t)+sumi1, log(P11t));
        elseif(receivedSignal(i) == 3)
            sumi21 = maxStar(log(P00g)+sumi2, log(P01g)) - maxStar(log(P10g)+sumi2, log(P11g));
            sumi12 = maxStar(log(P00g)+sumi1, log(P10g)) - maxStar(log(P01g)+sumi1, log(P11g));
        elseif(receivedSignal(i) == 0)
            sumi21 = maxStar(log(P00a)+sumi2, log(P01a)) - maxStar(log(P10a)+sumi2, log(P11a));
            sumi12 = maxStar(log(P00a)+sumi1, log(P10a)) - maxStar(log(P01a)+sumi1, log(P11a));
        end
        LQ1(i) = LQ1(i)+sumi21;
        LQ2(i) = LQ2(i)+sumi12;
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
        Lq1(nonZerosElementi1,i) = sumi1 - Lr1(nonZerosElementi1,i)+1.0e-10;
        Lq2(nonZerosElementi2,i) = sumi2 - Lr2(nonZerosElementi2,i)+1.0e-10;
        if(receivedSignal(i) == 1)
            sumi21 = maxStar(log(P00c)+sumi2, log(P01c)) - maxStar(log(P10c)+sumi2, log(P11c));
            sumi12 = maxStar(log(P00c)+sumi1, log(P10c)) - maxStar(log(P01c)+sumi1, log(P11c));
        elseif(receivedSignal(i) == 2)
            sumi21 = maxStar(log(P00t)+sumi2, log(P01t)) - maxStar(log(P10t)+sumi2, log(P11t));
            sumi12 = maxStar(log(P00t)+sumi1, log(P10t)) - maxStar(log(P01t)+sumi1, log(P11t));
        elseif(receivedSignal(i) == 3)
            sumi21 = maxStar(log(P00g)+sumi2, log(P01g)) - maxStar(log(P10g)+sumi2, log(P11g));
            sumi12 = maxStar(log(P00g)+sumi1, log(P10g)) - maxStar(log(P01g)+sumi1, log(P11g));
        elseif(receivedSignal(i) == 0)
            sumi21 = maxStar(log(P00a)+sumi2, log(P01a)) - maxStar(log(P10a)+sumi2, log(P11a));
            sumi12 = maxStar(log(P00a)+sumi1, log(P10a)) - maxStar(log(P01a)+sumi1, log(P11a));
        end
        Lq1(nonZerosElementi1,i) = Lq1(nonZerosElementi1,i)+sumi21;
        Lq2(nonZerosElementi2,i) = Lq2(nonZerosElementi2,i)+sumi12;
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


