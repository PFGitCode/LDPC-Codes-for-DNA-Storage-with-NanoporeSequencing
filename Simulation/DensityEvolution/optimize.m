clear
decodeMethod = 'method3'; % baseline, method1, method2, method3, method4
channelError = 1.27;
fprintf("%s, channelError: %d\n",decodeMethod, channelError);

errorModel = "data/testData" + num2str(channelError)+".mat";
load(errorModel,'mu','matrix');  % mu = [muA,muC,muT,muG]; matrix = [mA,mC,mT,mG];
acc = 0.01;
edge = 30/acc;
capSoft = capSoft(matrix,mu);

[chan1,chan2,P1,P2,checkMatrix] = pdfmSoft(matrix,mu,acc);
ext = [-30 acc 2*edge+1];
mapping = [-10 0.0002 50000];
% vard = [0	0.452989470443867	0.0100000000000000	0	0	0	0	0	0.0100000000000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.0100000000000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.517010529556133	0	0	0	0	0];
% vard = [0,0.457051136858101,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.126257734001549,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.413370978067786,0,0,0,0,0.00332015107256263];
% vard = [0,0.457063126714577,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.125520487724686,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.417416385560738,0,0,0,0,0];
% vard = [0	0.445142877710581	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.100000000000000	0.100000000000000	0	0	0	0.100000000000000	0	0	0	0	0	0.254857122289419	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
% vard = [0;0.445060896151953;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.269359271975545;0;0;0;0;0;0;0;0;0;0.285579831872501;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]';
% vard = [0,0.450780190974063,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.261283438993563,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.287936370032375,0];
vard = [0	0.462706514204586	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.0100000000000000	0	0	0.0100000000000000	0	0	0.0100000000000000	0.507293485795414];
chkd = [zeros(1,7),0.7188,0.2812];
dv = 50;
dc = 9;
iter = 100;
stop_pe = 1e-10;
sumvard = 0;
sumchkd = 0;
for i = 1: numel(vard)
    sumvard = sumvard+ vard(i)/i;
end
for i = 1: numel(chkd)
    sumchkd = sumchkd+ chkd(i)/i;
end
rate = 1 - sumchkd/sumvard;

[result_pe1, result_pe2, A1, A2] = de_irregular2D_opt(chan1,chan2,iter,ext,mapping,stop_pe,vard,chkd,dv,dc,P1,P2,decodeMethod);
bestV = vard;
bestResult = result_pe2(end);
bestA1 = A1;
bestA2 = A2;
for iterOp = 1:1000
    pl1 = zeros(iter,1);
    pl2 = zeros(iter,1);
    h1 = zeros(1,dv);
    h2 = zeros(1,dv);
    for j = 2:dv
        pl1 = pl1 + A1(:,j)*vard(j);
    end
    for j = 2:dv
        pl2 = pl2 + A2(:,j)*vard(j);
    end
    lambda = zeros(size(vard));
    aeq = ones(size(vard));
    beq = 1;
    
    ra = 0;
    aeq2 = zeros(size(vard));
    for i = 1:dv
        ra = ra + vard(i)/i;
        aeq2(i) = 1/i;
    end
    
    beq = [beq,ra];
    aeq = [aeq;aeq2];
    
    pl1Delta = zeros(iter,1);
    for i = 2:iter
        pl1Delta(i) = pl1(i-1) - pl1(i);
    end
    
    sigma = 0.001;
    AA1 = -A1(2:end,:);
    b1 = pl1Delta(2:end,:)*sigma - pl1(2:end);
    %         b1(1) = 0;
    AA1_2 = A1(2:end,:);
    b1_2 = pl1Delta(2:end,:)*sigma + pl1(2:end);
    AA1 = [AA1;AA1_2];
    b1 = [b1;b1_2];
    %
    AA1 = [AA1;AA1_2];
    b1 = [b1;pl1(1:end-1)];
    f1 = zeros(size(vard));
%             c = 0;
    for i = 2:length(pl1Delta(pl1Delta~=0))+1
        f1 =  A1(i,:)/pl1Delta(i);
%                     c = c + pl1(i)/pl1Delta(i);
    end
%         f1 = [f1, -c];
    lb = zeros(dv,1);
    up = ones(dv,1);
    up(1) = 0;
    x = linprog(f1,AA1,b1,aeq,beq,lb,up);
    if isempty(x)
        continue;
    end
    [result_pe1, result_pe2, A1,A2] = de_irregular2D_opt(chan1,chan2,iter,ext,mapping,stop_pe,x',chkd,dv,dc,P1,P2,decodeMethod);
    if result_pe2 < bestResult
        bestV = x'
        bestResult = result_pe2;
        bestA1 = A1;
        bestA2 = A2;
    else
        1;
    end
%     vard = x';
end