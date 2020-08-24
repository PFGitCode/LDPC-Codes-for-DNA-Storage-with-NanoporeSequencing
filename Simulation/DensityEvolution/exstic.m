channelError = 1.27;
errorModel = "data/testData" + num2str(channelError)+".mat";
decodeMethod = 'method3';
load(errorModel);
acc = 0.01;
edge = 30/acc;
[chan1,chan2,P1,P2,checkMatrix] = pdfmSoft(matrix,mu,acc);
cap = capSoft(matrix,mu);
ext = [-30 acc 2*edge+1];
mapping = [-10 0.0002 50000];
%
chkd = [zeros(1,7),0.7188,0.2812];
vard = [0,0.2340,0.3051,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1570,0.3038];
vardLen = 50;
elementNum = 5;
choose = combnk(1:vardLen,elementNum);
choose_rand = choose(randperm(length(choose)),:);
bestResult = 1;
% bestV = [0,0.476785119531555,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.523214880468445];
ra = 0;
for i = 1:length(vard)
    ra = ra + vard(i)/i;
end
dv = vardLen;
dc = 9;
for ii = 1:nchoosek(vardLen,elementNum)
    
    chose = choose_rand(ii,:);
    %     chose = find(vard~=0);
    if chose(1) == 1
        continue;
    end
    aeq = ones(1,elementNum);
    aeq2 = ones(1,elementNum)./chose;
    beq = 1;
    beq = [beq,ra];
    aeq = [aeq;aeq2];
    b1 = 0;
    f = zeros(2,dv);
    f(1,chose(1)) = -1;
    f(2,chose(2)) = -1;
    
    b = -ones(1,2)*0.1;
    %         b(2) = -0.2;
    x = linprog([],[],[],aeq,beq,ones(1,elementNum)*0.01,ones(1,elementNum));
    if min(x) < 0
        continue;
    end
    if isempty(x)
        1;
        continue;
    end
    vard = zeros(1,vardLen);
    vard(chose) = x';
    dv = find(vard~=0, 1, 'last' );
    iter = 10;
    stop_pe = 1e-6;
    sumvard = 0;
    sumchkd = 0;
    for i = 1: numel(vard)
        sumvard = sumvard+ vard(i)/i;
    end
    for i = 1: numel(chkd)
        sumchkd = sumchkd+ chkd(i)/i;
    end
    rate = 1 - sumchkd/sumvard;
    [result_pe1, result_pe2, A1,A2] = de_irregular2D_opt(chan1,chan2,iter,ext,mapping,stop_pe,vard,chkd,dv,dc,P1,P2,decodeMethod);
    if result_pe2 < bestResult
        bestV = vard;
        bestResult = result_pe2(end);
        bestA1 = A1;
        bestA2 = A2;
        1;
    end
end