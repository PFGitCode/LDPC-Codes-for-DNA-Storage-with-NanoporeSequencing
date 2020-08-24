clear
decodeMethod = 'baseline'; % baseline, method1, method2, method3, method4
channelError = 1.2;
fprintf("%s, channelError: %d\n",decodeMethod, channelError);

errorModel = "data/testData" + num2str(channelError)+".mat";
load(errorModel,'mu','matrix');  % mu = [muA,muC,muT,muG]; matrix = [mA,mC,mT,mG];
acc = 0.01;
edge = 30/acc;
capSoft = capSoft(matrix,mu);

[chan1,chan2,P1,P2,checkMatrix] = pdfmSoft(matrix,mu,acc);
ext = [-30 acc 2*edge+1];
mapping = [-10 0.0002 50000];
% vard = [0,0,1];
% chkd = [0,0,0,0,0,1];
% dv = 3;
% dc = 6;

vard = [0,0.4734,zeros(1,36),5.6316e-05,zeros(1,22),0.0020,zeros(1,37),0.5239];
% vard = [0	0.452989470443867	0.0100000000000000	0	0	0	0	0	0.0100000000000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.0100000000000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.517010529556133	0	0	0	0	0];
% vard = [0	0.462706514204586	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.0100000000000000	0	0	0.0100000000000000	0	0	0.0100000000000000	0.507293485795414];
% vard = [0	0.462703768710953	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.0170000234259019	0	0	0	0	0	0	0.520296207863145];
chkd = [zeros(1,7),0.7188,0.2812];
dv = 100;
dc = 9;
iter = 10000;
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

result_pe1 = de_irregular2D(chan1,chan2,iter,ext,mapping,stop_pe,vard,chkd,dv,dc, P1,P2,decodeMethod);
