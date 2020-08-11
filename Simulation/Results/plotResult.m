factor = [0.5,0.55,0.60, 0.65,0.7];
% time0 = [0.992961,1.138987,1.198494, 1.371803, 1.380006];
% time1 = [0.857447, 1.213813, 1.298658, 1.314779];
time2 = [0.175687,0.216331,0.2779,0.429261094,0.711139130056237];
time3 = [0.179212,0.235267,0.329232337196443,0.480870201300001,0.758963252746631];
time4 = [0.322751,0.435206,0.624231,0.927640,2.108794];
timeq = [2.512537,2.856301,3.650136,6.706798,16.100932];

% base = plot(factor,time0);
% hold on;
method4p = plot(factor,time4,'-+');
hold on;
method3p = plot(factor,time3,'-^');
hold on;
method2p = plot(factor,time2,'-*');
hold on;
% method1p = plot(factor,time1,':');
% hold on;
quaterp = plot(factor,timeq,'-o');
legend([method2p, method3p,method4p,quaterp], ...
["method2","method3","method4","Quaternary"]);
xlabel('Variance Factor');
ylabel('Bit Error Rate');
set(gca, 'YScale', 'log');

% hard
% factor = [0.5,0.55,0.57,0.58,0.60, 0.65,0.7, 0.75];
% error0 = [0.036304,0.085444,0.108599,0.110732,0.130383, 0.175307, 0.208228, 0.230341];
% error1 = [0.000134,0.004357,0.014965,0.023557,0.053081, 0.144877, 0.196834, 0.225259]; 
% error2 = [0,0,1.06e-06,4.44e-06,5.472e-05,0.0072048,0.0945926,0.222921];
% error3 = [0,0,5.32e-06,1.57e-05,1.016e-04,0.0061018,0.0701196,0.178733];
% error4 = [0,0,0,6.00e-07,1.450e-05,0.0027268, 0.0553712,0.171086];
% errorq = [0,0,0,0,6.9995333e-06,0.0016133, 0.046436,0.1445572];
% 
% base = plot(factor,error0,'-s');
% hold on;
% method4p = plot(factor,error4,'-*');
% hold on;
% method3p = plot(factor,error3,'-+');
% hold on;
% method2p = plot(factor,error2,'-x');
% hold on;
% method1p = plot(factor,error1,'-^');
% hold on;
% quaterp = plot(factor,errorq,'-o');
% % legend([base2 method2 quaternary2 quaternary4], ...
% % ["Baseline,p_1=2{\lambda}","Baseline,p_1=4{\lambda}","Algorithm 3,p_1=2{\lambda}","Algorithm 3,p_1=4{\lambda}","quaternary,p_1=2{\lambda}","quaternary,p_1=4{\lambda}"]);
% legend([base,method1p,method2p, method3p,method4p,quaterp], ...
% ["base","method1","method2","method3","method4","Quaternary"]);
% xlabel('Variance Factor');
% ylabel('Bit Error Rate');
% set(gca, 'YScale', 'log');

% soft
% factor = [0.70,0.75,0.78,0.80,0.85,0.90,0.95];
% error0 = [0.0424,0.095011,0.123297,0.138911,0.180662,0.212731,0.234923];
% error1 = [5.886e-05,0.002939,0.015905,0.034965,0.117265,0.190178,0.225995];
% error2 = [0,7.260000000000000e-07,6.840000000000000e-06,3.48e-05,6.555e-04,0.0079401,0.041338533333333];
% error3 = [0,0,1.22e-06,7.56e-06,1.602e-04,0.0022829,0.021];
% error4 = [0,0,0,2.14e-06,6.36e-05,0.001419766666667,0.014226200000000];
% 
% 
% base = plot(factor,error0,'-s');
% hold on;
% method4p = plot(factor,error4,'-*');
% hold on;
% method3p = plot(factor,error3,'-+');
% hold on;
% method2p = plot(factor,error2,'-x');
% hold on;
% method1p = plot(factor,error1,'-^');
% legend([base,method1p,method2p, method3p,method4p], ...
% ["base","method1","method2","method3","method4"]);
% xlabel('Variance Factor');
% ylabel('Bit Error Rate');
% set(gca, 'YScale', 'log');

% singleTime = [0.021884,0.022672,0.022061,0.058358,0.431393];
% X = categorical({'method1','method2','method3','method4','quaternary'});
% X = reordercats(X,{'method1','method2','method3','method4','quaternary'});
% bar(X,singleTime)