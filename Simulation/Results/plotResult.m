% hard
% factor = [0.5,0.55,0.57,0.58,0.60, 0.65,0.7, 0.75];
% error0 = [0.036304,0.085444,0.108599,0.110732,0.130383, 0.175307, 0.208228, 0.230341];
% error1 = [0.000134,0.004357,0.014965,0.023557,0.053081, 0.144877, 0.196834, 0.225259]; 
% error2 = [0,0,1.06e-06,4.44e-06,5.472e-05,0.0072048,0.0945926,0.222921];
% error3 = [0,0,5.32e-06,1.57e-05,1.016e-04,0.0061018,0.0701196,0.178733];
% error4 = [0,0,0,6.00e-07,1.450e-05,0.0027268, 0.0553712,0.171086];
% errorq = [0,0,0,0,6.9995333e-06,0.0016133, 0.046436,0.1445572];
% 
% time0 = [1.259748, 1.220378, 1.243236, 1.112529];
% time1 = [0.857447, 1.213813, 1.298658, 1.314779];
% time2 = [0.219883,0.418113,0,0];
% time3 = [0.239829,0.456138,0,0];
% time4 = [0.907315,1.563351,2.845757,4.664305];
% timeq = [3.975534,6.509376,0,0,];

% soft
factor = [0.70,0.75,0.78,0.80,0.85,0.90,0.95];
error0 = [0.0424,0.095011,0,0.138911,0.180662,0.212731,0.234923];
error1 = [5.886e-05,0.002939,0,0.034965,0.117265,0.190178,0.225995];
error2 = [0,2.26e-06,0,3.48e-05,6.555e-04,0.0079401,0.041338533333333];
error3 = [0,0,2.22e-06,7.56e-06,1.602e-04,0.0022829,0.021];
error4 = [0,0,0,2.14e-06,6.36e-05,0.001419766666667,0.014226200000000];
errorq = [0,0,0,0,0,0,0];



base = plot(factor,error0);
hold on;
method4p = plot(factor,error4,'--');
hold on;
method3p = plot(factor,error3,'-.');
hold on;
method2p = plot(factor,error2,'-.');
hold on;
method1p = plot(factor,error1,':');
hold on;
quaterp = plot(factor,errorq,'-.');
% legend([base2 method2 quaternary2 quaternary4], ...
% ["Baseline,p_1=2{\lambda}","Baseline,p_1=4{\lambda}","Algorithm 3,p_1=2{\lambda}","Algorithm 3,p_1=4{\lambda}","quaternary,p_1=2{\lambda}","quaternary,p_1=4{\lambda}"]);
legend([base,method1p,method2p, method3p,method4p,quaterp], ...
["base","method1","method2","method3","method4","Quaternary"]);
xlabel('Variance Factor');
ylabel('Bit Error Rate');
set(gca, 'YScale', 'log');
