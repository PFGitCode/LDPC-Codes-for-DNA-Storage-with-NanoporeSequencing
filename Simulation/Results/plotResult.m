factor = [0.60, 0.65,0.7, 0.75];
error0 = [0.130383, 0.175307, 0.208228, 0.230341];
error1 = [0.053081, 0.144877, 0.196834, 0.225259]; 
error2 = [0,0.005165,0,0];
error3 = [0,0.004306,0,0];
error4 = [5.0200000e-05,0.0027268, 0.0553712,0];
errorq = [6.9995333e-06,0.0016133, 0.1445572,0];

time0 = [1.259748, 1.220378, 1.243236, 1.112529];
time1 = [0.857447, 1.213813, 1.298658, 1.314779];
time2 = [0.219883,0.418113,0,0];
time3 = [0.239829,0.456138,0,0];
time4 = [0.907315,1.563351,2.845757,4.664305];
timeq = [3.975534,6.509376,0,0,];

base = plot(factor,error0);
hold on;
method4p = plot(factor,error4,'--');
hold on;
method3p = plot(factor,error3,'-.');
hold on;
method2p = plot(factor,error2,'-.');
hold on;
method1p = plot(factor,error1,'-.');
hold on;
quaterp = plot(factor,errorq,'.');
% legend([base2 method2 quaternary2 quaternary4], ...
% ["Baseline,p_1=2{\lambda}","Baseline,p_1=4{\lambda}","Algorithm 3,p_1=2{\lambda}","Algorithm 3,p_1=4{\lambda}","quaternary,p_1=2{\lambda}","quaternary,p_1=4{\lambda}"]);
legend([base,method1p,method2p, method3p,method4p,quaterp], ...
["base","method1","method2","method3","method4","Quaternary"]);
xlabel('Variance Factor');
ylabel('Bit Error Rate');
set(gca, 'YScale', 'log');
