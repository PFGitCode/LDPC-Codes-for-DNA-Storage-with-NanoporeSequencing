% r = 0.8;
%
% x1 = -0.1:0.01:1.6;
% x2 = -0.1:0.01:1.6;
% [X1, X2] = meshgrid(x1,x2);
% X = [X1(:) X2(:)];
%
% muA = [0.5,0.62];
% muT = [0.09,0.45];
% muC = [0.06,0.35];
% muG = [0.15,0.85];
%
% thetaA = [0.1,0.1]*r;
% thetaT = [0.2,0.2]*r;
% thetaC = [0.15,0.07]*r;
% thetaG = [0.1,0.3]*r;
%
% angleA = -pi/6;
% angleT = 0;
% angleC = 0;
% angleG = -pi/3;
%
% mA = generateMatrix(angleA,thetaA);
% mT = generateMatrix(angleT,thetaT);
% mC = generateMatrix(angleC,thetaC);
% mG = generateMatrix(angleG,thetaG);
% i = 1;
% for x1 = -0.1:0.01:1.6
%     for x2 = -0.1:0.01:1.6
%         yA = mvnpdf([x1,x2],muA,mA);
%         yC = mvnpdf([x1,x2],muC,mC);
%         yT = mvnpdf([x1,x2],muT,mT);
%         yG = mvnpdf([x1,x2],muG,mG);
%         [val, pos] = max([yA, yC, yT, yG]);
% %         scatter(x1, x2, 'ro', 'MarkerSize', 3);
%         group(i) = pos;
%         i = i+1;
% %         hold on
% %         if pos == 1
% %             y
% %         elseif pos == 2
% %             scatter(x1,x2, 'r');
% %         elseif pos == 3
% %             scatter(x1,x2, 'y');
% %         else
% %             scatter(x1,x2, 'g');
% %         end
%     end
% end
% % x = [-0.1:0.1:1.6];
% % y = [-0.1:0.1:1.6];
% gscatter(X(:,1),X(:,2),group);
%
% % X = [mvnrnd(muA,mA,2000);mvnrnd(muT,mT,2000);mvnrnd(muC,mC,2000);mvnrnd(muG,mG,2000)];
% XA = [mvnrnd(muA,mA,2000)];
% XC = mvnrnd(muC,mC,2000);
% % X = [mvnrnd(muT,mT,2000)];
% % GMModel = fitgmdist(X,4);
%
%
% figure
% % y = [zeros(2000,1);ones(2000,1);ones(2000,1)+1;ones(2000,1)+2];
%  yA = [zeros(2000,1)];
%  yC = [zeros(2000,1)];
% h = gscatter(XA(:,1),XA(:,2),yA);
% hold on
% gscatter(XC(:,1),XC(:,2),yC,'b');
%
% gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
% group = gca;
% % fcontour(gmPDF,[g.XLim g.YLim])
% % title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
% legend(h,'A','T','C','G')
% xlabel('Dwell Time');
% ylabel('Current Drop');
% hold off
%
%
% function m = generateMatrix(sigma,theta)
% theta_X = theta(1);
% theta_Y = theta(2);
% %
% s = [theta_X^2, 0; 0,theta_Y^2];
% B = [cos(sigma), sin(sigma);-sin(sigma), cos(sigma)];
% m = B'*s*B;
% end

r = 0.8;

x1 = -0.1:0.01:1.6;
x2 = -0.1:0.01:1.6;
[X1, X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

muA = [0.5,0.62];
muT = [0.09,0.45];
muC = [0.06,0.35];
muG = [0.15,0.85];

thetaA = [0.1,0.1]*r;
thetaT = [0.2,0.2]*r;
thetaC = [0.15,0.07]*r;
thetaG = [0.1,0.3]*r;

angleA = -pi/6;
angleT = 0;
angleC = 0;
angleG = -pi/3;

mA = generateMatrix(angleA,thetaA);
mT = generateMatrix(angleT,thetaT);
mC = generateMatrix(angleC,thetaC);
mG = generateMatrix(angleG,thetaG);

X = [mvnrnd(muA,mA,20000);mvnrnd(muT,mT,20000);mvnrnd(muC,mC,20000);mvnrnd(muG,mG,20000)];
GMModel = fitgmdist(X,4);


figure
% y = [zeros(2000,1);ones(2000,1);ones(2000,1)+1;ones(2000,1)+2];
y = zeros(80000,1);
for i = 1:80000
    x1 = X(i,1);
    x2 = X(i,2);
    yA = mvnpdf([x1,x2],muA,mA);
    yC = mvnpdf([x1,x2],muC,mC);
    yT = mvnpdf([x1,x2],muT,mT);
    yG = mvnpdf([x1,x2],muG,mG);
    [val, pos] = max([yA, yC, yT, yG]);
%     [vals, poss] = sort([yA, yC, yT, yG]);
%     if vals(4) - vals(3) < 1e-1
%         y(i) = 0;
%     else
%         y(i) = pos;
%     end
%     if (abs(yA-yT) < 0.15 || abs(yG-yT) < 0.15) && (x1 < 0.45 && x1 > -0.3 && x2 < 0.7 && x2 > 0.4 && x1+x2 <1)
%         y(i) = 0;
%     else
        y(i) = pos;
%     end
end
y(y == 3) = y(y == 3) -1;
y(y == 4) = 1;
% for 
% end

h = gscatter(X(:,1),X(:,2),y);
hold on
gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
g = gca;
% fcontour(gmPDF,[g.XLim g.YLim])
% title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
legend(h,'A and G','C and T')
xlabel('Dwell Time');
ylabel('Current Drop');
hold off

% yA = reshape(mvnpdf(X,muA,mA),length(x2),length(x1));
% yT = reshape(mvnpdf(X,muT,mT),length(x2),length(x1));
% yC = reshape(mvnpdf(X,muC,mC),length(x2),length(x1));
% yG = reshape(mvnpdf(X,muG,mG),length(x2),length(x1));
%
% fA = @(x,y) mvnpdf([x,y],muA,mA);
% fT = @(x,y) mvnpdf([x,y],muT,mT);
% fC = @(x,y) mvnpdf([x,y],muC,mC);
% fG = @(x,y) mvnpdf([x,y],muG,mG);
%
% hold on
% fcontour(fA,[-0.2 1 -0.1 1.6]);
% fcontour(fT,[-0.2 1 -0.1 1.6]);
% fcontour(fC,[-0.2 1 -0.1 1.6]);
% fcontour(fG,[-0.2 1 -0.1 1.6]);
%
% h = gscatter(x1,x2,yA);

function m = generateMatrix(sigma,theta)
theta_X = theta(1);
theta_Y = theta(2);
%
s = [theta_X^2, 0; 0,theta_Y^2];
B = [cos(sigma), sin(sigma);-sin(sigma), cos(sigma)];
m = B'*s*B;
end