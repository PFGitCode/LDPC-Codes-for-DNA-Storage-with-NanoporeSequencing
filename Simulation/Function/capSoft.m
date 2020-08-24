function I = capSoft(matrix,mu)
muA = mu(1:2);
muC = mu(3:4);
muT = mu(5:6);
muG = mu(7:8);

mA = matrix(:,1:2);
mC = matrix(:,3:4);
mT = matrix(:,5:6);
mG = matrix(:,7:8);

funy = @(yTime,yCurrent) ((mvnpdf([yTime,yCurrent],muA,mA)/4 + mvnpdf([yTime,yCurrent],muC,mC)/4 + ...
mvnpdf([yTime,yCurrent],muT,mT)/4 + mvnpdf([yTime,yCurrent],muG,mG)/4)) *...
log2((mvnpdf([yTime,yCurrent],muA,mA)/4 + mvnpdf([yTime,yCurrent],muC,mC)/4 + ...
mvnpdf([yTime,yCurrent],muT,mT)/4 + mvnpdf([yTime,yCurrent],muG,mG)/4));
% funyy = @(yTime,yCurrent) funy*log2(funy);
% funy = @(yTime,yCurrent) 1/(sqrt(det(mA)*(2*pi)^2))*exp(-0.5*([yTime,yCurrent]-muA)*inv(mA)*([yTime,yCurrent]-muA)');
% format long
q = integral2(@(yTime,yCurrent)arrayfun(funy,yTime,yCurrent),-2,2,-2,2);
hy = -q;
fy = (hyxFun(mA)+hyxFun(mC)+hyxFun(mT)+hyxFun(mG))/4;
hyx = fy;
I = hy - hyx;

% fy = -(hyxFunNew(muA,mA)+hyxFunNew(muC,mC)+hyxFunNew(muT,mT)+hyxFunNew(muG,mG))/4;
% hyx = fy;
% I = hy - hyx;

end

function fy = hyxFun(m)
e = 2.7182818;
fy = 0.5*log2((2*pi*e)^2*det(m));
end

function fxy = hyxFunNew(mu,m)
funyx =  @(yTime,yCurrent) (mvnpdf([yTime,yCurrent],mu,m))*log2((mvnpdf([yTime,yCurrent],mu,m)));
fxy = integral2(@(yTime,yCurrent)arrayfun(funyx,yTime,yCurrent),-0.6,1,-0.2,1.5);
end