function [chan1,chan2,P1,P2,checkMatrix] = pdfmSoft(matrix,mu,acc)
muA = mu(1:2);
muC = mu(3:4);
muT = mu(5:6);
muG = mu(7:8);

mA = matrix(:,1:2);
mC = matrix(:,3:4);
mT = matrix(:,5:6);
mG = matrix(:,7:8);

edge = 30/acc;
p1A = zeros(1,2*edge+1);
p1C = zeros(1,2*edge+1);
p1T = zeros(1,2*edge+1);
p1G = zeros(1,2*edge+1);
p2A = zeros(1,2*edge+1);
p2C = zeros(1,2*edge+1);
p2T = zeros(1,2*edge+1);
p2G = zeros(1,2*edge+1);

p1A_N = zeros(1,2*edge+1);
p1C_N  = zeros(1,2*edge+1);
p1T_N  = zeros(1,2*edge+1);
p1G_N  = zeros(1,2*edge+1);
p2A_N  = zeros(1,2*edge+1);
p2C_N  = zeros(1,2*edge+1);
p2T_N  = zeros(1,2*edge+1);
p2G_N  = zeros(1,2*edge+1);

chan1 = zeros(1,2*edge+1);
chan2 = zeros(1,2*edge+1);
% chan_lu = zeros(1,2*edge+1);
startDwellTime = -1;
endDwellTime = 2;
startCur = -1;
endDwellCur = 2;
sampleNumber = 0;
Lu = 0;
Lu1 = 0;
checkSize = (endDwellTime - startDwellTime)/acc;
checkMatrix = zeros(checkSize+1,checkSize+1);
pro = 0;
% for i = startDwellTime:acc:endDwellTime
%     for j = startCur:acc:endDwellCur
%         yTime = i;
%         yCurrent = j;
%
%         pdfA = mvnpdf([yTime,yCurrent],muA,mA);
%         pdfT = mvnpdf([yTime,yCurrent],muT,mT);
%         pdfC = mvnpdf([yTime,yCurrent],muC,mC);
%         pdfG = mvnpdf([yTime,yCurrent],muG,mG);
%
%         llr1 = log((pdfA + pdfC)/(pdfT+pdfG));
%         llr2 = log((pdfA + pdfT)/(pdfC+pdfG));
%         lr1 = round(llr1/acc) + edge +1;
%         lr1 = bound(lr1, edge);
%         p1A(lr1) =  pdfA;
%         p1C(lr1) =  pdfC;
%         p1T(lr1) =  pdfT;
%         p1G(lr1) =  pdfG;
%     end
% end
counter = 1;

for i = startDwellTime:acc:endDwellTime
    for j = startCur:acc:endDwellCur
        sampleNumber = sampleNumber + 1;
        yTime = i;
        yCurrent = j;
        
        pdfA = mvnpdf([yTime,yCurrent],muA,mA);
        pdfT = mvnpdf([yTime,yCurrent],muT,mT);
        pdfC = mvnpdf([yTime,yCurrent],muC,mC);
        pdfG = mvnpdf([yTime,yCurrent],muG,mG);
        
        llr1 = log((pdfA + pdfC)/(pdfT+pdfG));
        llr2 = log((pdfA + pdfT)/(pdfC+pdfG));
        %         llr_lu = log((pdfA + pdfG)/(pdfC+pdfT));
        
        lr1 = round(llr1/acc) + edge +1;
        lr1 = bound(lr1, edge);
        nlr1 = round(-llr1/acc) + edge +1;
        nlr1 = bound(nlr1, edge);
        
        lr2 = round(llr2/acc) + edge +1;
        lr2 = bound(lr2, edge);
        nlr2 = round(-llr2/acc) + edge +1;
        nlr2 = bound(nlr2, edge);
        
        pro = pro + (pdfA + pdfC + pdfT + pdfG)*acc^2/4;
        chan1(lr1) = chan1(lr1) + (pdfA + pdfC)*acc^2/4;
        chan1(nlr1) = chan1(nlr1) + (pdfT + pdfG)*acc^2/4;
        chan2(lr2) = chan2(lr2) + (pdfA + pdfT)*acc^2/4;
        chan2(nlr2) = chan2(nlr2) + (pdfC + pdfG)*acc^2/4;
        
        p1A(lr1) =  p1A(lr1) + pdfA*acc^2/4;
        p1C(lr1) =  p1C(lr1) +pdfC*acc^2/4;
        p1T(lr1) =  p1T(lr1) +pdfT*acc^2/4;
        p1G(lr1) =  p1G(lr1) +pdfG*acc^2/4;
        
        p1A(nlr1) =  p1A(nlr1) + pdfG*acc^2/4;
        p1C(nlr1) =  p1C(nlr1) + pdfT*acc^2/4;
        p1T(nlr1) =  p1T(nlr1) + pdfC*acc^2/4;
        p1G(nlr1) =  p1G(nlr1) + pdfA*acc^2/4;
        
        p2A(lr2) =  p2A(lr2) + pdfA*acc^2/4;
        p2C(lr2) =  p2C(lr2) + pdfC*acc^2/4;
        p2T(lr2) =  p2T(lr2) + pdfT*acc^2/4;
        p2G(lr2) =  p2G(lr2) + pdfG*acc^2/4;
  
        p2A(nlr2) =  p2A(nlr2) + pdfG*acc^2/4;
        p2C(nlr2) =  p2C(nlr2) + pdfT*acc^2/4;
        p2T(nlr2) =  p2T(nlr2) + pdfC*acc^2/4;
        p2G(nlr2) =  p2G(nlr2) + pdfA*acc^2/4;
        
%         p1A_N(lr1) =  p1A_N(lr1) + 1;
%         p1C_N(lr1) =  p1C_N(lr1) + 1;
%         p1T_N(lr1) =  p1T_N(lr1) + 1;
%         p1G_N(lr1) =  p1G_N(lr1) + 1;
%         
%         p1A_N(nlr1) =  p1A_N(nlr1) + 1;
%         p1C_N(nlr1) =  p1C_N(nlr1) + 1;
%         p1T_N(nlr1) =  p1T_N(nlr1) + 1;
%         p1G_N(nlr1) =  p1G_N(nlr1) + 1;
%         
%         p2A_N(lr2) =  p2A_N(lr2) + 1;
%         p2C_N(lr2) =  p2C_N(lr2) + 1;
%         p2T_N(lr2) =  p2T_N(lr2) + 1;
%         p2G_N(lr2) =  p2G_N(lr2) + 1;
%   
%         p2A_N(nlr2) =  p2A_N(nlr2) + 1;
%         p2C_N(nlr2) =  p2C_N(nlr2) + 1;
%         p2T_N(nlr2) =  p2T_N(nlr2) + 1;
%         p2G_N(nlr2) =  p2G_N(nlr2) + 1;
        
        p1A(lr1) =  pdfA;
        p1C(lr1) =  pdfC;
        p1T(lr1) = pdfT;
        p1G(lr1) =  pdfG;
        
        p1A(nlr1) = pdfG;
        p1C(nlr1) =  pdfT;
        p1T(nlr1) = pdfC;
        p1G(nlr1) =  pdfA;
        
        p2A(lr2) = pdfA;
        p2C(lr2) =   pdfC;
        p2T(lr2) =   pdfT;
        p2G(lr2) =  pdfG;
        %         p2A(nlr2) =  pdfA;
        %         p2C(nlr2) =  pdfC;
        %         p2T(nlr2) =  pdfT;
        %         p2G(nlr2) =  pdfG;
        p2A(nlr2) =   pdfG;
        p2C(nlr2) =  pdfT;
        p2T(nlr2) =  pdfC;
        p2G(nlr2) =  pdfA;
        %         chan_lu(lrlu) = chan_lu(lrlu) + (pdfA + pdfG)*acc^2/4;
        %         chan_lu(nlrlu) = chan_lu(nlrlu) + (pdfC + pdfT)*acc^2/4;
        
        %         if (pdfT > pdfA && pdfT > pdfG && yCurrent < muG(2) &&  yTime < muA(1))
        %             checkMatrix(round(i/acc)+101,round(j/acc)+101) = 1;
        %         else
        %             1;
        %         end
    end
end
% Lu = Lu/sampleNumber;
% Lu1 = Lu1/sampleNumber;
% LuLLR = log(Lu/Lu1);
% p1A = p1A./p1A_N;
% p1C = p1C./p1C_N;
% p1T = p1T./p1T_N;
% p1G = p1G./p1G_N;
% 
% p2A = p2A./p2A_N;
% p2C = p2C./p2C_N;
% p2T = p2T./p2T_N;
% p2G = p2G./p2G_N;

P1 = [p1A; p1C; p1T;p1G];
P2 = [p2A; p2C; p2T;p2G];
end

function llr = bound(llr, edge)
if llr > (2*edge+1)
    llr = 2*edge+1;
elseif llr < 1
    llr = 1;
end
end

% [muA,muC,muT,muG] = matsplit(mu);
% [thetaA,thetaC,thetaT,thetaG ] = matsplit(theta);

% c = normpdf(muA,muC,thetaC);
% a = normpdf(muA,muA,thetaA);
% g = normpdf(muA,muG,thetaG);
% t = normpdf(muA,muT,thetaT);
% llr1 = log((c+a)/(g+t));
% llr2= log((t+a)/(g+c));