function [chan1,chan2,chan_lu,checkMatrix] = pdfmSoft(matrix,mu,acc)
muA = mu(1:2);
muC = mu(3:4);
muT = mu(5:6);
muG = mu(7:8);

mA = matrix(:,1:2);
mC = matrix(:,3:4);
mT = matrix(:,5:6);
mG = matrix(:,7:8);

edge = 30/acc;
% chan1A = zeros(1,2*edge+1);
% chan2A = zeros(1,2*edge+1);
% chan1C = zeros(1,2*edge+1);
% chan2C = zeros(1,2*edge+1);
% chan1T = zeros(1,2*edge+1);
% chan2T = zeros(1,2*edge+1);
% chan1G = zeros(1,2*edge+1);
% chan2G = zeros(1,2*edge+1);
chan1 = zeros(1,2*edge+1);
chan2 = zeros(1,2*edge+1);
chan_lu = zeros(1,2*edge+1);
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
        llr_lu = log((pdfA + pdfG)/(pdfC+pdfT));
        
        lr1 = round(llr1/acc) + edge +1;
        lr1 = bound(lr1, edge);
        nlr1 = round(-llr1/acc) + edge +1;
        nlr1 = bound(nlr1, edge);
        
        lr2 = round(llr2/acc) + edge +1;
        lr2 = bound(lr2, edge);
        nlr2 = round(-llr2/acc) + edge +1;
        nlr2 = bound(nlr2, edge);
        
        lrlu = round(llr_lu/acc) + edge +1;
        lrlu = bound(lrlu, edge);
        nlrlu = round(-llr_lu/acc) + edge +1;
        nlrlu = bound(nlrlu, edge);
        %         if(abs(lr1) > 2*edge+1 || lr1 < 1 || abs(lr2) > 2*edge+1 || lr2 < 1)
        %             continue;
        %         end
        pro = pro + (pdfA + pdfC + pdfT + pdfG)*acc^2/4;
        chan1(lr1) = chan1(lr1) + (pdfA + pdfC)*acc^2/4;
        chan1(nlr1) = chan1(nlr1) + (pdfT + pdfG)*acc^2/4;
        chan2(lr2) = chan2(lr2) + (pdfA + pdfT)*acc^2/4;
        chan2(nlr2) = chan2(nlr2) + (pdfC + pdfG)*acc^2/4;
        chan_lu(lrlu) = chan_lu(lrlu) + (pdfA + pdfG)*acc^2/4;
        chan_lu(nlrlu) = chan_lu(nlrlu) + (pdfC + pdfT)*acc^2/4;
        
%         if (pdfT > pdfA && pdfT > pdfG && yCurrent < muG(2) &&  yTime < muA(1))
%             checkMatrix(round(i/acc)+101,round(j/acc)+101) = 1;
%         else
%             1;
%         end
    end
end
Lu = Lu/sampleNumber;
Lu1 = Lu1/sampleNumber;
LuLLR = log(Lu/Lu1);
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