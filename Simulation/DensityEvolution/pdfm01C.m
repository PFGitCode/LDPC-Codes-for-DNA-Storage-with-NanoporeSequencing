function [chan1,chan2] = pdfm01C(mu,theta,acc,Pct,Pat,Pcg,Pag,Ptg,Pac)
[muA,muC,muT,muG] = matsplit(mu);
[thetaA,thetaC,thetaT,thetaG ] = matsplit(theta);
edge = 30/acc;
chan1 = zeros(1,2*edge+1);
chan2 = zeros(1,2*edge+1);
llr1 = zeros(1,140);
llr2 = zeros(1,140);
i = 1;
UsedThrehold = muT + thetaT;
for y = 0.1:0.01:1.3
    c = normpdf(y,muC,thetaC);
    a = normpdf(y,muA,thetaA);
    g = normpdf(y,muG,thetaG);
    t = normpdf(y,muT,thetaT);
    cp = normcdf(y,muC,thetaC);
    ap = normcdf(y,muA,thetaA);
    gp = normcdf(y,muG,thetaG);
    tp = normcdf(y,muT,thetaT);
    llr1(i) = log((c+a)/(g+t));
    llr2(i)= log((t+a)/(g+c));
    if llr1(i) < -30
        llr1(i) = -30;
    elseif llr1(i) > 30
        llr1(i) = 30;
    end
    llr2(i)= log((t+a)/(g+c));
    if llr2(i) < -30
        llr2(i) = -30;
    elseif llr2(i) > 30
        llr2(i) = 30;
    end
    chan1_index = round(llr1(i)/acc) + edge +1;
    nchan1_index = round(-llr1(i)/acc) + edge +1;
    chan1_index_pro = a;
%     nchan1_index_pro = t+g;
    
    chan2_index = round(llr2(i)/acc) + edge +1;
    nchan2_index = round(-llr2(i)/acc) + edge +1;
    chan2_index_pro = t+a;
    nchan2_index_pro = c+g;
    chan1(chan1_index) = chan1_index_pro*0.01;
%     chan1(nchan1_index) = nchan1_index_pro*0.01/4;
    chan1(chan2_index) = chan2_index_pro*0.01/4;
    chan1(nchan2_index) = nchan2_index_pro*0.01/4;
    i = i+1;
end






% when y = A
llrA1 = log((1-Pat-Pag)/(Pat+Pag));
llrA2 = log((1-Pac-Pag)/(Pac+Pag));
%set positive and negative position for A
lA1 = round(llrA1/acc) + edge +1;
nlA1 = round(-llrA1/acc) + edge +1;
lA2 = round(llrA2/acc) + edge +1;
nlA2 = round(-llrA2/acc) + edge +1;
% when y = C
llrC1 = log((1-Pcg-Pct)/(Pcg+Pct));
llrC2 = log((Pac+Pct)/(1-Pac-Pct));
%set positive and negative position for C
lC1 = round(llrC1/acc) + edge +1;
nlC1 = round(-llrC1/acc) + edge +1;
lC2 = round(llrC2/acc) + edge +1;
nlC2 = round(-llrC2/acc) + edge +1;
% when y = T
llrT1 = log((Pat+Pct)/(1-Pat-Pct));
llrT2 = log((1-Pag-Pct)/(Pag+Pct));
%set positive and negative position for T
lT1 = round(llrT1/acc) + edge +1;
nlT1 = round(-llrT1/acc) + edge +1;
lT2 = round(llrT2/acc) + edge +1;
nlT2 = round(-llrT2/acc) + edge +1;
% when y = G
llrG1 = log((Pcg+Pag)/(1-Pcg-Pag));
llrG2 = log((Ptg+Pag)/(1-Ptg-Pag));
%set positive and negative position for T
lG1 = round(llrG1/acc) + edge +1;
nlG1 = round(-llrG1/acc) + edge +1;
lG2 = round(llrG2/acc) + edge +1;
nlG2 = round(-llrG2/acc) + edge +1;

%when x = A
chan1A(lA1) = (1 - Pac - Pag - Pat) + Pac;%A or C for channel1 positive
chan1A(nlA1) = Pag + Pat; %T or G for channel1 negative
chan2A(lA2) = (1 - Pac - Pag - Pat) + Pat; %A or T for channel1 positive
chan2A(nlA2) = Pac + Pag;%C or G for channel1 positive
%when x = C
chan1C(lC1) = (1 - Pac - Pcg - Pct) + Pac;
chan1C(nlC1) = Pct + Pcg;
chan2C(lC2) = Pac + Pct;
chan2C(nlC2) = (1 - Pac - Pcg - Pct) + Pcg;
%when x = T
chan1T(lT1) = Pat + Pct;
chan1T(nlT1) = (1 - Pct - Ptg - Pat) + Ptg;
chan2T(lT2) = Pat + (1 - Pct - Ptg - Pat);
chan2T(nlT2) = Pct + Ptg;
%when x = G
chan1G(lG1) = Pag + Pcg;
chan1G(nlG1) = Ptg + (1 - Pcg - Pag - Ptg);
chan2G(lG2) = Pag + Ptg;
chan2G(nlG2) = Pcg + (1 - Pcg - Pag - Ptg);

chan1 = (chan1A + chan1C + chan1T + chan1G)/4;
chan2 = (chan2A + chan2C + chan2T + chan2G)/4;
end


% [muA,muC,muT,muG] = matsplit(mu);
% [thetaA,thetaC,thetaT,thetaG ] = matsplit(theta);

% c = normpdf(muA,muC,thetaC);
% a = normpdf(muA,muA,thetaA);
% g = normpdf(muA,muG,thetaG);
% t = normpdf(muA,muT,thetaT);
% llr1 = log((c+a)/(g+t));
% llr2= log((t+a)/(g+c));