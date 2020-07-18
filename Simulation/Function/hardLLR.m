function llr = hardLLR(r, P) % P = [Pac,Pat,Pag,Pct,Pcg,Ptg]
Pac = P(1);
Pat = P(2);
Pag = P(3);
Pct = P(4);
Pcg = P(5);
Ptg = P(6);

n = length(r);
llr = zeros(2*n,1);

Paa = 1 - Pac - Pat - Pag;
Pcc = 1 - Pac - Pcg - Pct;
Ptt = 1 - Pat - Ptg - Pct;
Pgg = 1 - Pag - Pcg - Ptg;

for i = 1:n
    if(r(i) == 1)
        llr(i) = log((Pac + Pcc)/(Pcg + Pct));
        llr(i+n) = log((Pac + Pct)/(Pcg + Pcc));
    elseif(r(i) == 2)
        llr(i) = log((Pat + Pct)/(Ptg + Ptt));
        llr(i+n) = log((Pat + Ptt)/(Ptg + Pct));
    elseif(r(i) == 3)
        llr(i) = log((Pag + Pcg)/(Pgg + Ptg));
        llr(i+n) = log((Pag + Ptg)/(Pgg + Pcg));
    elseif(r(i) == 0)
        llr(i) = log((Paa + Pac)/(Pag + Pat));
        llr(i+n) = log((Paa + Pat)/(Pag + Pac));
    end
end
end