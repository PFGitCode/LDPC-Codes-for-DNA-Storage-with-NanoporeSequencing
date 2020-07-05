function llr = hardLLR(r, P) % P = [Pac,Pat,Pag,Pct,Pcg,Ptg]
Pac = P(1);
Pat = P(2);
Pag = P(3);
Pct = P(4);
Pcg = P(5);
Ptg = P(6);

n = length(r);
llr = zeros(2*n,1);
for i = 1:n
    if(r(i) == 1)
        llr(i) = log((1-Pcg-Pct)/(Pcg+Pct));
        llr(i+n) = log((Pac+Pct)/(1-Pac-Pct));
    elseif(r(i) == 2)
        llr(i) = log((Pat+Pct)/(1-Pat-Pct));
        llr(i+n) = log((1-Ptg-Pct)/(Ptg+Pct));
    elseif(r(i) == 3)
        llr(i) = log((Pcg+Pag)/(1-Pcg-Pag));
        llr(i+n) = log((Ptg+Pag)/(1-Ptg-Pag));
    elseif(r(i) == 0)
        llr(i) = log((1-Pat-Pag)/(Pat+Pag));
        llr(i+n) = log((1-Pac-Pag)/(Pac+Pag));
    end
end
end