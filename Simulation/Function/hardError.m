function P = hardError(mu,matrix)
muA = mu(1:2);
muC = mu(3:4);
muT = mu(5:6);
muG = mu(7:8);

mA = matrix(:,1:2);
mC = matrix(:,3:4);
mT = matrix(:,5:6);
mG = matrix(:,7:8);

Pac = 0;
Pat = 0;
Pag = 0;
Pct = 0;
Pcg = 0;
Ptg = 0;

gap = 0.01;
x1 = -1:gap:1;
x2 = -1:gap:1.6;

for i = 1:length(x1)
    for j = 1:length(x2)
        x = x1(i);
        y = x2(j);
        
        a = mvnpdf([x,y],muA,mA);
        c = mvnpdf([x,y],muC,mC);
        t = mvnpdf([x,y],muT,mT);
        g = mvnpdf([x,y],muG,mG);
        
        if a < c
            Pac = Pac + a * gap^2;
        else
            Pac = Pac + c * gap^2;
        end
        
        if a < t
            Pat = Pat + a * gap^2;
        else
            Pat = Pat + t * gap^2;
        end
        
        if a < g
            Pag = Pag + a * gap^2;
        else
            Pag = Pag + g * gap^2;
        end
        
        if c < t
            Pct = Pct + c * gap^2;
        else
            Pct = Pct + t * gap^2;
        end
        
        if c < g
            Pcg = Pcg + c * gap^2;
        else
            Pcg = Pcg + g * gap^2;
        end
        
        if t < g
            Ptg = Ptg + t * gap^2;
        else
            Ptg = Ptg + g * gap^2;
        end
    end
end
Pac = Pac/2;
Pat = Pat/2;
Pag = Pag/2;
Pct = Pct/2;
Pcg = Pct/2;
Ptg = Ptg/2;

P = [Pac,Pat,Pag,Pct,Pcg,Ptg];
end

