function data = hardChannel(data, P) % P = [Pac,Pat,Pag,Pct,Pcg,Ptg]
Pac = P(1);
Pat = P(2);
Pag = P(3);
Pct = P(4);
Pcg = P(5);
Ptg = P(6);

for i = 1:length(data)
    if data(i) == 0
        data(i) = randsrc(1,1,[0,1,2,3;1-Pac-Pat-Pag,Pac,Pat,Pag]);
    elseif data(i) == 1
        data(i) = randsrc(1,1,[0,1,2,3;Pac,1-Pac-Pcg-Pct,Pct,Pcg]);
    elseif data(i) == 2
        data(i) = randsrc(1,1,[0,1,2,3;Pat,Pct,1-Pat-Pct-Ptg,Ptg]);
    elseif data(i) == 3
        data(i) = randsrc(1,1,[0,1,2,3;Pag,Pcg,Ptg,1-Ptg-Pcg-Pag]);
    end
end

end