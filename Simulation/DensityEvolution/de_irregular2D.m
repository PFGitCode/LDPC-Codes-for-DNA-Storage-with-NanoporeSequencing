function result_pe1 = de_irregular2D(chan1,chan2,iter,ext,mapping,stop_pe,vard,chkd,dv,dc,LuLLR)
z1 = chan1;
z2 = chan2;
% Node average for the extrisic info
sumV = 0;
sumC = 0;
for i = 1:length(vard)
    sumV = sumV + vard(i)/i;
end
for i = 1:length(chkd)
    sumC = sumC + chkd(i)/i;
end
vardN = (vard./[1:length(vard)])/sumV;

c = 0;
pe1 = 0.5;
pe2 = 0.5;
result_pe1 = zeros(1,iter);
result_pe2 = zeros(1,iter);
round(0.3);
while ((c < iter) & (pe1 > stop_pe)& (pe2 > stop_pe))
    c = c + 1;
    y_ave1=0;
    y_ave2=0;
    for i=2:dc
        if chkd(i)~=0
            y1 = new_xchk(ext, z1, i-1, mapping);
            y2 = new_xchk(ext, z2, i-1, mapping);
            y1 = y1 / (sum(y1)*ext(2));
            y2 = y2 / (sum(y2)*ext(2));
            y_ave1=y_ave1+chkd(i)*y1;
            y_ave2=y_ave2+chkd(i)*y2;
        end
    end
    
    %extrinsic infor
    xvar_aveEX1=0;
    xvar_aveEX2=0;
    for j=1:dv
        if vardN(j)~=0
            zEx1 = new_xvarEx(y_ave1, j, ext);
            zEx2 = new_xvarEx(y_ave2, j, ext);
            xvar_aveEX1=xvar_aveEX1+vardN(j)*zEx1;
            xvar_aveEX2=xvar_aveEX2+vardN(j)*zEx2;
        end
    end
    xvar_ave1=0;
    xvar_ave2=0;
    for j=2:dv
        if vard(j)~=0
            z1 = new_xvar(chan1, y_ave1,j-1, ext, ext);
            z2 = new_xvar(chan2, y_ave2,j-1, ext, ext);
            xvar_ave1=xvar_ave1+vard(j)*z1;
            xvar_ave2=xvar_ave2+vard(j)*z2;
        end
    end
    if decodeMethod == "method3"
        z1 = method3Ex(xvar_ave1,xvar_aveEX2,ext,LuLLR);
        z2 = method3Ex(xvar_ave2,xvar_aveEX1,ext,LuLLR);
    elseif decodeMethod == "method4"
        z1 = method4Ex(xvar_ave1,xvar_aveEX2,ext,LuLLR);
        z2 = method4Ex(xvar_ave2,xvar_aveEX1,ext,LuLLR);
    end
    %     z1=xvar_ave1;
    %     z2=xvar_ave2;
    pe1 = sum(z1(1:round((ext(3)-1)/2 )))*ext(2);
    pe2 = sum(z2(1:round((ext(3)-1)/2 )))*ext(2);
    result_pe1(c) = pe1;
    result_pe2(c) = pe2;
    fprintf('%.10f\n',pe1);
    fprintf('%.10f\n',pe2);
end
