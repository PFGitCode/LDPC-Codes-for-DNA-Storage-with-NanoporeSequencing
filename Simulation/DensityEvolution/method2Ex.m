function y = method2Ex(func1,func2, ext, P,bit)
pA = P(1,:);
pC = P(2,:);
pT = P(3,:);
pG = P(4,:);

sf1 = length(func1);
acc = 0.01;
edge = abs(ext(1)/acc);

extrinsic21 = zeros(1,sf1);
iter = 1;
counter = 1;
for  i = ext(1):ext(2):-ext(1)
    funPos = round(i/acc)+ edge +1;
    %     if pA(funPos) == 0 || pC(funPos) == 0 || pT(funPos) == 0 || pG(funPos) == 0
    %         pos = edge+1;
    %         extrinsic21(pos) = extrinsic21(pos) + func2(funPos);
    %         continue;
    %     end
    %     if pA(funPos) == 0
    %         pdfA = 0.00000001;
    %     else
    pdfA = pA(funPos);
    %     end
    %
    %     if pC(funPos) == 0
    %         pdfC = 0.00000001;
    %     else
    pdfC = pC(funPos);
    %     end
    %     if pT(funPos) == 0
    %         pdfT = 0.00000001;
    %     else
    pdfT = pT(funPos);
    %     end
    %     if pG(funPos) == 0
    %         pdfG = 0.00000001;
    %     else
    pdfG = pG(funPos);
    %     end
    sumi2 = -i;
    Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
    Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
    if (pdfT > pdfA && pdfT > pdfG && pdfC > pdfA && pdfC > pdfG)
        if bit == 2
            sumi(iter) = log((pdfA/(pdfA+pdfT)*Pz0L2+(pdfC/(pdfC+pdfG))*Pz1L2)/(pdfT/(pdfA+pdfT)*Pz0L2+(pdfG/(pdfC+pdfG))*Pz1L2));
        elseif bit == 1
            sumi(iter) = log((pdfA/(pdfA+pdfC)*Pz0L2+(pdfT/(pdfG+pdfT))*Pz1L2)/(pdfC/(pdfC+pdfA)*Pz0L2+(pdfG/(pdfG+pdfT))*Pz1L2));
        end
        pos = round(sumi(iter)/acc)+ edge +1;
        pos = bound(pos, edge);
        extrinsic21(pos) = extrinsic21(pos) + func2(funPos);
        iter = iter + 1;
%         counter = counter+1
    else
        pos = edge+1;
        extrinsic21(pos) = extrinsic21(pos) + func2(funPos);
        sumi(iter) = sumi2;
        iter = iter + 1;
        continue;
    end
    
end
zeropad1 = zeros(1,sf1 + sf1);
zeropad2 = zeropad1;

zeropad1(1:sf1) = func1;
zeropad2(1:sf1) = extrinsic21;
% ext2 = (-(sf2new-1)/2)*acc;
ext2 = ext(1);
F_zeropad1 = fft(zeropad1);
F_zeropad2 = fft(zeropad2);
foo1 = F_zeropad1;
foo2 = F_zeropad2;

foo = foo1 .* foo2;
IF_zeropad = ifft(foo);

% extract function of the appropriate length from the middle ...
% question: where is the zero?
% clearly, the minimum index is the sum of all input minimum indices

minx = ext(1) + ext2;

% let's hope that state(2) and ext(2) are the same

ext_minx_index = round((ext(1) - minx)/ext(2)) + 1;
ext_maxx_index = ext_minx_index + ext(3) - 1;

%thud = IF_zeropad/sum(IF_zeropad);
ufl = abs(sum(IF_zeropad(1:(ext_minx_index-1))));
ofl = abs(sum(IF_zeropad((1+ext_maxx_index):length(IF_zeropad))));

IF_zeropad(ext_minx_index) = IF_zeropad(ext_minx_index) + ufl;
IF_zeropad(ext_maxx_index) = IF_zeropad(ext_maxx_index) + ofl;

y = abs(IF_zeropad(ext_minx_index:ext_maxx_index));
end
function llr = bound(llr, edge)
if llr > (2*edge+1)
    llr = 2*edge+1;
elseif llr < 1
    llr = 1;
end
end
%

% function y = method2Ex(func1,func2, ext, P,bit)
% pdfA = P(1,:);
% pdfC = P(2,:);
% pdfT = P(3,:);
% pdfG = P(4,:);
%
% sf1 = length(func1);
% acc = 0.01;
% edge = abs(ext(1)/acc);
%
% extrinsic21 = zeros(1,sf1);
% iter = 1;
% for  i = ext(1):ext(2):-ext(1)
%     if i == 0
%         1;
%     end
%     funPos = round(i/acc)+ edge +1;
%     sumi2 = i;
% %     if pA(funPos) == 0
% %         %         pa = pA(2*edge+2 - funPos);
% %         pa = 0.00000001;
% %     else
% %         pa = pA(funPos);
% %     end
% %
% %     if pC(funPos) == 0
% %         pc = 0.00000001;
% %     else
% %         pc = pC(funPos);
% %     end
% %     if pT(funPos) == 0
% %         pt = 0.00000001;
% %     else
% %         pt = pT(funPos);
% %     end
% %     if pG(funPos) == 0
% %         pg = 0.00000001;
% %     else
% %         pg = pG(funPos);
% %     end
%     Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
%     Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
% %     Pz1L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp(-(sumi1)/2);
% %     Pz0L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp((sumi1)/2);
%
%
% %     if bit == 1
%         sumi(iter) = log(((pdfA(funPos)/(pdfA(funPos) + pdfC(funPos)))*Pz0L2+(pdfT(funPos)/(pdfT(funPos) + pdfG(funPos)))*Pz1L2)/((pdfC(funPos)/((pdfC(funPos) + pdfA(funPos))))*Pz0L2+(pdfG(funPos)/(pdfG(funPos)+pdfT(funPos)))*Pz1L2));
% %     else
% %         sumi(iter)  = log(((pdfA(funPos)/(pdfA(funPos) + pdfT(funPos)))*Pz0L2+(pdfC(funPos)/(pdfC(funPos) + pdfG(funPos)))*Pz1L2)/((pdfT(funPos)/((pdfT(funPos) + pdfA(funPos))))*Pz0L2+(pdfG(funPos)/(pdfG(funPos)+pdfC(funPos)))*Pz1L2));
% %     end
%     pos = round(-sumi(iter)/acc)+ edge +1;
%     pos = bound(pos, edge);
%     extrinsic21(pos) = extrinsic21(pos) + func2(funPos);
%     iter = iter + 1;
% end
% % extrinsic21 = extrinsic21*0.5;
% % extrinsic21(edge+1) = extrinsic21(edge+1) + 0.5;
% zeropad1 = zeros(1,sf1 + sf1);
% zeropad2 = zeropad1;
%
% zeropad1(1:sf1) = func1;
% zeropad2(1:sf1) = extrinsic21;
% % ext2 = (-(sf2new-1)/2)*acc;
% ext2 = ext(1);
% F_zeropad1 = fft(zeropad1);
% F_zeropad2 = fft(zeropad2);
% foo1 = F_zeropad1;
% foo2 = F_zeropad2;
%
% foo = foo1 .* foo2;
% IF_zeropad = ifft(foo);
%
% % extract function of the appropriate length from the middle ...
% % question: where is the zero?
% % clearly, the minimum index is the sum of all input minimum indices
%
% minx = ext(1) + ext2;
%
% % let's hope that state(2) and ext(2) are the same
%
% ext_minx_index = round((ext(1) - minx)/ext(2)) + 1;
% ext_maxx_index = ext_minx_index + ext(3) - 1;
%
% %thud = IF_zeropad/sum(IF_zeropad);
% ufl = abs(sum(IF_zeropad(1:(ext_minx_index-1))));
% ofl = abs(sum(IF_zeropad((1+ext_maxx_index):length(IF_zeropad))));
%
% IF_zeropad(ext_minx_index) = IF_zeropad(ext_minx_index) + ufl;
% IF_zeropad(ext_maxx_index) = IF_zeropad(ext_maxx_index) + ofl;
%
% y = abs(IF_zeropad(ext_minx_index:ext_maxx_index));
% end
% function llr = bound(llr, edge)
% if llr > (2*edge+1)
%     llr = 2*edge+1;
% elseif llr < 1
%     llr = 1;
% end
% end
