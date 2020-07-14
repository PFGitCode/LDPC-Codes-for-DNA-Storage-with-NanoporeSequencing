function y = method2Ex(func1,func2, ext, P,bit)
pA = P(1,:);
pC = P(2,:);
pT = P(3,:);
pG = P(4,:);

sf1 = length(func1);
acc = 0.01;
edge = abs(ext(1)/acc);

extrinsic21 = zeros(1,sf1);
for  i = ext(1):ext(2):-ext(1)
    funPos = round(i/acc)+ edge +1;
    sumi2 = i;
    if pA(funPos) == 0
        %         pa = pA(2*edge+2 - funPos);
        pdfA = 0.00000001;
    else
        pdfA = pA(funPos);
    end
    
    if pC(funPos) == 0
        pdfC = 0.00000001;
    else
        pdfC = pC(funPos);
    end
    if pT(funPos) == 0
        pdfT = 0.00000001;
    else
        pdfT = pT(funPos);
    end
    if pG(funPos) == 0
        pdfG = 0.00000001;
    else
        pdfG = pG(funPos);
    end
%     Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
%     Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
    %     Pz1L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp(-(sumi1)/2);
    %     Pz0L1 = exp(-(sumi1)/2)/(1+exp(-(sumi1)))*exp((sumi1)/2);
    
    
    if bit == 1
        Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
        Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
        sumi = log(((pdfA/(pdfA + pdfC))*Pz0L2+(pdfT/(pdfT + pdfG))*Pz1L2)/((pdfC/((pdfC + pdfG)))*Pz0L2+(pdfG/(pdfG+pdfT))*Pz1L2));
    else
        Pz1L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp((sumi2)/2);
        Pz0L2 = exp(-(sumi2)/2)/(1+exp(-(sumi2)))*exp(-(sumi2)/2);
        sumi = log(((pdfA/(pdfA + pdfT))*Pz0L2+(pdfC/(pdfC + pdfG))*Pz1L2)/((pdfT/((pdfT + pdfG)))*Pz0L2+(pdfG/(pdfG+pdfC))*Pz1L2));
    end
    pos = round(sumi/acc)+ edge +1;
    pos = bound(pos, edge);
    extrinsic21(pos) = extrinsic21(pos) + func2(funPos);
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

