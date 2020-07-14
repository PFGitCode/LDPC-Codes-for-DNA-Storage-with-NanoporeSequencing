function y = method1Ex(func1,func2, ext)
sf1 = length(func1);
acc = 0.01;
edge = abs(ext(1)/acc);
alpha = 0.1;

extrinsic21 = zeros(1,sf1);
for  i = ext(1):ext(2):-ext(1)
    funPos = round(i/acc)+ edge +1;
    sumi2 = i;
    sumi = alpha*sumi2;
    pos = round(sumi/acc)+ edge +1;
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

