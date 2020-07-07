function y = new_xvarEx(func1, num, ext)

sf1 = length(func1);

zeropad1 = zeros(1, sf1*num - num +1);

% direct probabilities
zeropad1(1:sf1) = func1*ext(2);

F_zeropad1 = fft(zeropad1);
foo1 = F_zeropad1.^num;
IF_zeropad = ifft(foo1);

% extract function of the appropriate length from the middle ...
% question: where is the zero?
% clearly, the minimum index is the sum of all input minimum indices

minx = num*ext(1);

% let's hope that state(2) and ext(2) are the same

ext_minx_index = round((ext(1) - minx)/ext(2)) + 1;
ext_maxx_index = ext_minx_index + ext(3) - 1;

%thud = IF_zeropad/sum(IF_zeropad);
ufl = abs(sum(IF_zeropad(1:(ext_minx_index-1))));
ofl = abs(sum(IF_zeropad((1+ext_maxx_index):length(IF_zeropad))));

IF_zeropad(ext_minx_index) = IF_zeropad(ext_minx_index) + ufl;
IF_zeropad(ext_maxx_index) = IF_zeropad(ext_maxx_index) + ofl;

y = abs(IF_zeropad(ext_minx_index:ext_maxx_index));


