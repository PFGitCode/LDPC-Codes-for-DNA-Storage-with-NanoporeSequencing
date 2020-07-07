function [f_n_log_pos,f_n_log_neg,p_result_zero] = new_fft_convolve_chk_Extend(f_log_pos1, f_log_neg1, f_log_pos2, f_log_neg2, n_log1, n_log2, excess1, excess2, p_zero);
num = 1;
% perform FFTs over R x GF(2) to obtain probabilities with sign information

% I will *assume* that the range of f_log_pos and f_log_neg are over the range
% min:increment:-increment, with min*increment elements
%
% the input ofl_pos is the probability of the input "flog" message
% being zero ... can directly include this in the sum
%
% furthermore, by the symmetry property, "flog" can only be zero
% if the sign is positive

% find the probabilities of being positive or negative ... used in calculating the
% conditional probability

% these form the extended (by 1) messages, which give room for an extra zero message
q_pos1 = zeros(1,length(f_log_pos1)+1);
q_pos1(1:length(f_log_pos1)) = f_log_pos1;
q_neg1 = zeros(1,length(f_log_neg1)+1);
q_neg1(1:length(f_log_neg1)) = f_log_neg1;

q_pos2 = zeros(1,length(f_log_pos2)+1);
q_pos2(1:length(f_log_pos2)) = f_log_pos2;
q_neg2 = zeros(1,length(f_log_neg2)+1);
q_neg2(1:length(f_log_neg2)) = f_log_neg2;
% sf = size of the extended message
sf = length(f_log_pos1)+1;

% find length to nearest higher power of 2 to help with fft
zp2_size = 2^(ceil(log2(sf*num - num + 1)));
zeropad_pos1 = zeros(1, zp2_size);
zeropad_neg1 = zeropad_pos1;

zeropad_pos2 = zeros(1, zp2_size);
zeropad_neg2 = zeropad_pos2;

% preserve *direct* probabilities under convolution
% also note that overflow is included in zeropad_pos
zeropad_pos1(1:sf) = q_pos1*n_log1(2);
zeropad_neg1(1:sf) = q_neg1*n_log1(2);
zeropad_pos1(sf) = excess1(1);
zeropad_neg1(sf) = excess1(3);

zeropad_pos2(1:sf) = q_pos2*n_log2(2);
zeropad_neg2(1:sf) = q_neg2*n_log2(2);
zeropad_pos2(sf) = excess2(1);
zeropad_neg2(sf) = excess2(3);

% find the probabilities of being positive or negative and finite ...
% as well as marginal and conditional probabilities
p_pos_fin1 = sum(zeropad_pos1);
p_pos1 = sum(zeropad_pos1) + excess1(2);

p_pos_fin2 = sum(zeropad_pos2);
p_pos2 = sum(zeropad_pos2) + excess2(2);
% note: excess(1), excess(3) are already contained in zeropad

p_neg_fin1 = sum(zeropad_neg1);
p_neg1 = sum(zeropad_neg1) + excess1(4);

p_neg_fin2 = sum(zeropad_neg2);
p_neg2 = sum(zeropad_neg2) + excess2(4);
% fprintf('%f \n ',zeropad_pos/p_pos_fin);
% fprintf('%f \n ',zeropad_neg/p_pos_fin);
% here we take the FFT of the *conditional* density
F_pos_fin1 = fft(zeropad_pos1/p_pos_fin1);
F_neg_fin1 = fft(zeropad_neg1/p_neg_fin1);

F_pos_fin2 = fft(zeropad_pos2/p_pos_fin2);
F_neg_fin2 = fft(zeropad_neg2/p_neg_fin2);

result_pos = zeros(1, zp2_size);
result_neg = zeros(1, zp2_size);
p_result_zero = 0;

% now combine the magnitude and sign components to get the result
% have to do proper handing for -\infty messages (i.e., wrap=0)

foo = zeros(num+1,zp2_size);

% now marginalize with respect to positive or negative events
% (i.e., c even or odd)
%
% note that result_pos and result_neg are joint probabilities
% of the amplitude and being positive or negative

for c = 0:3
    if c == 0
        v = (1 - (p_pos1/(p_pos1+p_zero))*(p_pos2/(p_pos2+p_zero))) * (p_pos1+p_zero)*(p_pos2+p_zero);
        p_result_zero = p_result_zero + v;    
        temp2=(1-excess1(2))*(1-excess2(2));
        temp3= p_pos1*p_pos2;
        v = (F_pos_fin1.*F_pos_fin2);
    elseif c == 1
        v = (1 - (p_pos1/(p_pos1+p_zero))) * p_neg2 * (p_pos1+p_zero);
        p_result_zero = p_result_zero + v;  
        temp2=(1-excess2(4))*(1-excess1(2));
        temp3= p_pos1*p_neg2;
        v = (F_neg_fin1.*F_pos_fin2);
    elseif c == 3
        v = (1 - (p_pos2/(p_pos2+p_zero))) * p_neg1 * (p_pos2+p_zero);
        p_result_zero = p_result_zero + v;  
        temp2=(1-excess1(4))*(1-excess2(2));
        temp3= p_neg1*p_pos2;
        v = (F_neg_fin1.*F_pos_fin2);
    elseif c == 2
        v =  p_neg1 * p_neg2;
        p_result_zero = p_result_zero + v;  
        temp2=(1-excess1(4))*(1-excess2(4));
        temp3= p_neg2*p_neg1;
        v = (F_neg_fin1.*F_neg_fin2);
    end
    %     v = (F_pos_fin1.^(num-c)).*(F_neg_fin1.^c);
    v = v * temp2;
    v = v * temp3;
    
    if (mod(c,2)==0)
        % even negatives -- result is positive
        result_pos = result_pos + v;
    else
        % odd negatives -- result is negative
        result_neg = result_neg + v;
    end
    
    % at least one underflow
    % regardless of sign, attached to p_result_zero
    
    % if there are c negative messages, then there are at most num-c positive messages
%     
%     if c==0
%         for d = 0:(num-c)
%             % calculate no underflow first ... then complement
%             w = (1-excess1(4))^d;
%             v = (1-w)*prod(1:num)/(prod(1:c)*prod(1:d)*prod(1:(num-c-d)));
%             v = v * p_neg1^c * p_pos1^d * p_zero^(num-c-d);
%             
%             p_result_zero = p_result_zero + v;
%             
%         end
%     end
    for d = 0:(num-c)
        
        % calculate no underflow first ... then complement
        w = (1-excess1(2))^c * (1-excess1(4))^d;
        v = (1-w)*prod(1:num)/(prod(1:c)*prod(1:d)*prod(1:(num-c-d)));
        v = v * p_neg1^c * p_pos1^d * p_zero^(num-c-d);
        
        p_result_zero = p_result_zero + v;
        
    end
    
end

result_pos=ifft(result_pos);
result_neg=ifft(result_neg);
result_pos = abs(result_pos)/n_log1(2);
result_neg = abs(result_neg)/n_log1(2);

f_n_log_pos = result_pos(1:(sf*num - num + 1));
f_n_log_neg = result_neg(1:(sf*num - num + 1));


