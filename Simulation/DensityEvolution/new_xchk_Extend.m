function y = new_xchk_Extend(ext, f_ext1, f_ext2, mapping)
num = 2;
% for some reason, "wrap" = normal log-likelihood ratio
% and "flog" = log(tanh(x/2)), where x is the llr

% excess = [ofl_pos ufl_pos ofl_neg ufl_neg]
[f_log_pos1,f_log_neg1,excess1] = new_wrap2flog(ext, f_ext1, mapping);
[f_log_pos2,f_log_neg2,excess2] = new_wrap2flog(ext, f_ext2, mapping);

p_zero = f_ext1(round((ext(3)-1)/2 + 1))*ext(2);
% n_log contains the new bins for the convolved function

n_log = [mapping(1)*num mapping(2) ((mapping(3)+1)*num - num + 1)];

[f_n_log_pos,f_n_log_neg,p_result_zero] = new_fft_convolve_chk_Extend(f_log_pos1, f_log_neg1, f_log_pos2, f_log_neg2, n_log, n_log, excess1, excess2, p_zero);

[f_n_ext,ofl_pos,ofl_neg] = new_flog2wrap(n_log, f_n_log_pos, f_n_log_neg, ext, p_result_zero);

% note: the overflows should be placed at the value corresponding to
% 2*atanh(exp(-mapping(2)/2)), because this is the uppermost value that
% can be represented using the quantization given by mapping

f_n_ext = new_chk_overflow(f_n_ext,ofl_pos,ofl_neg,ext,mapping); 

y = f_n_ext;
end

