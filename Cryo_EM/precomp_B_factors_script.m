addpath('../kam_cryo')
addpath('../SHT')
addpath(genpath('../pswf_t_mod'))
addpath(genpath('../easyspin-5.2.18'))
addpath('../aspire')
initpath
addpath('../nfft')

maxL = 5; L = 31; r_cut = 1/2; r_select_ratio = 1;

s_list = gen_s_list(maxL, r_cut, r_select_ratio, floor(L/2));
W = precomp_wigner_weights(maxL);
parpool('local', maxNumCompThreads)
[B, L1_list, L2_list, L3_list, blk_id] = precomp_bispectrum_coeffs_matMult(maxL, s_list, L, W);

save('/scratch/network/eitanl/B_factors_bispect_maxL5_L31', '-v7.3')
