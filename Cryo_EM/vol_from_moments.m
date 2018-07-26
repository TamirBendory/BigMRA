clear all; close all;
addpath('~/Documents/kam_cryo')
addpath('../SHT')
addpath('../aspire')
initpath
addpath('~/EMDB_maps')
addpath(genpath('../pswf_t_mod'))
addpath(genpath('../nfft'))
% addpath(genpath('../manopt'))

load('/scratch/eitanl/pure_noise_test_res.mat', 'm1', 'm2', 'm3', 'maxL', 'L', 'r_cut')

load('/scratch/eitanl/B_factors_bispect_maxL5_L20.mat', 'B', 'L1_list', 'L2_list', 'L3_list', 'blk_id', 'W')

B_lists.L1 = L1_list;
B_lists.L2 = L2_list;
B_lists.L3 = L3_list;
B_lists.blk_id = blk_id;

if isempty(gcp('nocreate')), parpool(maxNumCompThreads); end
[a_lms, gamma] = recover_vol_coeffs_from_moments(m3, m2, m1, maxL, L, B, W, B_lists, r_cut, []);

info.maxL = maxL; % truncation for spherical harmonics expansion
info.L0 = L; % size of volume
info.r_cut = 1/2; % bandlimit
info.N=floor(info.L0/2);

[Psilms, ~, jball, ~] = precompute_spherical_basis(info);
vol_rec = recover_from_ALM_v4_given_Psilms(a_lms, info.N, jball, Psilms, maxL, L);
vol_rec = real(icfftn(vol_rec));

clear I B L1_list L2_list L3_list blk_id W
save('/scratch/eitanl/pure_noise_test_res_vol_from_moms.mat', '-v7.3')
