
clear all; close all;
addpath('~/Documents/kam_cryo')
addpath('../SHT')
addpath(genpath('../easyspin-5.2.18'))
% addpath(genpath('../Prol'))
% addpath(genpath('../pswf_t_fast_nfft_v3.3.2'))
addpath(genpath('../pswf_t_mod'))
addpath('~/EMDB_maps')
addpath('../aspire')
initpath
addpath(genpath('../manopt'))
addpath('../nfft')

maxL = 5;
L = 20;
r_cut = 1/2;
r_select_ratio = 1;

s_lens = gen_s_list(maxL, r_cut, r_select_ratio, floor(L/2));
[L_list, m_list, s_list, p_inds, n_inds] = gen_vec_coeff_lists(maxL, s_lens);

info.mol = 'emd_8117';
info.maxL = maxL; % truncation for spherical harmonics expansion
info.L0 = L; % size of volume
info.r_cut = r_cut; % bandlimit
info.N=floor(info.L0/2);

[Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(info);
vol = double(ReadMRC([info.mol '.map'])); %read volume
vol = cryo_downsample(vol, info.L0); 
[a_lms, vol_rec] = expand_vol_spherical_basis(vol, info, Psilms, jball);
vol_rec = real(icfftn(vol_rec));

[~, PSWF_quad_int, points_inside_the_circle] = precomp_pswf_t(L-1, 1, 1e-1, 1);
% Rots = genRotationsGrid( 60, 0, 0);
Rots = randrot(3, 5e4);
% Rots = zeros(3,3,5e4);
% for ii = 1:size(Rots,3), Rots(:,:,ii) = randrot(3); end
SNR = 1/5;
[M1, M2, M3, sigma, gamma] = moments_from_harmonics_projs_PSWFs(Rots, a_lms, L, SNR, PSWF_quad_int, points_inside_the_circle);
% [M1, M2, M3] = moments_from_harmonics_projs(Rots, a_lms, Psilms_2D, jball_2D, L, PSWF_Nn_p, PSWF_quad_int);
M3 = cellfun(@(x) x(:), M3, 'UniformOutput', 0);
M3 = vertcat(M3{:});

% load('/scratch/network/eitanl/B_factors_bispect_maxL5_L20.mat', 'B', 'L1_list', 'L2_list', 'L3_list', 'blk_id', 'W')
W = precomp_wigner_weights(maxL);
tic, [B, L1_list, L2_list, L3_list, blk_id] = precomp_bispectrum_coeffs(maxL, s_lens, L, W); time_precomp = toc;
M3_analytic = bispectrum_from_harmonics_par_noKloop(a_lms, L, B, W, L1_list, L2_list, L3_list, s_lens, blk_id);

B_lists.L1 = L1_list;
B_lists.L2 = L2_list;
B_lists.L3 = L3_list;
B_lists.blk_id = blk_id;

q0 = sum(PSWF_quad_int.ang_freq == 0);
[r,w] = lgwt(20*q0, 0, 1);
w = w.*r;
j_l = generate_spherical_bessel_basis(length(a_lms)-1, s_lens, 1/2, (1/2)*r);

[R_0n, alpha_Nn_2D] = PSWF_radial_2D(0, q0-1, pi*(L-1), r); % generate 2D radial prolates
R_0n = bsxfun(@times, R_0n, 2./alpha_Nn_2D(:).');
a_lms_vec = cellfun(@(x)x(:), a_lms, 'UniformOutput', 0);
a_lms_vec = vertcat(a_lms_vec{:});
M2_analytic = power_spectrum_from_harmonics(a_lms_vec, j_l, R_0n, w, L);

j_0 = cell2mat(generate_spherical_bessel_basis(0, s_lens, 1/2, 0));
M1_analytic = mean_from_harmonics(a_lms_vec, j_0, s_lens(1), L);

a_lms_true = a_lms;
tic, [a_lms_analytic, gamma_analytic] = recover_vol_coeffs_from_moments(M3_analytic, M2_analytic, M1_analytic, maxL, L, B, W, B_lists, r_cut, []); time_rec_analytic = toc
tic, [a_lms, gamma_rec] = recover_vol_coeffs_from_moments(M3, M2, maxL, L, B, W, B_lists, r_cut, []); time_rec = toc
