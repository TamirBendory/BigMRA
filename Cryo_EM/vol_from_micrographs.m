%% Set paths:
clear all; close all;
addpath('~/Documents/kam_cryo')
addpath('../SHT')
addpath('../aspire')
initpath
addpath('~/EMDB_maps')
addpath(genpath('../pswf_t_mod'))
addpath(genpath('../nfft'))
% addpath(genpath('../manopt'))

%% Set parameters, load micrograph
maxL = 5;
r_cut = 1/2;
SNR = 1/128;
seed = 10; % seed for noise generation

load('/scratch/network/eitanl/1qlq_crop_micros_100_M_7420.mat') % contains: I, L, M, gamma, num_micros, projs_per_micrograph, vol

% Add noise:
[I, ~, ~, sigma_true] = cryo_addnoise(I, SNR, 'gaussian', seed);

%% Recover volume:
if isempty(gcp('nocreate')), parpool('local', maxNumCompThreads); end
[SNR_est, ~, var_n]=cryo_estimate_snr(cfft2(I)/M); % estimate SNR and noise variance from corners
sigma = sqrt(var_n);

% Compute moments and debias:
% batch_size = 5e5;
% [m1, m2, m3] = moments_from_micrograph_steerable_parMicro(I, L, batch_size);
[m1, m2, m3] = moments_from_micrograph_steerable_windows(I, L);
[m2, m3] = debias_moments_steerable(m1, m2, m3, sigma, L);
m3 = cellfun(@(x) x(:), m3, 'UniformOutput', 0); 
m3 = vertcat(m3{:});

% Load precomputed factors for bispectrum:
load('/scratch/network/eitanl/B_factors_bispect_maxL5_L31.mat', 'B', 'L1_list', 'L2_list', 'L3_list', 'blk_id', 'W')

B_lists.L1 = L1_list;
B_lists.L2 = L2_list;
B_lists.L3 = L3_list;
B_lists.blk_id = blk_id;

% Optimize for volume expansion coefficients:
[a_lms, gamma_rec] = recover_vol_coeffs_from_moments(m3, m2, m1, maxL, L, B, W, B_lists, r_cut, []);

%% Compute true expansion coefficients:
info.maxL = maxL; % truncation for spherical harmonics expansion
info.L0 = L; % size of volume
info.r_cut = 1/2; % bandlimit
info.N=floor(info.L0/2);
info.pixA = 1+2/3; 

[Psilms, ~, jball, ~] = precompute_spherical_basis(info);
[a_lms_true, vol_true_trunc] = expand_vol_spherical_basis(vol, info, Psilms, jball);
vol_true_trunc = real(icfftn(vol_true_trunc));

%% Align to ground truth, compute FSC:
trials_align = 1e5; % guesses for rotation between volumes
par_flag = 1; % parallelize guesses
[a_lms_aligned, R_rec, cost_align] = align_vol_coeffs(a_lms, a_lms_true, trials_align, par_flag);
vol_rec_aligned = recover_from_ALM_v4_given_Psilms(a_lms_aligned, info.N, jball, Psilms, info.maxL, info.L0);

% Compute FSC, plot:
cutoff_FSC = 0.5; % use 0.5 criterion to compare to ground truth
[res_trunc, fig_trunc]=plotFSC(vol_true_trunc, vol_rec_aligned, cutoff_FSC, info.pixA);
[res_full, fig_full]=plotFSC(vol, vol_rec_aligned, cutoff_FSC, info.pixA);

%% Save the results:

% Save FSC plots:
savefig(fig_trunc, '1qlq_trunc_FSC.fig')
savefig(fig_full, '1qlq_full_FSC.fig')

% Save vars:
save('1qlq_micro_rec_results.mat', 'gamma_rec', 'a_lms_aligned', 'res_trunc', 'res_full', 'm1', 'm2', 'm3')

