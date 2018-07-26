function [M1, M2, M3, sigma, gamma] = moments_from_harmonics_projs_PSWFs(Rots, a_lms, L, SNR, PSWF_quad_int, points_inside_the_circle)
% Evaluate the three first moments of the dataset on a L x L cartesian grid

R = floor(L/2);
maxL = length(a_lms)-1;
[psi_Nn, n_list] = PSWF_2D_full_cart(maxL, R);
s_list = cellfun('size', a_lms, 1);
beta = sph_Bessel_to_2D_PSWF_factors(maxL, n_list, size(a_lms{1},1), R);

psi_lNs = cell(maxL+1,1);
for l = 0:maxL
    psi_lNs{l+1} = cell(2*l+1, 1); % N
    for N = -l:l
        psi_lNs{l+1}{N+l+1} = (-1)^(N*(N<0))*psi_Nn{N+maxL+1}*beta{l+1}{abs(N)+1}(1:s_list(l+1),:).'; % r x s
    end
end

trials = size(Rots, 3);
projs = zeros(3*L-2, 3*L-2,  trials); 
% projs = zeros(L, L,  trials); 
for t = 1:trials
    [~, projs(:,:,t)] = gen_proj_given_rot_PSWFs(Rots(:,:,t), a_lms, psi_lNs, L);
end
if ~isinf(SNR)
    [projs, ~, ~, sigma] = cryo_addnoise(projs, SNR, 'gaussian');
else
    sigma = 0;
end

% beta_PSWF = 1; T = 1e-1; realFlag = 0;
% [PSWF_Nn_p, ~, pts_in_disc] = precomp_pswf_t(L-1, beta_PSWF, T, realFlag);
% ang_freq = PSWF_Nn_p.ang_freq;
% Psi = 2*[real(PSWF_Nn_p.samples), -imag(PSWF_Nn_p.samples(:, ang_freq>0))];
% Psi(:,ang_freq == 0) = Psi(:,ang_freq == 0)/2;
% 
% [M1, M2, M3] = moments_from_micrograph_steerable_LS(projs, L, Psi, ang_freq, pts_in_disc);
[M1, M2, M3] = moments_from_micrograph_steerable(projs, L, PSWF_quad_int, points_inside_the_circle);
gamma = (L/(3*L-2))^2;

% W = precomp_wigner_weights(maxL);
% [B, L1_list, L2_list, L3_list, blk_id] = precomp_bispectrum_coeffs_LS(maxL, s_list, L, W, psi_lNs, Psi, ang_freq, pts_in_disc);
% M3_analytic = bispectrum_from_harmonics_par_noKloop(a_lms, L, B, W, L1_list, L2_list, L3_list, s_list, blk_id);
