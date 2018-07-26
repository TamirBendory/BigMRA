%% Paths

clear all; close all;
addpath('~/Documents/kam_cryo')
addpath('../SHT')
% addpath(genpath('../Prol'))
addpath(genpath('../pswf_t_mod'))
addpath(genpath('../pswf_t_fast_nfft_v3.3.2'))
addpath('~/EMDB_maps')
addpath('../aspire')
initpath

%% Basic parameters and volume
info.mol = 'emd_8117'; % molecule id
info.maxL = 5; % truncation for spherical harmonics expansion
info.r_cut = 1/2;
info.L0 = 31; % size of volume
info.N=floor(info.L0/2);

vol = double(ReadMRC([info.mol '.map'])); %read volume
vol = cryo_downsample(vol, info.L0); 

% Expand volume:
[Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(info);
[a_lms, vol_rec] = expand_vol_spherical_basis(vol, info, Psilms, jball);

%% Test 3D Fourier -> 2D PSWFs transform

% Generate random projection:
[I, D_l] = gen_rand_proj(a_lms, Psilms_2D, jball_2D, info.L0);

L = info.L0;
beta_PSWF = 1;
T = 1e-1;
[psi_Nn, n_list] = PSWF_2D_full_cart(length(a_lms)-1, floor(L/2), beta_PSWF, T);
% psi_Nn = PSWF_2D_full_cart_OLD(length(a_lms)-1, n_list, floor(L/2));
s_list = cellfun('size', a_lms, 1);
beta = sph_Bessel_to_2D_PSWF_factors(info.maxL, n_list, size(a_lms{1},1), info.N);

% Compare analytic and numeric expansion coefficients in 2D prolates:
[coeffs_numer2, PSWF_Nn_p, PSWF_quad_int, points_inside_the_circle] = pswf_t_f_fast(I, info.N, beta_PSWF, T, 1, [], [], 4, '../nfft');

[Mt, ang_freq] = precomp_pswf_t_mat(info.N, beta_PSWF, T);
coeffs_numer = Mt*I(:);

rec_img_numer = pswf_t_b( coeffs_numer, PSWF_Nn_p, 1 );
coeffs_numer = PSWF_coeff_convert(coeffs_numer, PSWF_Nn_p.ang_freq);

coeffs_analytic = cell(2*info.maxL+1,1);
a_lms_rot = cellfun(@mtimes, a_lms, D_l, 'UniformOutput', 0);
for N = -info.maxL:info.maxL
    coeffs_analytic{N+info.maxL+1} = zeros(n_list(abs(N)+1),1);
    for l = abs(N):info.maxL
        coeffs_analytic{N+info.maxL+1} = coeffs_analytic{N+info.maxL+1} + ...
            (-1)^(N*(N<0))*beta{l+1}{abs(N)+1}(1:size(a_lms_rot{l+1},1),:).'*a_lms_rot{l+1}(:,N+l+1);
    end
end

rec_img_analytic = zeros(L);
for N = -info.maxL:info.maxL
    rec_img_analytic(:) = rec_img_analytic(:) + psi_Nn{N+info.maxL+1}*coeffs_analytic{N+info.maxL+1};
end

psi_lNs = cell(size(psi_Nn));
for l = 0:info.maxL
    psi_lNs{l+1} = cell(2*l+1, 1); % N
    for N = -l:l
        psi_lNs{l+1}{N+l+1} = (-1)^(N*(N<0))*psi_Nn{N+info.maxL+1}*beta{l+1}{abs(N)+1}(1:s_list(l+1),:).'; % r x s
    end
end

rec_img_analytic_2 = zeros(L);
for l = 0:info.maxL
    for N = -l:l
        rec_img_analytic_2(:) = rec_img_analytic_2(:) + psi_lNs{l+1}{N+l+1}*a_lms_rot{l+1}(:,N+l+1);
    end
end

norm(rec_img_analytic_2(1:end-1,1:end-1) - rec_img_numer,'fro')/norm(rec_img_numer(:))

%% Test expansion about a different center

maxN = info.maxL;
b_Nn = coeffs_analytic;

[PSWF_Nn_p, PSWF_quad_int, pts_in_disc] = precomp_pswf_t(L-1, beta_PSWF, T, 0);
ang_freq = PSWF_Nn_p.ang_freq;
Psi = 2*[real(PSWF_Nn_p.samples), -imag(PSWF_Nn_p.samples(:, ang_freq>0))];
Psi(:,ang_freq == 0) = Psi(:,ang_freq == 0)/2;

sz_img = size(I);
row = randi(size(I, 1)); col = randi(size(I,2));
L = info.L0;
img = extract_patch(rec_img_analytic, row, col, L-1);
img = img(1:end-1, 1:end-1);
x = Psi\img(pts_in_disc);
a_numer = x(1:length(ang_freq)); 
a_numer(ang_freq>0) = a_numer(ang_freq>0) + 1i*x(length(ang_freq)+1:end);

a_patchExp = 0;
for N = -info.maxL:info.maxL
    patch = extract_patch(reshape(psi_Nn{N+maxN+1},L,L,[]), row, col, L-1);
    patch = reshape(patch(1:end-1,1:end-1,:), 4*(L-1)^2, []);
    x = Psi\patch(pts_in_disc,:);
    a_patch = x(1:length(ang_freq),:); 
    a_patch(ang_freq>0,:) = a_patch(ang_freq>0,:) + 1i*x(length(ang_freq)+1:end,:);
    a_patchExp = a_patchExp + a_patch*coeffs_analytic{N+maxN+1};
end

norm(a_numer - a_patchExp)/norm(a_numer)
% a_numer = pswf_t_short(img, L-1, 1, 0, PSWF_quad_int, points_inside_the_circle);
img_rec_numer = pswf_t_b( a_numer, PSWF_Nn_p, 1 );
keyboard

% a_numer_rec = pswf_t_f_fast(img, L-1, 1, 1e-1, 0, PSWF_Nn_p, PSWF_quad_int, '../nfft');
a_numer_rec = pswf_t_short(img_rec_numer, L-1, 1, 0, PSWF_quad_int, points_inside_the_circle);
norm(a_numer-a_numer_rec)

a_numer = PSWF_coeff_convert(a_numer, PSWF_Nn_p.ang_freq);

q_list = cellfun('size', a_numer, 1);
a_analytic = cell(length(q_list),1);
a_analytic = cellfun(@(x) 0, a_analytic, 'UniformOutput', 0);
for N = -maxN:maxN
    patch = extract_patch(psi_Nn{N+maxN+1}, row, col, L-1);
    T_N  = pswf_t_f_fast(patch, L-1, 1, 1e-1, 0, PSWF_Nn_p, PSWF_quad_int, '../nfft');
    T_N = PSWF_coeff_convert(T_N, PSWF_Nn_p.ang_freq);
    
    a_analytic = cellfun(@(x,y) x + y*b_Nn{N+maxN+1}, a_analytic, T_N, 'UniformOutput', 0);
end

img_rec_analytic = zeros(2*L-2);
psi_Kq = PSWF_2D_full_cart(length(a_analytic)-1, L-1);
psi_Kq = psi_Kq(length(a_analytic):end);
for N = 0:length(a_analytic)-1
    if N == 0
        img_rec_analytic(:) = img_rec_analytic(:) + real(psi_Kq{N+1}*a_analytic{N+1});
    else
        img_rec_analytic(:) = img_rec_analytic(:) + 2*real(psi_Kq{N+1}*a_analytic{N+1});
    end
end


%% Test moment expansion:

% Compute numeric steerable moments:
[m1, m2, m3] = moments_from_micrograph_steerable(rec_img_analytic, L);

% Compare with anayltic:
c = pi*info.N;
n_list = cellfun('size', coeffs_analytic, 1);
blk_id = []; q_list = {};
for ii = 1:max(PSWF_Nn_p.ang_freq)+1
    num_freqs = length(find(PSWF_Nn_p.ang_freq == ii-1));
    q_list = [q_list; num_freqs];
end

m3_anal = cellfun(@(x) 0, q_list, 'UniformOutput', 0);
for ii = 1:(2*maxN+1)^3
    [N1, N2, N3] = ind2sub((2*maxN+1)*[1,1,1], ii);
    N1 = N1-maxN-1; N2 = N2-maxN-1; N3 = N3-maxN-1;
    
    I1 = psi_Nn{N1+maxN+1}; n1_len = size(I1,3);
    I2 = psi_Nn{N2+maxN+1}; n2_len = size(I2,3);
    I3 = psi_Nn{N3+maxN+1}; n3_len = size(I3,3);
    blk_id = [];
    for ii = 1:max(PSWF_Nn_p.ang_freq)+1
        num_freqs = length(find(PSWF_Nn_p.ang_freq == ii-1));
        blk_id = [blk_id; ii*ones(num_freqs^2*n1_len*n2_len*n3_len, 1)];
    end
    
    T = gen_triple_prod(I1, I2, I3, PSWF_Nn_p, PSWF_quad_int, blk_id, q_list);
    T = cellfun(@(x) x(:,:).'*b_Nn{N1+maxN+1}, T, 'UniformOutput', 0); % (q1, q2, n3, n2)
    T = cellfun(@(x) reshape(x, [], n2_len)*b_Nn{N2+maxN+1}, T, 'UniformOutput', 0); % (q1, q2, n3)
    T = cellfun(@(x) reshape(x, [], n3_len)*conj(b_Nn{N3+maxN+1}), T, 'UniformOutput', 0); % (q1, q2)
    m3_anal = cellfun(@plus, m3_anal, T, 'UniformOutput', 0);
end

