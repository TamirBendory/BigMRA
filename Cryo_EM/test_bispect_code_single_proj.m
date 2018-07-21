%% Paths

clear all; close all;
addpath('../kam_cryo')
addpath('../SHT')
addpath(genpath('../pswf_t_mod'))
addpath(genpath('../easyspin-5.2.18'))
addpath('~/EMDB_maps')
addpath('../aspire')
initpath

%% Basic parameters and volume
info.mol = 'emd_8117'; % molecule id
info.maxL = 1; % truncation for spherical harmonics expansion
info.r_cut = 1/2;
info.L0 = 20; % size of volume
info.N=floor(info.L0/2);

vol = double(ReadMRC([info.mol '.map'])); %read volume
vol = cryo_downsample(vol, info.L0); 

% Expand volume:
[Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis(info);
[a_lms, vol_rec] = expand_vol_spherical_basis(vol, info, Psilms, jball);

%% Generate basis & projection
L = info.L0;
R = floor(L/2);
maxL = length(a_lms)-1;

[psi_Nn, n_list] = PSWF_2D_full_cart(maxL, R, 1, 1e-5);
s_list = cellfun('size', a_lms, 1);
beta = sph_Bessel_to_2D_PSWF_factors(maxL, n_list, size(a_lms{1},1), R);

psi_lNs = cell(maxL+1,1);
for l = 0:maxL
    psi_lNs{l+1} = cell(2*l+1, 1); % N
    for N = -l:l
        psi_lNs{l+1}{N+l+1} = (-1)^(N*(N<0))*psi_Nn{N+maxL+1}*beta{l+1}{abs(N)+1}(1:s_list(l+1),:).'; % r x s
    end
end

Rot = randR(3);
R = RN2RL( getSHrotMtx( Rot, length(a_lms)-1, 'complex' ) );
a_lms_rot = cellfun(@mtimes, a_lms, R, 'UniformOutput', 0);

I = zeros(L);
for l = 0:maxL
    for N = -l:l
        I(:) = I(:) + psi_lNs{l+1}{N+l+1}*a_lms_rot{l+1}(:,N+l+1);
    end
end

%% Basis for moments:

beta_PSWF = 1; T = 1e-5; realFlag = 0;
[PSWF_Nn_p, ~, pts_in_disc] = precomp_pswf_t(L-1, beta_PSWF, T, realFlag);
ang_freq = PSWF_Nn_p.ang_freq;
Psi = 2*[real(PSWF_Nn_p.samples), -imag(PSWF_Nn_p.samples(:, ang_freq>0))];
Psi(:,ang_freq == 0) = Psi(:,ang_freq == 0)/2;

maxN = max(ang_freq);
blk_id = []; q_list = [];
for ii = 0:maxN
    num_freqs = sum(ang_freq == ii);
    blk_id(end+1:end+num_freqs^2, 1) = (ii+1)*ones(num_freqs^2, 1);
    q_list(ii+1) = num_freqs;
end

%% Compute moments directly:

m2 = zeros(length(find(ang_freq == 0)), 1);
m3 = zeros(size(blk_id));

sz_img = size(I);
parfor ii = 1:numel(I)
    [row, col] = ind2sub(sz_img, ii);
%     form image to be expanded:
    img = extract_patch(I, row, col, L-1);
    assert(norm(squeeze(img(L, L, :) - I(row, col, :)))==0)
    img = reshape(img(1:end-1, 1:end-1, :), 4*(L-1)^2, size(img,3));
    
%     Expand in PSWFs:
    x = Psi\img(pts_in_disc, :);
    coeff = x(1:length(ang_freq),:); 
    coeff(ang_freq>0, :) = coeff(ang_freq>0, :) + 1i*x(length(ang_freq)+1:end, :);
    coeff = PSWF_coeff_convert(coeff, ang_freq);
    
%     Compute contribution to moments:
    m2 = m2 + coeff{1}*squeeze(I(row,col,:));
    m3_add = cellfun(@(x) real(x*diag(squeeze(I(row,col,:)))*x'), coeff, 'UniformOutput', 0);
    m3_add = cellfun(@(x) x(:), m3_add, 'UniformOutput', 0);
    m3 = m3 + vertcat(m3_add{:});
end

m2 = m2./numel(I);
m3 = m3./numel(I);

%% Compute moments from expansion:

% generate vectorization indices:
L_inds = []; N_inds = [];
ii = 1;
for l = 0:maxL
    for N = -l:l
        L_inds(end+1) = l;
        N_inds(end+1) = N;
    end
end

vec = @(x) x(:);
q_cell = vec(num2cell(q_list));

%% Partial separation:

m3_part_sep = zeros(size(blk_id));
parfor ii = 1:length(L_inds)^3*L^2
    [ii_1, ii_2, ii_3, row, col] = ind2sub([length(L_inds), length(L_inds), length(L_inds), L, L], ii);
    L1 = L_inds(ii_1); N1 = N_inds(ii_1);
    L2 = L_inds(ii_2); N2 = N_inds(ii_2);
    L3 = L_inds(ii_3); N3 = N_inds(ii_3);
    
    I1 = reshape(psi_lNs{L1+1}{N1+L1+1}*a_lms_rot{L1+1}(:,N1+L1+1), L, L, []);
    I2 = reshape(psi_lNs{L2+1}{N2+L2+1}*a_lms_rot{L2+1}(:,N2+L2+1), L, L, []);
    I3 = reshape(psi_lNs{L3+1}{N3+L3+1}*a_lms_rot{L3+1}(:,N3+L3+1), L, L, []);
    
    patch = extract_patch(I2, row, col, L-1);
    patch(:,:,end+1) = extract_patch(I3, row, col, L-1);
    patch = reshape(patch(1:end-1,1:end-1,:), 4*(L-1)^2, size(patch,3));
    
    x = Psi\patch(pts_in_disc, :);
    T = x(1:length(ang_freq),:);
    T(ang_freq>0, :) = T(ang_freq>0, :) + 1i*x(length(ang_freq)+1:end, :);
    T = mat2cell(T, q_list, 2);
    
    T = cellfun(@(x) vec(x(:,1))*diag(squeeze(I1(row,col,:)))*vec(x(:,2))', T, 'UniformOutput', 0); % q1 x q2
    
    % Contract tensor:
    for k = 1:length(T)
        T{k} = T{k}(:);
    end
    
    m3_part_sep = m3_part_sep + vertcat(T{:});
end
m3_part_sep = m3_part_sep./L^2;

%% Reshape psi_lNs back:
for l = 0:maxL
    for N = -l:l
        psi_lNs{l+1}{N+l+1} = reshape(psi_lNs{l+1}{N+l+1}, L, L, []);
    end
end

%% Full separation - no Wigner:

m3_exp = zeros(size(blk_id));
parfor ii = 1:length(L_inds)^3*L^2
    [ii_1, ii_2, ii_3, row, col] = ind2sub([length(L_inds), length(L_inds), length(L_inds), L, L], ii);
    L1 = L_inds(ii_1); N1 = N_inds(ii_1);
    L2 = L_inds(ii_2); N2 = N_inds(ii_2);
    L3 = L_inds(ii_3); N3 = N_inds(ii_3);
    
    patch = extract_patch(psi_lNs{L2+1}{N2+L2+1}, row, col, L-1);
    patch(:,:,end+1:end+s_list(L3+1)) = extract_patch(psi_lNs{L3+1}{N3+L3+1}, row, col, L-1);
    patch = reshape(patch(1:end-1,1:end-1,:), 4*(L-1)^2, size(patch,3));
    
    x = Psi\patch(pts_in_disc, :);
    T = x(1:length(ang_freq),:);
    T(ang_freq>0, :) = T(ang_freq>0, :) + 1i*x(length(ang_freq)+1:end, :);
    T = mat2cell(T, q_list, s_list(L2+1) + s_list(L3+1));
    
    T = cellfun(@(x) vec(x(:,1:s_list(L2+1)))*vec(x(:,s_list(L2+1)+1:end))', T, 'UniformOutput', 0); % (q1, s2) x (q2, s3)
    
    T = cellfun(@(x,y) reshape(permute(reshape(x, y, s_list(L2+1), y, s_list(L3+1)), [1,3,4,2]), y^2, []), T, q_cell, 'UniformOutput', 0); %(q1, q2) x (s3, s2)
    
    T = vec(vertcat(T{:}))*reshape(psi_lNs{L1+1}{N1+L1+1}(row, col, :), 1, s_list(L1+1)); % (q1, q2, s3, s2) x s1
    
    T = T*a_lms{L1+1}; % (q1, q2, s3, s2) x m1
    T = reshape(T.', [] ,s_list(L2+1))*a_lms{L2+1}; % (m1, q1, q2 ,s3) x m2
    T = reshape(T.', [], s_list(L3+1))*conj(a_lms{L3+1}); % (m2, m1, q1, q2) x m3
    T = reshape(T.', 2*L3+1, 2*L2+1, 2*L1+1, []); % (m3, m2, m1, q1, q2)
    
    m3_add = 0;
    for m1 = -L1:L1
        for m2 = -L2:L2
            m3 = m1+m2;
            if abs(m3) > L3, continue; end
            m3_add = m3_add + R{L1+1}(m1+L1+1, N1+L1+1)*R{L2+1}(m2+L2+1, N2+L2+1)*conj(R{L3+1}(m3+L3+1, N3+L3+1))...
                *squeeze(T(m3+L3+1, m2+L2+1, m1+L1+1, :));
        end
    end
    
    m3_exp = m3_exp + m3_add;
end
m3_exp = m3_exp./L^2;

%% Check convergence of triple Wigner D product

L1 = randi(maxL+1)-1; L2 = randi(maxL+1)-1;
L3_vals = abs(L1-L2):min(L1+L2, maxL);
L3_ind = randi(length(L3_vals));
L3 = L3_vals(L3_ind);

N1 = randi(2*L1+1)-L1-1; m1 = randi(2*L1+1)-L1-1;
N2_vals = max(-L2, -L3-N1):min(L2, L3-N1); m2_vals = max(-L2, -L3-m1):min(L2, L3-m1);
N2 = N2_vals(randi(length(N2_vals))); m2 = m2_vals(randi(length(m2_vals)));
m3 = m1+m2; N3 = N1+N2;

D_sample = 0;
% Rots = genRotationsGrid( 800, 0, 0);
% trials = size(Rots, 3);
trials =1e9;
Rots = randrot(3, trials);
parfor ii = 1:trials
    Rot = Rots(:,:,ii);
%     R = analytic_WignerD(maxL, eulang(Rot));
    R = get_WignerD_EasySpin(maxL, Rot);
%     Rot = Rots(:,:,ii);
%     R = RN2RL( getSHrotMtx( Rot, maxL, 'complex' ) );
    
    D_sample = D_sample + R{L1+1}(N1+L1+1, m1+L1+1)*R{L2+1}(N2+L2+1, m2+L2+1)*conj(R{L3+1}(N3+L3+1, m3+L3+1));
end
D_sample = D_sample/trials;


W = precomp_wigner_weights(maxL);
D_anal = W{L1+1,L2+1}(N2+L2+1, N1+L1+1, L3_ind)*W{L1+1,L2+1}(m2+L2+1, m1+L1+1, L3_ind);
