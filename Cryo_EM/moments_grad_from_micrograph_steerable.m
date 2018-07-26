function [m1, m2, m3, m3_g] = moments_grad_from_micrograph_steerable(I, L, L_list, PSWF_Nn_p, PSWF_quad_int)
% L is size of single projection
m1 = mean(I(:));

L = L-1;        % Conform to PSWF convention where images are 2*L+1
beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values stand for greater oversampling
T = 1e-1;       % Truncation parameter
realFlag = 1;   % Flag indicating whether the data is assumed to be real-valued

if ~exist('PSWF_quad_int', 'var') || isempty(PSWF_quad_int) || ~exist('points_inside_the_circle', 'var') || isempty(points_inside_the_circle)
    [PSWF_Nn_p, PSWF_quad_int, ~] = precomp_pswf_t(L, beta, T, realFlag);
end

m2 = zeros(length(find(ang_freq == 0)), 1);

blk_id = [];
for ii = 1:length(ang_freq_u)
    num_freqs = length(find(ang_freq == ang_freq_u(ii)));
    blk_id = [blk_id; ii*ones(num_freqs^2, 1)];
end
m3 = zeros(size(blk_id));

sz_img = size(I(:,:,1));
for ii = 1:numel(I(:,:,1))
    [row, col] = ind2sub(sz_img, ii);
%     form image to be expanded:
    img = extract_patch(I, row, col, L);
    
%     Expand in PSWFs:
    coeff = pswf_t_f_fast(img, L, beta, T, realFlag, PSWF_Nn_p, PSWF_quad_int, '../nfft');
    
%     Compute contribution to moments:
    m2 = m2 + coeff{1}*squeeze(img(L+1,L+1,:));
    m3_add = cellfun(@(x) real(x*diag(squeeze(img(L+1,L+1,:)))*x'), coeff, 'UniformOutput', 0);
    m3_add = cellfun(@(x) x(:), m3_add, 'UniformOutput', 0);
    m3 = m3 + vertcat(m3_add{:});
    
    img_lms = extract_patch(I_lms, row, col, L);
    coeff_lms = pswf_t_f_fast(img_lms, L, beta, T, realFlag, PSWF_Nn_p, PSWF_quad_int, '../nfft');
    G_add = cellfun(@(x) zeros(size(x,1), size(x,2), length(L_list)), m3_add, 'UniformOutput', 0);
    for k = 1:length(G_add)
        for lms = 1:length(L_list)
            G_add{k} = real(coeff{k}*diag(squeezeimg_lms(L+1,L+1,:,lms))*coeff{k}'); % q1 x q2
        G_add{k}
end

m2 = m2./numel(I);
m3 = m3./numel(I);

m3 = accumarray(blk_id, m3, [length(ang_freq_u), 1], @(x) {reshape(x, sqrt(length(x)), [])});