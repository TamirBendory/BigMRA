function [m1, m2, m3, PSWF_quad_int, points_inside_the_circle, plan] = moments_from_micrograph_steerable(I, L, PSWF_quad_int, points_inside_the_circle, plan)
% L is size of single projection
m1 = mean(I(:));

L = L-1;        % Conform to PSWF convention where images are 2*L+1
beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values stand for greater oversampling
c = beta*pi*L;  % Bandlimit 
T = 1e-1;       % Truncation parameter
realFlag = 1;   % Flag indicating whether the data is assumed to be real-valued

if ~exist('PSWF_quad_int', 'var') || isempty(PSWF_quad_int) || ~exist('points_inside_the_circle','var') || isempty(points_inside_the_circle)
    [~, PSWF_quad_int, points_inside_the_circle, plan] = precomp_pswf_t(L, beta, T, realFlag);
end

ang_freq = PSWF_quad_int.ang_freq;
ang_freq_u = unique(ang_freq);

m2 = zeros(length(find(ang_freq == 0)), 1);

blk_id = []; q_list = [];
for ii = 1:length(ang_freq_u)
    num_freqs = sum(ang_freq == ang_freq_u(ii));
    q_list(end+1) = num_freqs;
    blk_id(end+1:end+num_freqs^2, 1) = ii*ones(num_freqs^2, 1);
end
m3 = zeros(size(blk_id));

sz_img = size(I(:,:,1));
for col = 1:sz_img(2)
    for row = 1:sz_img(1)
%       form image to be expanded:
        img = extract_patch(I, row, col, L);
    
%       Expand in PSWFs:
        coeff = pswf_t_short(img, L, beta, realFlag, PSWF_quad_int, points_inside_the_circle, plan);
        coeff = mat2cell(coeff, q_list, 1);
    
%       Compute contribution to moments:
        m2 = m2 + coeff{1}*reshape(img(L+1,L+1,:), [], 1);
        m3_add = cellfun(@(x) x*diag( reshape(img(L+1,L+1,:), [], 1) )*x', coeff, 'UniformOutput', 0);
        m3_add = cellfun(@(x) real(x(:)), m3_add, 'UniformOutput', 0);
        m3 = m3 + vertcat(m3_add{:});
    end
end

m2 = m2./numel(I);
m3 = m3./numel(I);

m3 = accumarray(blk_id, m3, [length(ang_freq_u), 1], @(x) {reshape(x, sqrt(length(x)), [])});
