function [m1, m2, m3] = moments_from_micrograph_steerable_naive(I, L)
% L is size of single projection
m1 = mean(I(:))/L;

beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values stand for greater oversampling
T = 1e-1;       % Truncation parameter

[Mt, ang_freq] = precomp_pswf_t_mat(L-1, beta, T);
maxN = max(ang_freq);
q_list = zeros(maxN+1, 1);
for ii = 0:maxN, q_list(ii+1) = sum(ang_freq == ii); end
q_cumsum = cumsum([0; q_list(:)]);
q_sq_cumsum = cumsum([0; q_list(:).^2]);

m2 = zeros(sum(ang_freq == 0), 1);

blk_id = [];
for ii = 0:maxN
    blk_id(end+1:end+q_list(ii+1)^2, 1) = ii+1;
end
m3 = zeros(size(blk_id));

[x,y] = meshgrid(-L+1:L-1, -L+1:L-1); pts_notin_disc = sqrt(x.^2 + y.^2) > L-1;

I = padarray(I, [L-1, L-1]);
sz_img = size(I(:,:,1));

for t = 1:size(I,3)
    for row = L:sz_img(1)-L+1
        for col = L:sz_img(2)-L+1
        patch_curr = I(row-L+1:row+L-1, col-L+1:col+L-1, t);
        patch_curr(pts_notin_disc) = 0;
        
%       Expand in PSWFs:
        coeff = Mt*patch_curr(:);
        coeff = mat2cell(coeff, q_list, 1);
    
%       Compute contribution to moments:
        m2 = m2 + coeff{1}*patch_curr(L,L);
        m3_add = cellfun(@(x) patch_curr(L,L)*real(x*x'), coeff, 'UniformOutput', 0);
        m3_add = cellfun(@(x) x(:), m3_add, 'UniformOutput', 0);
        m3 = m3 + vertcat(m3_add{:});
        end
    end
end

m2 = m2./numel(I);
m3 = m3./numel(I);

m3 = accumarray(blk_id, m3, [max(ang_freq)+1, 1], @(x) {reshape(x, sqrt(length(x)), [])});
