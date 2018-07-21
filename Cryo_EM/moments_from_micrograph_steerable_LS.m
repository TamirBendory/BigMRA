function [m1, m2, m3] = moments_from_micrograph_steerable_LS(I, L, Psi, ang_freq, pts_in_disc)
% L is size of single projection

m1 = mean(I(:));

maxN = max(ang_freq);
blk_id = [];
for ii = 0:maxN
    num_freqs = sum(ang_freq == ii);
    blk_id(end+1:end+num_freqs^2, 1) = (ii+1)*ones(num_freqs^2, 1);
end

m2 = zeros(length(find(ang_freq == 0)), 1);
m3 = zeros(size(blk_id));

sz_img = size(I(:,:,1));
parfor ii = 1:numel(I(:,:,1))
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

m3 = accumarray(blk_id, m3, [maxN+1, 1], @(x) {reshape(x, sqrt(length(x)), [])});