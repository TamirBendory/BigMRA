function [r_inds, c_inds, blk_id, sz, nz] = get_blkdiag_inds(ang_freq)

ang_freq_u = unique(ang_freq);
blk_sizes = zeros(size(ang_freq_u));
for ii = 1:length(ang_freq_u)
    blk_sizes(ii) = length(find(ang_freq == ang_freq_u(ii)));
end

sz = sum(blk_sizes)*[1, 1]; % size of sparse blkdiag matrix
nz = sum(blk_sizes.^2);     % # of nonzero elements

blk_sizes_cumsum = [0 ; cumsum(blk_sizes)];
r_inds = []; c_inds = []; blk_id = [];
for ii = 1:length(ang_freq_u)
    curr_inds = (blk_sizes_cumsum(ii) + 1):blk_sizes_cumsum(ii+1);
    [rows, cols] = ndgrid(curr_inds, curr_inds);
    r_inds = [r_inds; rows(:)]; c_inds = [c_inds; cols(:)];
    blk_id = [blk_id; ii*ones(blk_sizes(ii)^2, 1)];
end

