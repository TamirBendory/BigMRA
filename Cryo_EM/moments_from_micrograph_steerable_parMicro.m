function [m1, m2, m3] = moments_from_micrograph_steerable_parMicro(I, L, batch_size)
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
cntr_idx = sub2ind([2*L-1, 2*L-1], L, L); % linear index of center of patch

pixels_per_img = size(I,1)*size(I,2); % pixels per micrograph
num_micros = size(I,3);
batch_num = ceil(pixels_per_img/batch_size); % number of batches per micrograph

I = padarray(I, [L-1, L-1]);
idx_micro = cart2inds(size(I), L:size(I,1)-L+1, L:size(I,2)-L+1); % linear indices of original micrograph
sz_img = size(I(:,:,1));

if isempty(gcp('nocreate')), parpool('local', maxNumCompThreads); end
parfor t = 1:num_micros
    I_curr = I(:,:,t);
    
    for batch = 1:batch_num
        idx_vals = idx_micro(batch_size*(batch-1) + 1 : min(batch_size*batch, pixels_per_img)); % linear indices of current batch
        row_vals = rem(idx_vals-1, sz_img(1)) + 1; % row indices
        col_vals = (idx_vals - row_vals)/sz_img(1) + 1; % col indices
        batch_size_actual = length(idx_vals); % size of current batch
        
        batch_curr = zeros(2*L-1, 2*L-1, batch_size_actual, 'double');
        for ii = 1:batch_size_actual
            row = row_vals(ii); col = col_vals(ii);
            batch_curr(:, :, ii)...
                = I_curr(row-L+1:row+L-1, col-L+1:col+L-1);
        end
        batch_curr = reshape(batch_curr, (2*L-1)^2, batch_size_actual);
        batch_curr(pts_notin_disc, :) = 0;
    
%       Expand in PSWFs:
        coeff = Mt*batch_curr;
        m3_add = zeros(length(blk_id), 1, 'double');
        for N = 0:maxN
            if N == 0
                m2 = m2 + real(coeff(1:q_list(1), :))*reshape(batch_curr(cntr_idx,:), batch_size_actual, 1);
            end
            
            coeff_N = coeff(q_cumsum(N+1)+1: q_cumsum(N+2), :);
            tmp = bsxfun(@times, coeff_N, batch_curr(cntr_idx, :));
            tmp = real(tmp*coeff_N');
            m3_add(q_sq_cumsum(N+1)+1: q_sq_cumsum(N+2)) = tmp(:);
        end
        m3 = m3 + m3_add;
    end
    
end

m2 = m2./(pixels_per_img*num_micros);
m3 = m3./(pixels_per_img*num_micros);

m3 = accumarray(blk_id, m3, [max(ang_freq)+1, 1], @(x) {reshape(x, sqrt(length(x)), [])});
