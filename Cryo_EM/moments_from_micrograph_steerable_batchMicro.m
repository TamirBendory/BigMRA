function [m1, m2, m3] = moments_from_micrograph_steerable_batchMicro(I, L, batch_size)
% L is size of single projection
m1 = mean(I(:))/L;

micro_batch = batch_micrograph_patches(I, L, batch_size);
num_pixels = numel(I);
clear I
disp('Done batching')

if isempty(gcp('nocreate')), parpool('local', maxNumCompThreads); end

beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values stand for greater oversampling
T = 1e-1;       % Truncation parameter
realFlag = 1;   % Flag indicating whether the data is assumed to be real-valued

[~, PSWF_quad_int, points_inside_the_circle] = precomp_pswf_t(L-1, beta, T, realFlag);
plan_par = WorkerObjWrapper(@precomp_nfft_par, {PSWF_quad_int, L-1, beta}, @nfft_finalize);
disp('Done precomputing')

ang_freq = PSWF_quad_int.ang_freq;
q_list = PSWF_quad_int.n_list;

m2 = zeros(length(find(ang_freq == 0)), 1);

blk_id = [];
for ii = 0:max(ang_freq)
    blk_id(end+1:end+q_list(ii+1)^2, 1) = ii+1;
end
m3 = zeros(size(blk_id));

sz_img = [size(I,1), size(I,2)];
pixels_per_img = prod(sz_img);
num_micros = size(I,3);
batch_num = ceil(pixels_per_img/batch_size);
for batch = 1:batch_num
    idx_vals = (batch_size*(batch-1) + 1) : min(batch_size*batch, pixels_per_img);
    row_vals = rem(idx_vals-1, sz_img(1)) + 1;
    col_vals = (idx_vals - row_vals)/sz_img(1) + 1;
    batch_size_actual = length(idx_vals);
    batch_curr = zeros(2*L-1, 2*L-1, batch_size_actual*num_micros, 'double');
    for ii = 1:batch_size_actual
        row = row_vals(ii); col = col_vals(ii);
        inds_row = max(row-L+1, 1):min(row+L-1, sz_img(1));
        inds_col = max(col-L+1, 1):min(col+L-1, sz_img(2));
        batch_curr(inds_row - row + L+1, inds_col - col + L+1, (ii-1)*num_micros+1:ii*num_micros)...
            = I(inds_row, inds_col, :);
    end
    micro_batch{batch} = batch_curr;
end

% spmd, mpiprofile('on'); end
parfor t = 1:batch_num
    
    % Generate batch:
    idx_vals = (batch_size*(batch-1) + 1) : min(batch_size*batch, pixels_per_img);
    row_vals = rem(idx_vals-1, sz_img(1)) + 1;
    col_vals = (idx_vals - row_vals)/sz_img(1) + 1;
    batch_size_actual = length(idx_vals);
    batch_curr = zeros(2*L-1, 2*L-1, batch_size_actual*num_micros, 'double');
    for ii = 1:batch_size_actual
        row = row_vals(ii); col = col_vals(ii);
        inds_row = max(row-L+1, 1):min(row+L-1, sz_img(1));
        inds_col = max(col-L+1, 1):min(col+L-1, sz_img(2));
        batch_curr(inds_row - row + L+1, inds_col - col + L+1, (ii-1)*num_micros+1:ii*num_micros)...
            = I(inds_row, inds_col, :);
    end
    
    plan = plan_par.Value;
    
%   Expand in PSWFs:
    coeff = pswf_t_short(micro_batch{t}, L-1, beta, realFlag, PSWF_quad_int, points_inside_the_circle, plan);
    coeff = mat2cell(coeff, q_list, size(coeff,2));
    
    center_pts = reshape(micro_batch{t}(L, L, :), [], 1);
%   Compute contribution to moments:
    m2 = m2 + coeff{1}*center_pts;
    m3_add = cellfun(@(x) real(x*diag(center_pts)*x'), coeff, 'UniformOutput', 0);
    m3_add = cellfun(@(x) x(:), m3_add, 'UniformOutput', 0);
    m3 = m3 + vertcat(m3_add{:});
end
% spmd, pInfo = mpiprofile('info'); end
% pInfo = [pInfo{:}];

m2 = m2./num_pixels;
m3 = m3./num_pixels;

m3 = accumarray(blk_id, m3, [max(ang_freq)+1, 1], @(x) {reshape(x, sqrt(length(x)), [])});
