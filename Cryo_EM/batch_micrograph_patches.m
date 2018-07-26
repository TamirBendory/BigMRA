function micro_batch = batch_micrograph_patches(I, L, batch_size)
% Splits micrograph into batches of patches, one batch per column

micro_batch = cell(ceil(size(I,1)*size(I,2)/batch_size), 1);
pixels_per_img = size(I,1)*size(I,2);
num_micros = size(I,3);

I = padarray(I, [L-1, L-1]);
idx_micro = cart2inds(size(I), L:size(I,1)-L+1, L:size(I,2)-L+1); % linear indices of original micrograph
sz_img = [size(I,1), size(I,2)];

batch_num = ceil(pixels_per_img/batch_size);
for batch = 1:batch_num
    idx_vals = idx_micro(batch_size*(batch-1) + 1 : min(batch_size*batch, pixels_per_img));
    row_vals = rem(idx_vals-1, sz_img(1)) + 1;
    col_vals = (idx_vals - row_vals)/sz_img(1) + 1;
    batch_size_actual = length(idx_vals);
    batch_curr = zeros(2*L-1, 2*L-1, batch_size_actual*num_micros, 'double');
    for ii = 1:batch_size_actual
        row = row_vals(ii); col = col_vals(ii);
        inds_row = row + (-(L-1):(L-1));
        inds_col = col + (-(L-1):(L-1));
        batch_curr(:, :, (ii-1)*num_micros+1:ii*num_micros)...
            = I(inds_row, inds_col, :);
    end
    micro_batch{batch} = batch_curr;
end