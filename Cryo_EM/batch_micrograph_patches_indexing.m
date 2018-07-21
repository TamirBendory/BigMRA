function micro_batch = batch_micrograph_patches_indexing(I, L, batch_size)
% Splits micrograph into batches of patches, one batch per column

micro_batch = cell(ceil(size(I,1)*size(I,2)/batch_size), 1);
pixels_per_img = size(I,1)*size(I,2);
num_micros = size(I,3);

I = padarray(I, [L-1, L-1]);
idx_micro = cart2inds(size(I), L:size(I,1)-L+1, L:size(I,2)-L+1); % linear indices of original micrograph
sz_img = [size(I,1), size(I,2)];
I = reshape(I, [], num_micros);

batch_num = ceil(pixels_per_img/batch_size);
for batch = 1:batch_num
    idx_vals = idx_micro(batch_size*(batch-1) + 1 : min(batch_size*batch, pixels_per_img));
    batch_size_actual = length(idx_vals);
    
    row_vals = rem(idx_vals-1, sz_img(1)) + 1;
    col_vals = (idx_vals - row_vals)/sz_img(1) + 1;
    
    row_vals = bsxfun(@plus, row_vals, -(L-1):(L-1));
    col_vals = bsxfun(@plus, col_vals, -(L-1):(L-1));
    
    patch_inds = bsxfun(@plus, row_vals, permute((col_vals - 1).*sz_img(1), [1,3,2]));
    patch_inds = reshape(patch_inds, batch_size_actual, (2*L-1)^2);
    
    micro_batch{batch} = reshape(I(patch_inds(:), :), batch_size_actual, 2*L-1, 2*L-1, num_micros);
end