function img = extract_patch(I, row, col, L)

img = zeros(2*L+1, 2*L+1, size(I,3), 'double');
    inds_row = max(row-L, 1):min(row+L, size(I,1));
    inds_col = max(col-L, 1):min(col+L, size(I,2));
    img(inds_row - row + L+1, inds_col - col + L+1, :)...
        = I(inds_row, inds_col, :);
