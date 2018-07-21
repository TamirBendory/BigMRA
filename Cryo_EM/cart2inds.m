function inds = cart2inds(sz, row_vals, col_vals)

row_vals = row_vals(row_vals > 0);
col_vals = col_vals(col_vals > 0);
inds = zeros(length(row_vals)*length(col_vals), 1);
for ii = 1:length(col_vals)
    inds( (ii-1)*length(row_vals) + (1:length(row_vals)) ) = (col_vals(1)-1)*sz(2) + (ii-1)*sz(1) + row_vals;
end