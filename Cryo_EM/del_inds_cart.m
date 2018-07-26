function inds = del_inds_cart(inds, sz, row_vals, col_vals)

row_vals = row_vals(row_vals > 0);
col_vals = col_vals(col_vals > 0);
for ii = 1:length(col_vals)
    idx = find(inds == (col_vals(1)-1)*sz(2) + row_vals(1), 1);
    inds(idx:idx + length(row_vals)-1) = [];
%     inds( (col_vals(1)-1)*sz(2) + (ii-1)*sz(1) + row_vals ) = inf;
end