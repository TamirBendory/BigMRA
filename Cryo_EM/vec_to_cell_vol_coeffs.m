function a_lms = vec_to_cell_vol_coeffs(a_lms, L_list)

a_lms = accumarray(L_list+1, a_lms, [max(L_list)+1, 1], @(x) {x});
for l = 0:length(a_lms)-1
    a_lms{l+1} = reshape(a_lms{l+1}, [], l+1);
    a_lms{l+1}(:,1) = (a_lms{l+1}(:,1) + (-1)^l*conj(a_lms{l+1}(:,1)))/2;
    neg_m = bsxfun(@times, conj(a_lms{l+1}(:,end:-1:2)), (-1).^(l+(l:-1:1)));
    a_lms{l+1} = [neg_m, a_lms{l+1}];
end