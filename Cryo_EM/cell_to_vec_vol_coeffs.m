function a_lms = cell_to_vec_vol_coeffs(a_lms)

a_lms = cellfun(@(x) x(:, (size(x,2)-1)/2+1:end), a_lms, 'UniformOutput', 0);
a_lms = cellfun(@(x) x(:), a_lms, 'UniformOutput', 0);
a_lms = vertcat(a_lms{:});