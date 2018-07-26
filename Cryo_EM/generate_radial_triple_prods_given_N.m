function T = generate_radial_triple_prods_given_N(Lr, N1, N2, N3, n_list, r, phi, PSWF_Nn_p, PSWF_quad_int)

L = 2*Lr+1;
c = Lr*pi;
n1_len = n_list(abs(N1)+1);
n2_len = n_list(abs(N2)+1);
n3_len = n_list(abs(N3)+1);

psi_N1 = PSWF_2D_N(N1, n1_len-1, c, r(:), phi(:));
psi_N1 = reshape(psi_N1, size(r,1), size(r,2), []);
psi_N2 = PSWF_2D_N(N2, n2_len-1, c, r(:), phi(:));
psi_N2 = reshape(psi_N2, size(r,1), size(r,2), []);
psi_N3 = PSWF_2D_N(N3, n3_len-1, c, r(:), phi(:));
psi_N3 = reshape(psi_N3, size(r,1), size(r,2), []);

% Create a dr1 x dr2 image:
sz_img = size(r);
maxK = max(PSWF_Nn_p.ang_freq);
q_list = zeros(maxK+1);
for k = 0:maxK, q_list(k+1) = length(find(PSWF_Nn_p.ang_freq == k)); end
blk_id = [];
for k = 0:maxK
    blk_id = [blk_id; k*ones(q_list(k+1)^2*n1_len*n2_len*n3_len, 1)];
end

acc_vec = 0;
parfor ii = 1:prod(sz_img)
    [row, col] = ind2sub(sz_img, ii);
    img = zeros(2*L+1, 2*L+1, n2_len+n3_len);
    inds_row = max(row-L, 1):min(row+L, sz_img(1));
    inds_col = max(col-L, 1):min(col+L, sz_img(2));
    img(inds_row - row + L+1, inds_col - col + L+1, 1:n2_len)...
        = psi_N2(inds_row, inds_col, :);
    img(inds_row - row + L+1, inds_col - col + L+1, n2_len+1:end)...
        = psi_N3(inds_row, inds_col, :);
    
    T_N = pswf_t_f_fast(img, L, 1, 1e-1, 1, PSWF_Nn_p, PSWF_quad_int, '../nfft');
    T_N = PSWF_coeff_convert(T_N, PSWF_Nn_p.ang_freq);
    T_n2 = cellfun(@(x) x(:,1:size(psi_N2,2)), T_N, 'UniformOutput', 0);
    T_n3 = cellfun(@(x) x(:,1+size(psi_N2,2):end), T_N, 'UniformOutput', 0);
    T_N = cellfun(@(x,y) x(:)*y(:).', T_n2, T_n3, 'UniformOutput', 0); % (q1, n2) x (q2, n3)
    T_N = cellfun(@(x) squeeze(conj(psi_N1(row, col, :)))*x(:)', T_N, 'UniformOutput', 0); % n1 x (q1, n2, q2, n3)
    T_N = cellfun(@(x)x(:), T_N, 'UniformOutput', 0);
    acc_vec = acc_vec + vertcat(T_N{:});
end

T = accumarray(blk_id, acc_vec, [maxK+1, 1], @(x) {x});
for k = 0:maxK
    T{k+1} = reshape(T{k+1}, n1_len, q_list(k+1), n2_len, q_list(k+1), n3_len);
    T{k+1} = permute(T{k+1}, [1,3,5,2,4]);
end
    