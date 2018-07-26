function T = gen_triple_prod(I1, I2, I3, PSWF_Kq_p, PSWF_quad_int_Kq, blk_id, q_list)

L = size(I1,1);
T = 0;
n1_len = size(I1,3);
n2_len = size(I2,3);
n3_len = size(I3,3);
beta_PSWF = 1; Trunc = 1e-1; realFlag = 0;
vec = @(x)x(:);
parfor ii = 1:L^2
    [row, col] = ind2sub([L, L], ii);
    patch = extract_patch(I2, row, col, L-1);
    patch(:,:,end+1:end+n3_len) = extract_patch(I3, row, col, L-1);
    
    coeffs  = pswf_t_f_fast(patch, L-1, beta_PSWF, Trunc, realFlag, PSWF_Kq_p, PSWF_quad_int_Kq, '../nfft');
    coeffs = PSWF_coeff_convert(coeffs, PSWF_Kq_p.ang_freq);
    
    coeffs = cellfun(@(x) vec(x(:,1:n2_len))*vec(x(:,n2_len+1:end))', coeffs, 'UniformOutput', 0); % (q1, n2) x (q2, n3)
    coeffs = cellfun(@(x) vec(I1(row,col,:))*vec(x).', coeffs, 'UniformOutput', 0); % n1 x (q1, n2, q2, n3)
    coeffs = cellfun(vec, coeffs, 'UniformOutput', 0);
    
    T = T + vertcat(coeffs{:})./L^2;
end

T = accumarray(blk_id, T, [], @(x) {x});
T = cellfun(@(x,y) reshape(x, n1_len, y, n2_len, y, n3_len), T, q_list, 'UniformOutput', 0); 
T = cellfun(@(x) permute(x, [1,3,5,2,4]), T, 'UniformOutput', 0); % n1 x n2 x n3 x q1 x q2