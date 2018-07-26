function T = generate_radial_triple_prods(c, maxN, maxK, n_list, q_list)

% Generate grid for r integral:
[r_quad,w_r] = lgwt(1e2, 0, 1);
w_r = w_r.*r_quad;
phi_quad = linspace(-pi, pi, 1e2); 
w_phi = 2*pi/length(phi_quad)*ones(length(phi_quad),1);
[R,P] = ndgrid(r_quad, phi_quad);
w_q = kron(w_phi, w_r);
[X_q, Y_q] = pol2cart(P(:), R(:));

% Generate grid for dr integral:
[xx, w_x] = lgwt(1e2, -1,1);
[X,Y] = ndgrid(xx, xx);
[P_dr, R_dr] = cart2pol(X(:), Y(:));
w_dr = kron(w_x, w_x);

X_sum = bsxfun(@plus, X_q(:), X(:).'); % (x) x (dx)
Y_sum = bsxfun(@plus, Y_q(:), Y(:).'); % (y) x (dy)
[P_sum, R_sum] = cart2pol(X_sum, Y_sum);

rho_kN = cell(maxK+1, 2*maxN+1);
parfor ii = 1:(2*maxN+1)*(maxK+1)
    [k, N] = ind2sub([maxK+1, 2*maxN+1], ii);
    k = k-1; N = N-maxN-1;
    
    psi_Nn = PSWF_2D_N(N, n_list(abs(N)+1), c, R_sum(:), P_sum(:));
    [psi_Kq, alpha_Kq] = PSWF_2D_N(k, q_list(k+1), 2*c, R_dr(:), P_dr(:)); % dr x q
    psi_Kq = bsxfun(@times, psi_Kq, (2./alpha_Kq.').^2);
    
    psi_Nn = reshape(psi_Nn, numel(X_q), numel(X), []); % r x dr x n
    psi_Nn = permute(psi_Nn, [2,1,3]); % dr x r x n
    psi_Nn = reshape(psi_Nn, numel(X), []); % dr x (r, n)
    
    % dr integral:
    rho_kN{ii} = psi_Nn.'*diag(w_dr)*psi_Kq; % (r, n) x q
    rho_kN{ii} = reshape(rho_kN{ii}, numel(X_q), []); % r x (n, q)
end

T = cell(maxK+1, 2*maxN+1, 2*maxN+1); % k x N2 x N3
parfor ii = 1:(maxK+1)*(maxN+1)^2
    [k, N2, N3] = ind2sub([maxK+1, 2*maxN+1, 2*maxN+1], ii);
    k = k-1; N2 = N2-maxN-1; N3 = N3-maxN-1;
    N1 = N2+N3;
    if abs(N1) > maxN, continue; end
    psi_Nn = PSWF_2D_N(N1, n_list(abs(N1)+1), c, R(:), P(:)); % r x n1
    
    T{ii} = bsxfun(@times, rho_kN{k+1}{N2+maxN+1}, permute(rho_kN{k+1}{N3+maxN+1}, [1,3,2])); % r x (n2, q1) x (n3, q2)
    T{ii}= reshape(T{ii}, numel(X_q), []); % r x (n2, q1, n3, q2)
    T{ii} = psi_Nn.'*diag(w_q)*T{ii}; % n1 x (n2, q1, n3, q2)
    
    T{ii} = reshape(T{ii}, n_list(abs(N1)+1), n_list(abs(N2)+1), q_list(k+1), n_list(abs(N3)+1), q_list(k+1)); % n1 x n2 x q1 x n3 x q2
    T{ii} = reshape(T{ii}, n_list(abs(N1)+1), n_list(abs(N2)+1), q_list(k+1), n_list(abs(N3)+1), q_list(k+1)); % n1 x n2 x q1 x n3 x q2
    T{ii} = permute(T{ii}, [1,2,4,3,5]); % n1 x n2 x n3 x q1 x q2
end
    