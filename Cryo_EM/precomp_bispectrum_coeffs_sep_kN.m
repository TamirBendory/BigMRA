function B = precomp_bispectrum_coeffs_sep_kN(maxL, maxS, L)

Lr = floor(L/2);
n_list = precomp_Nn_list(Lr);
beta = sph_Bessel_to_2D_PSWF_factors(maxL, n_list, maxS, Lr); % l x N x s x n

[x,y] = meshgrid(-Lr:Lr, -Lr:Lr);
[phi, r] = cart2pol(x,y);
r = r./Lr;
beta = 1; 
T = 1e-1;
realFlag = 1;
[PSWF_Nn_p, PSWF_quad_int] = precomp_pswf_t(L, beta, T, realFlag);

B = cell(maxL+1, 1); % k x L2 x L3
for k = 0:maxL
    L23_vals = k:maxL;
    B{k+1} = cell(length(L23_vals), length(L23_vals));
    for ii_2 = 1:length(L23_vals)
        L2 = L23_vals(ii_2);
        for ii_3 = 1:length(L23_vals)
            L3 = L23_vals(ii_3);
            
            L1_vals = max(abs(L2-L3),k):min(L2+L3, maxL);
            B{k+1}{ii_2, ii_3} = cell(length(L1_vals), 1);
            for jj = 1:length(L1_vals)
                L1 = L1_vals(jj);
                tmp_acc = 0;
                parfor NN = 1:(2*L2+1)*(2*L3+1)
                    [N2, N3] = ind2sub([2*L2+1, 2*L3+1], NN);
                    N2 = N2-1; N3 = N3-1;
                    N1 = N2 + N3;
                    if abs(N1) > L1, continue; end
                    
                    T = generate_radial_triple_prods_given_N...
                        (Lr, N1, N2, N3, n_list, r, phi, PSWF_Nn_p, PSWF_quad_int);
                    
                    % Product and sum over n1:
                    tmp = conj(beta{L1+1}{1})*reshape(T, n_list(1), []); % s1 x (n2, n3, q1, q2)
                    tmp = reshape(tmp.', n_list(k+1), []); % n2 x (n3, q1, q2, s1)
                    
                    % Product and sum over n2:
                    tmp = beta{L2+1}{k+1}*tmp; % s2 x (n3, q1, q2, s1)
                    tmp = reshape(tmp.', n_list(k+1), []); % n3 x (q1, q2, s1, s2)
                    
                    % Product and sum over n3:
                    tmp = beta{L3+1}{k+1}*tmp; % s3 x (q1, q2, s1, s2)
                    tmp_acc = tmp_acc + ...
                        (-1)^(N1)*double(w3j(L2, L3, L1, N2, N3, N1))*...
                        reshape(tmp.', q_list(k+1), q_list(k+1), size(beta{L1+1}{1},1), size(beta{L2+1}{k+1},1), size(beta{L3+1}{k+1},1));
                    % q1 x q2 x s1 x s2 x s3
                end
                B{k+1}{ii_2, ii_3}{jj} = tmp_acc;
            end
        end
    end
end