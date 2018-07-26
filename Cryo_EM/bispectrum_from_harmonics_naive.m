function M = bispectrum_from_harmonics_naive(a_lms, L)

maxL = length(a_lms)-1;
maxn = size(a_lms{1},1);
maxI = maxn;
% beta = radial_PSWFs_3D_to_2D_factors(maxL, maxn, L); % l x N x s x n
beta = sph_Bessel_to_2D_PSWF_factors(maxL, maxn, maxI, L); % l x N x s x n
T = generate_radial_triple_prods(L, maxL, maxn); % k x n1 x n2 x n3 x q1 x q2)

M = cell(maxL+1, 1);
for k = 0:maxL
    M{k+1} = 0;
    for L2 = k:maxL
        for L3 = k:maxL
            L1_vals = abs(L2-L3):min(L2+L3, maxL);
            for jj = 1:length(L1_vals)
                L1 = L1_vals(jj);
                for m2 = -L2:L2
                    for m3 = -L3:L3
                        m1 = m2+m3;
                        if abs(m1) > L1, continue; end
                        for s1 = 1:maxI
                            for s2 = 1:maxI
                                for s3 = 1:maxI
                                    for n1 = 0:maxn
                                        for n2 = 0:maxn
                                            for n3 = 0:maxn
                                                M{k+1} = M{k+1} + conj(a_lms{L1+1}(s1, m1+L1+1))*a_lms{L2+1}(s2, m2+L2+1)*a_lms{L3+1}(s3, m3+L3+1)*...
                                                    sqrt(2*pi)*(-1)^(m1)*double(w3j(L2, L3, L1, k, -k, 0))*double(w3j(L2, L3, L1, m2, m3, -m1))*...
                                                    conj(beta{L1+1}{1}(s1,n1+1))*beta{L2+1}{k+1}(s2,n2+1)*beta{L3+1}{k+1}(s3,n3+1)*...
                                                    squeeze(T{k+1}(n1+1,n2+1,n3+1,:,:));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end