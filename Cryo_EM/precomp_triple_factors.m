% Precomputes factors for the expression of the bispectrum in terms of
% spherical harmonic coefficients.

addpath('../SHT')

L = 7;
C = cell(L+1, L+1); % indexed ell_1, ell_2, ell_3
% C{L1, L2, L3} is (2L1 + 1) x (2L2 + 1) and indexed by m1, m2.
Y_N = getSH(L, [0, pi/2], 'complex'); % precompute sphericals

if isempty(gcp('nocreate')), parpool('local', maxNumCompThreads); end
for L1 = 0:L
    for L2 = 0:L
        C{L1+1, L2+1} = cell(min(L1 + L2, L) - abs(L1-L2) + 1, 1);
        for L3 = abs(L1-L2):min(L1+L2, L)
            disp(['L1 = ' num2str(L1) '; L2 = ' num2str(L2) '; L3 = ' num2str(L3)])
            temp = zeros(2*L1+1, 2*L2+1, L1+1, L2+1);
            list1 = -L1:2:L1; list2 = -L2:2:L2;
            parfor ii = 1:(2*L1+1)*(L1+1)*(2*L2+1)*(L2+1)
                [m1, m2, m1p, m2p] = ind2sub([2*L1+1, 2*L2+1, L1+1, L2+1], ii);
                m1 = m1 - L1-1; m2 = m2 - L2-1;
                m1p = list1(m1p); m2p = list2(m2p);
                if abs(m1+m2) <= L3 && abs(m1p + m2p) <= L3
                    temp(ii) = ...
                        (-1)^(m1+m2+m1p+m2p)*w3j(L1, L2, L3, m1, m2, -m1-m2)*w3j(L1, L2, L3, m1p, m2p, -m1p-m2p)...
                        *Y_N(m1p+L1^2+L1+1)*Y_N(m2p+L2^2+L2+1)*Y_N(m1p+m2p+L3^2+L3+1);
                end
            end
            C{L1+1, L2+1}{L3 - abs(L1-L2) + 1} = temp;
        end
    end
end

save('triple_factors.mat','C','L')
            