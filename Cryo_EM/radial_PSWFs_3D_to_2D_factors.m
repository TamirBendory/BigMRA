function beta = radial_PSWFs_3D_to_2D_factors(maxL, maxn, L)

c = pi*L;
D = 3;
minEigenvalRatio = 10^-40;
matdim = 1000;
prolate_crea_options.isfixfirst = 1;
[r,w]=lgwt(20*maxn,0,1);
w = w.*r;

Y_l = YN2YL(getSH(maxL, [0, pi/2], 'complex'));

beta = cell(maxL+1,1);
for l = 0:maxL
    beta{l+1} = cell(l+1, 1);
    for N = 0:l % only compute for positive N
        [R_Nn, alpha_Nn_2D] = PSWF_radial_2D(N, maxn, c, r); % generate 2D radial prolates
        
        % generate 3D radial prolates:
        prolate_dat = prolate_crea(c,D,l,minEigenvalRatio, matdim, prolate_crea_options);
        Phi_ls = prolate_ev(prolate_dat, 0:prolate_dat.num_prols-1 , r);
        
        alpha_ratio = repmat((1/4)*alpha_Nn_2D(:).', size(Phi_ls,2), 1); % matrix of 2D eigenvalues
        
        beta{l+1}{N+1} = sqrt(2*pi)*Y_l{l+1}(l+1+N) .* alpha_ratio .* (Phi_ls.'*diag(w)*R_Nn); % s x n
    end
end
    