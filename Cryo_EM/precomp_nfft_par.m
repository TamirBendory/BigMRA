function plan = precomp_nfft_par(PSWF_quad_int, L, beta)

%% Precompute NFFT plan:
c = beta*pi*L;
sizeX = 2*L;
quadRulePtsX = PSWF_quad_int.quadRulePtsX; quadRulePtsY = PSWF_quad_int.quadRulePtsY;

usFftPts = c/L*[quadRulePtsX.' quadRulePtsY.'];
M=size(usFftPts,1);
x = usFftPts.'/pi/2;

plan = nfft_init_guru(2, sizeX, sizeX, M, 2*sizeX, 2*sizeX, 6, bitor(PRE_PHI_HUT,PRE_PSI), []);
nfft_set_x(plan,x);
nfft_precompute_psi(plan);

end

