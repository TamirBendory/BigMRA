function n_list = precomp_Nn_list(L, beta, T, realFlag)

PSWF_Nn_p = precomp_pswf_t(L, beta, T, realFlag);
ang_freq = PSWF_Nn_p.ang_freq;
n_list = zeros(max(ang_freq)+1, 1);
for N = 0:max(ang_freq)
    n_list(N+1) = length(find(ang_freq == N));
end
