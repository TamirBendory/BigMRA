function [m2_rec, m3_rec] = recon_moments_custom_PSWFs(m2, m3, L, beta_PSWF, T)

maxN = length(m3)-1;
psi_Nn = PSWF_2D_full_cart(maxN, L, beta_PSWF, T);
psi_Nn = psi_Nn(maxN+1:end);

m2_rec = zeros(2*L+1, 2*L+1);
m2_rec(:) = real(psi_Nn{1}*m2);

m3_rec = zeros((2*L+1)^2);
for ii = 1:length(psi_Nn)
    PSWF_funcs_curr = psi_Nn{ii};
    if ii == 1
        m3_rec = m3_rec + real(PSWF_funcs_curr*m3{ii}*PSWF_funcs_curr');
    else
        m3_rec = m3_rec + 2*real(PSWF_funcs_curr*m3{ii}*PSWF_funcs_curr');
    end
end

m3_rec = reshape(m3_rec, 2*L+1, 2*L+1, 2*L+1, 2*L+1);