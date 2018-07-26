function [m2_rec, m3_rec] = recon_moments_PSWFs(m2, m3, PSWF_Nn_p, points_inside_the_circle)

L = PSWF_Nn_p.L;
PSWF_funcs = PSWF_Nn_p.samples;
ang_freqs = PSWF_Nn_p.ang_freq;
ang_freqs_unique = unique(ang_freqs);

m2_rec = zeros(2*L, 2*L);
m2_rec(points_inside_the_circle) = real(PSWF_funcs(:,ang_freqs==0)*m2);

m3_rec = zeros((2*L)^2);
for ii = 1:length(ang_freqs_unique)
    PSWF_funcs_curr = PSWF_funcs(:,ang_freqs==ang_freqs_unique(ii)); 
    if ii == 1
        m3_rec(points_inside_the_circle, points_inside_the_circle) = m3_rec(points_inside_the_circle, points_inside_the_circle)...
            + real(PSWF_funcs_curr*m3{ii}*PSWF_funcs_curr');
    else
        m3_rec(points_inside_the_circle, points_inside_the_circle) = m3_rec(points_inside_the_circle, points_inside_the_circle)... 
            + 2*real(PSWF_funcs_curr*m3{ii}*PSWF_funcs_curr');
    end
end

m3_rec = reshape(m3_rec, 2*L, 2*L, 2*L, 2*L);