function [psi_Nn, alpha_Nn] = PSWF_2D_N(N, n, c, r_vals, phi_vals)

if n == 0
    [R_Nn, alpha_Nn] = PSWF_radial_2D(abs(N), 1, c, r_vals);
    R_Nn = R_Nn(:,1); alpha_Nn = alpha_Nn(1);
else
    [R_Nn, alpha_Nn] = PSWF_radial_2D(abs(N), n, c, r_vals);
end
R_Nn(r_vals > 1,:) = 0;
if N < 0
    R_Nn = (-1)^N*R_Nn; % account for multiplication by eigenvalue
end
psi_Nn = bsxfun(@times, R_Nn, exp(1i*N*phi_vals)/sqrt(2*pi));

    
    