function [Psi_lms, pts_in_ball, Psi_lms_2D, pts_in_disc_2D] = generate_hybrid_PSWFs(L, maxL)

c=pi*L;   % band limit

R = floor(L/2);
if mod(L,2) == 0
    [x,y,z] = meshgrid(-R:R-1, -R:R-1, -R:R-1);
else
    [x,y,z] = meshgrid(-R:R, -R:R, -R:R);
end
[phi,elev,r] = cart2sph(x,y,z);
r = r/R;
theta = pi/2-elev;
pts_in_ball = find(r <= 1);
pts_in_disc_2D = find(r(z==0) <= 1);

Psi_lms = cell(maxL+1, 1);
Psi_lms_2D = cell(maxL+1, 1);
Y_l = YN2YL(getSH(maxL, [phi(pts_in_ball), theta(pts_in_ball)], 'complex'));
n = L+1; T = 1e-1;
for l = 0:maxL
    [Phi_ls, alpha] = PSWF_radial_2D(l, n, c, r(pts_in_ball));
    lambda = (c/(2*pi))^2 * abs(alpha).^2;
    gamma = sqrt(abs(lambda./(1-lambda))); 

    n_end = find(gamma<=T,1,'first');
    
    if (~isempty(n_end))
        n = n_end;
    else
        % The case where the initial n isn't large enough -> Try again.        
        n = 2*n;
    end
    Psi_lms{l+1} = bsxfun(@times, Phi_ls, permute(Y_l{l+1},[2,3,1]));
    Psi_lms_2D{l+1} = Psi_lms{l+1}(z(pts_in_ball) == 0, :, 1:2:end); % restrict to xy-plane and take only m = l (mod 2)
end

end